#' Wrapper around `liana_wrap` to run liana for each sample.
#'
#' @param idents_col name of the cluster column
#'
#' @param sample_col name of the sample/context column
#'
#' @param verbose verbosity logical
#'
#' @param inplace logical (TRUE by default) if liana results are to be saved
#'  to the SingleCellExperiment object (`sce@metadata$liana_res`)
#'
#' @details takes a Seurat/SCE object and runs LIANA by sample/condition. The
#' key by which the samples are separated is build from the `condition_col` and
#' `sample_col`, separated by the `key_sep`.
#'
#' @inheritDotParams liana_wrap
#'
#' @export
#'
#' @returns If inplace is true returns a sce object with `liana_res` in
#' `sce@metadata`, else it returns a named list of tibble with liana results
#' per sample.
#'
liana_bysample <- function(sce,
                           idents_col,
                           sample_col,
                           verbose = TRUE,
                           inplace=TRUE,
                           ...){

    if(!is.factor(colData(sce)[[sample_col]])){
            liana_message(
                str_glue(
                    "`{sample_col}` was converted to a factor!"
                ), output = "message",
                verbose = verbose)
        }
    sce[[sample_col]] <- as.factor(sce[[sample_col]])


    # Map over key col
    liana_res <- map(levels(sce[[sample_col]]),
                     function(sample){

                         liana_message(str_glue("Current sample: {sample}"),
                                       output = "message",
                                       verbose = verbose)
                         # Subset to current sample
                         sce_temp <- sce[, sce[[sample_col]]==sample]

                         # Set cluster
                         colLabels(sce_temp) <- sce_temp[[idents_col]]

                         # Run LIANA on each
                         liana_wrap(sce=sce_temp, ...)

                     }) %>%
        setNames(levels(sce[[sample_col]]))

    if(!inplace){
        return(liana_res)
    } else{
        sce@metadata$liana_res <- liana_res
        return(sce)
    }
}

#' Join sample descriptor column from SCE to another table
#' @param sce SingleCellExperiment
#' @param right (df) to join to colData per sample
#' @param sample_col unique identifier column from colData
#' @param group_col sample_col descriptor column
#'
#' @noRd
.join_meta <- function(sce, right, sample_col, group_col){

    meta <- colData(sce) %>%
        as_tibble() %>%
        dplyr::select(!!sample_col, !!group_col) %>%
        distinct() %>%
        mutate( {{ group_col }} := as.factor(.data[[group_col]]))

    left_join(right, meta, by=sample_col)
}


#' Helper function to assign weights
#' @param lrs ligand_receptor tibble
#' @param resource resource
#' @param entity name of the entity
#'
#' @export
assign_lr_weights <- function(lrs,
                              resource,
                              entity="ligand"){

    entity_weight <- str_glue("{entity}_weight")
    entity_source <- str_glue("{entity}_source")

    entity_df <- lrs %>%
        pull(entity) %>%
        enframe(value="entity") %>%
        select(entity) %>%
        distinct() %>%
        liana::decomplexify("entity")

    entity_weights <- entity_df %>%
        left_join(resource, by = c("entity"="target")) %>%
        group_by(source, entity_complex) %>%
        # count number of subunits
        mutate(n_found = n()) %>%
        # count number of expected subunits x_y = 1 (+ 1)
        mutate(n_expected = str_count(entity_complex, "_") + 1) %>%
        filter(n_found==n_expected) %>%
        # only keep subunits which are sign-consistent
        summarise(weight = .sign_coh_mean(weight), .groups = "keep") %>%
        na.omit() %>%
        ungroup() %>%
        select({{ entity }} := entity_complex,
               {{ entity_source }} := source,
               {{ entity_weight }} := weight
        )


    return(entity_weights)
}


#' Filter nun-abundant cell types
#'
#' @inheritParams get_abundance_summary
#' @inheritParams plot_abundance_summary
#'
#' @export
#'
#' @return a filtered SingleCellExperiment Object
filter_nonabundant_celltypes <- function(sce,
                                         sample_col,
                                         idents_col,
                                         min_cells = 10,
                                         min_samples = 3,
                                         min_prop = 0.2,
                                         ctqc = NULL){

    # Calculate QC
    if(is.null(ctqc)){
        ctqc <- get_abundance_summary(sce,
                                      sample_col,
                                      idents_col,
                                      min_cells = min_cells,
                                      min_samples = min_samples,
                                      min_prop = min_prop
        )
    }

    plot_abundance_summary(ctqc)

    # Filter lowly-abundant celltypes from each sample
    keep_abundant <- ctqc %>%
        select(sample_col, idents_col, keep_min, keep_celltype) %>%
        filter(keep_min & keep_celltype) %>%
        ungroup()
    keep_abundant <- colData(sce) %>%
        as_tibble(rownames = "barcode") %>%
        select(barcode, sample_col = sample_col,idents_col = idents_col) %>%
        inner_join(keep_abundant, by = c("sample_col", "idents_col")) %>%
        pull(barcode)
    # Remove those cells which belong to cells not shared in enough samples
    sce <- sce[, keep_abundant]

    return(sce)
}


#' Function to get abundance summary
#'
#' @param sce SingleCellExperiment Object
#' @param sample_col column with sample ids
#' @param idents_col column with cell identity (cell type) ids
#' @param min_cells minimum cells per identity in each sample
#' @param min_samples minimum samples per cell identity
#' @param min_prop minimum proportion of samples in which a cell identity is
#' present (with at least `min_cells`)
#'
#' @returns a tibble
#'
#' @export
#'
#' @keywords internal
get_abundance_summary <- function(sce,
                                  sample_col,
                                  idents_col,
                                  min_cells = 10,
                                  min_samples = 3,
                                  min_prop = 0.2){
    colData(sce) %>%
        as_tibble() %>%
        dplyr::rename(sample_col = sample_col,
                      idents_col = idents_col) %>%
        # sample n
        mutate(sample_n = length(unique(sample_col))) %>%
        group_by(sample_col, idents_col) %>%
        # cells ct x sample
        mutate(cell_n = n()) %>%
        distinct_at(.vars = c("sample_col",
                              "sample_n",
                              "idents_col",
                              "cell_n")) %>%
        # keep cell type or not (per sample)
        mutate(keep_min = cell_n >= min_cells) %>%
        mutate(min_cells = min_cells) %>%
        ungroup() %>%
        # Keep celltype? (acc to minimum numbers across at least x samples
        # and at least x prop of the samples)
        group_by(idents_col) %>%
        mutate(keep_sum = sum(keep_min)) %>%
        mutate(sample_prop = keep_sum/sample_n) %>%
        mutate(keep_celltype = if_else((keep_sum >= min_samples) &
                                           (sample_prop >= min_prop),
                                       TRUE,
                                       FALSE)) %>%
        ungroup()
}

#' Function to Plot Abundance Summary
#'
#' @param ctqc cell type quality control summary obtained from
#' `get_abundance_summary`
#'
#' @param ncol number of columns for the facet wrap
#'
#' @return a ggplot2 object
#'
#' @export
#'
#' @keywords internal
plot_abundance_summary <- function(ctqc, ncol = 3){
    ctqc %>%
        ggplot(aes(x=sample_col, y=log10(cell_n), fill=keep_min)) +
        geom_bar(stat="identity", width=0.5) +
        scale_fill_manual(values = c('TRUE' = 'blue', 'FALSE' = 'red'),
                          guide="none") +
        geom_hline(yintercept = log10(unique(ctqc[["min_cells"]])),
                   linetype="dashed") +
        geom_text(
            mapping = aes(x = -Inf,
                          y = Inf,
                          label = str_glue("Prevalence: ",
                                           "{round(sample_prop, digits = 3)} ",
                                           "(N = {keep_sum})"),
                          color=keep_celltype),
            vjust=1.5, hjust=-0.1
        ) +
        scale_color_manual(values = c('TRUE' = 'black', 'FALSE' = 'red'),
                           guide="none") +
        facet_wrap(idents_col ~ ., ncol = ncol) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}
