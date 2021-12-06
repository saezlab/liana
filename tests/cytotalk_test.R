#' Compute cross-talk score from a Seurat Object
#'
#' This function executes all the required functions to extract and compute the
#' cross-talk scores as defined by the Cytotalk authors. .
#'
#' @param seurat_object Seurat object as input
#' @param op_resouce The resource of ligand-receptor interactions
#'
#' @return A data frame containing the resulting scores
#'
#' @import SingleCellExperiment SeuratObject tibble tidyr dplyr
#'
call_cytotalk <- function(seurat_object,
                          op_resouce,
                          assay.type = "logcounts") {

    set.seed(1004)
    # Seurat to SCE
    ligand_receptor_df <- op_resource[, c("source_genesymbol", "target_genesymbol")]
    colnames(ligand_receptor_df) <- c("ligand", "receptor")

    entity_genes <- union(op_resource$source_genesymbol, op_resource$target_genesymbol)
    sce <- seurat_to_sce(seurat_object, entity_genes = entity_genes, assay = "RNA")

    # remove any genes that are not in the sce object
    ligand_receptor_df %<>%
        filter(ligand %in% rownames(sce)) %>%
        filter(receptor %in% rownames(sce)) %>%
        arrange(ligand, receptor)

    # get LR pairs
    pairs <- expand_grid(source = unique(colLabels(sce)),
                         target = unique(colLabels(sce)))

    # apply cytotalk functions
    pem_scores <- compute_pem_scores(sce = sce,
                                     assay.type = assay.type) %>%
        rownames_to_column("gene") %>%
        pivot_longer(-gene, names_to = "celltype", values_to = "pem")

    # Compute -log10 non-self-talk scores
    nst_scores <- compute_nst_scores(sce = sce,
                                     ligand_receptor_df = ligand_receptor_df,
                                     assay.type = assay.type,
                                     seed = 1004)

    # compute all possible pairs of cell types
    pairs <- combn(x = unique(colLabels(sce)), m = 2, simplify = FALSE)
    result_df <- lapply(pairs, function(cell_types) {

        out_df <- data.frame(source = cell_types, target = rev(cell_types)) %>%
            left_join(nst_scores, by = c("source" = "celltype")) %>%
            dplyr::rename(nst_score_ligand = nst_score) %>%
            left_join(nst_scores, by = c("target" = "celltype", "ligand", "receptor")) %>%
            dplyr::rename(nst_score_receptor = nst_score) %>%
            left_join(pem_scores, by = c("ligand" = "gene", "source" = "celltype")) %>%
            dplyr::rename(pem_ligand = pem) %>%
            left_join(pem_scores, by = c("receptor" = "gene", "target" = "celltype")) %>%
            dplyr::rename(pem_receptor = pem) %>%
            # compute cross-talk
            mutate(es = (pem_ligand + pem_receptor) / 2, # calculate Expression Score
                   nst = (nst_score_ligand + nst_score_receptor) / 2 # combine -log10 NST score
            ) %>%
            mutate(Nes = minmax(es), Nnst = minmax(nst)) %>%
            mutate(crosstalk_score = Nes * Nnst) %>%
            select(source, ligand, target, receptor, everything())

        return(out_df)

    }) %>%
        bind_rows()

    return(result_df)
}



# CytoTalk Test -----
### LIANA
liana_path <- system.file(package = "liana")
seurat_object <- readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))
assay.type = "logcounts"
op_resource <- select_resource("OmniPath")[[1]] %>% decomplexify()


# LIANA pipe out
liana_out <- liana_pipe(seurat_object = seurat_object,
                        op_resource = op_resource,
                        expr_prop = 0)

# Basic CCLR scaffold
liana_res <- liana_out %>%
    select(source, ligand, receptor, target,
           ligand.expr, receptor.expr,
           ligand.pem, receptor.pem)



liana_cytotalk <- liana_wrap(seurat_object,
                             method=c("cytotalk"),
                             expr_prop = 0,
                             resource = "OmniPath") %>%
    # we remove autocrine as original cytotalk does not return them
    # in liana we use the inverse of the NST for autocrine signalling
    # to make cytotalk comparable to the other methods
    filter(source!=target) %>%
    select(source, ligand, target, receptor, crosstalk_score) %>%
    arrange(crosstalk_score, source, ligand, target, receptor)


# Get cytotalk scores
martin_crosstalk <- call_cytotalk(seurat_object = seurat_object,
                                  op_resouce = op_resource) %>%
    # set crosstalk to 0 if either pem is 0
    mutate(crosstalk_score :=
               if_else(pem_ligand==0 | pem_receptor==0,
                       0,
                       crosstalk_score)) %>%
    select(source, ligand, target, receptor, crosstalk_score) %>%
    as_tibble() %>%
    arrange(crosstalk_score, source, ligand, target, receptor)

diff <- setdiff(liana_cytotalk, martin_crosstalk)
diff




### To integrate
# LIANA call natmi as example
liana_call("natmi",
           seurat_object = seurat_object,
           op_resource = op_resource,
           lr_res = liana_out,
           decomplexify = TRUE)

# call cytotalk
liana_call("cytotalk",
           seurat_object = seurat_object,
           op_resource = op_resource,
           lr_res = liana_out,
           sce = sce)

# NST - in LIANA wrap
nst_scores <- compute_nst_scores(sce = sce,
                                 ligand_receptor_df = lr_res %>%
                                     select(ligand,receptor) %>%
                                     distinct(),
                                 assay.type = assay.type,
                                 seed = 1234)

liana_cts <- lr_res %>%
    left_join(nst_scores, by = c("target" = "celltype", "ligand", "receptor")) %>%
    dplyr::rename(target.nst = nst_score) %>%
    left_join(nst_scores, by = c("source" = "celltype", "ligand", "receptor")) %>%
    dplyr::rename(source.nst = nst_score) %>%
    # CytoTalk Scores are computed by celltype pairs
    group_by(source, target) %>%
    # Expression and Non-self-talk scores
    mutate(es = (ligand.pem + receptor.pem)/ 2) %>% # calculate Expression Score
    mutate(nst = (source.nst + target.nst)/ 2) %>% # combine -log10 NST scores
    # normalize ES and NST
    mutate(Nes = minmax(es), Nnst = minmax(nst)) %>%
    # compute cross-talk scores
    mutate(crosstalk_score = if_else(source!=target,
                                     Nes * Nnst, # Paracrine as CytoTalk
                                     Nes * (1-Nnst) # Autocrine we inverse NST
                                     )) %>%
    ungroup()


# Test with a real dataset
seurat_object <- readRDS("~/Repos/ligrec_decouple/data/input/citeseq/5k_pbmcs/5k_pbmcs_seurat.RDS")
seurat_object@meta.data$seurat_clusters
test_wrap <- liana_wrap(seurat_object,
                        squidpy.params = list(cluster_key = "seurat_clusters"))

test_agg <- liana_aggregate(test_wrap)
