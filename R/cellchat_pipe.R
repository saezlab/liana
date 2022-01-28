#' Run CellChat with OmniPath function [[DEPRECATED]]
#'
#' @param op_resource OmniPath Intercell Resource DN
#' @param sce Seurat object as input
#' @param exclude_anns Annotation criteria to be excluded
#' @param nboot number of bootstraps to calculate p-value
#' @param .format bool whether to format output
#' @param .normalize # bool whether to normalize non-normalized data with
#' @param .raw_use whether use the raw data or gene expression data projectected
#'    to a ppi (should be kept to TRUE)
#' @param expr_prop minimum proportion of gene expression per cell type (0 by default),
#'  yet perhaps one should consider setting this to an appropriate value between 0 and 1,
#'  as an assumptions of these method is that communication is coordinated at the cluster level.
#' @param organism Obtain CellChatDB for which organism ('mouse' or 'human')
#' @param de_thresh diff expression of genes p-value
#'
#' @inheritDotParams CellChat::subsetCommunication
#'
#' @return A DF of intercellular communication network
#'
#' @importFrom Seurat Idents GetAssayData
#' @importFrom purrr pmap
#' @importFrom magrittr %>% %<>%
#' @importFrom dplyr select mutate mutate_at distinct_at filter
#' @importFrom tibble column_to_rownames enframe tibble
#' @importFrom tidyr unite unnest separate
#'
#' @export
#'
#' @details CellChat's objects are not lazily documented/exported thus the
#'   whole package has to be imported.
call_cellchat <- function(sce,
                          op_resource,
                          .format = TRUE,
                          exclude_anns = c(),
                          nboot = 100,
                          assay = "RNA",
                          .seed = 1004,
                          .normalize = FALSE,
                          .do_parallel = FALSE,
                          .raw_use = TRUE,
                          expr_prop = 0,
                          organism = "human",
                          thresh = 1,
                          de_thresh = 0.05,
                          ...
                          ){
    stringsAsFactors <- options('stringsAsFactors')[[1]]
    options(stringsAsFactors = FALSE)

    # Convert sce to seurat
    if(class(sce) == "SingleCellExperiment"){
        sce %<>% .liana_convert(., assay=assay)
    }

    # create a dataframe of the cell labels
    labels <- Seurat::Idents(sce)
    meta <- data.frame(group = labels, row.names = names(labels))

    cellchat.omni <- CellChat::createCellChat(
        object = `if`(!.normalize,
                      GetAssayData(sce,
                                   assay = assay,
                                   slot = "data"), # works with lognorm data
                      CellChat::normalizeData(
                          GetAssayData(sce,
                                       assay = assay,
                                       slot = "data"))
                      ),
        meta = meta,
        group.by = "group")
    cellchat.omni <- CellChat::addMeta(cellchat.omni,
                                       meta = meta)
    cellchat.omni <- CellChat::setIdent(cellchat.omni,
                                        ident.use = "group")

    if(.do_parallel){
        future::plan("multiprocess")
    }

    # load CellChatDB
    if(organism == "human"){
        ccDB <- CellChat::CellChatDB.human
    } else if(organism == "mouse") {
        ccDB <- CellChat::CellChatDB.mouse
    }

    if(!is.null(op_resource)){ # OmniPath resource conversion
        ccDB <- cellchat_formatDB(ccDB,
                                  op_resource,
                                  exclude_anns)

    } else{ # default CellChatDB
        ccDB$interaction <- ccDB$interaction %>%
            filter(!(annotation %in% exclude_anns))
    }

    ## set the used database in the object
    cellchat.omni@DB <- ccDB

    ## subset the expression data of signaling genes
    cellchat.omni <- CellChat::subsetData(cellchat.omni)

    # Infer the cell state-specific communications
    cellchat.omni <-
        CellChat::identifyOverExpressedGenes(cellchat.omni,
                                             thresh.pc = expr_prop,
                                             thresh.p = de_thresh)
    cellchat.omni <- CellChat::identifyOverExpressedInteractions(cellchat.omni)

    ## Compute the communication probability and infer cellular communication network
    if(!.raw_use){
        cellchat.omni <- CellChat::projectData(cellchat.omni,
                                               CellChat::PPI.human)
    }

    cellchat.omni <- CellChat::computeCommunProb(cellchat.omni,
                                                 raw.use = .raw_use,
                                                 seed.use = .seed,
                                                 do.fast = TRUE,
                                                 nboot = nboot)

    # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
    cellchat.omni <- CellChat::filterCommunication(cellchat.omni, min.cells = 1)

    # Extract the inferred cellular communication network
    df.omni <- CellChat::subsetCommunication(cellchat.omni, thresh = thresh, ...)

    if(.format){
        df.omni <- df.omni %>%
            select(source,
                   target,
                   ligand,
                   receptor,
                   prob,
                   pval) %>%
            as_tibble()
    }

    options(stringsAsFactors = stringsAsFactors)

    return(df.omni)

}


#' Helper Function to Format CellChatDB
#'
#' @param ccDB Inbuilt cellchatDB object
#'
#' @inheritParams call_cellchat
#' @importFrom tibble enframe column_to_rownames
#' @importFrom magrittr %>%
#' @importFrom stringr str_glue str_detect str_replace str_replace_all
#' @importFrom tidyr unite separate
#' @importFrom dplyr mutate mutate_at select filter
#' @importFrom tidyselect everything
#' @importFrom purrr pmap
#'
#' @export
cellchat_formatDB <- function(ccDB, op_resource, exclude_anns){
    # get complexes and interactions from omnipath
    complex_interactions <- op_resource %>%
        select(
            "ligand" = source_genesymbol,
            "receptor" = target_genesymbol,
            "evidence" = sources,
            category_intercell_source,
            category_intercell_target,
            is_directed,
            is_stimulation,
            is_inhibition
        ) %>%
        unite("annotation",
              c(category_intercell_source, category_intercell_target),
              sep="-") %>%
        unite("interaction_name", c(ligand, receptor), remove = FALSE) %>%
        mutate(pathway_name = "",
               agonist = "",
               antagonist = "",
               co_A_receptor = "",
               co_I_receptor = "") %>%
        mutate_at(vars(everything()), ~ replace(., is.na(.), ""))

    # Get OmniPath directed info
    omni_directions <- complex_interactions %>%
        select(interaction_name,
               is_directed,
               is_stimulation,
               is_inhibition,
               co_A_receptor,
               co_I_receptor
        ) %>%
        mutate(direction = pmap(., function(interaction_name,
                                            is_directed,
                                            is_stimulation,
                                            is_inhibition,
                                            co_A_receptor,
                                            co_I_receptor
        ){

            if(co_A_receptor=="" & co_I_receptor==""){

                if(is_directed==1){
                    if(is_stimulation==1 & is_inhibition==1){
                        return("Both")
                    } else if(is_stimulation==1){
                        return("Stimulation")
                    } else if(is_inhibition==1){
                        return("Inhibition")
                    }
                }
            } else{
                return(NA)
            }
        }
        )) %>% unnest(direction)

    # get complex type of interaction from OmniPath
    omni_interactions <- complex_interactions %>%
        # set LR interactions as rowname
        left_join(., (omni_directions %>%
                          select(interaction_name,
                                 direction)),
                  by = "interaction_name") %>%
        mutate_at(vars(everything()), ~ replace(., is.na(.), "")) %>%
        mutate(co_A_receptor = ifelse(.data$co_A_receptor == "" & (direction == "Stimulation") | (direction == "both"),
                                      "Stimulation", .data$co_A_receptor),
               co_I_receptor = ifelse(.data$co_I_receptor == "" & (direction == "Inhibition") | (direction == "both"),
                                      "Inhibition", .data$co_I_receptor)) %>%
        # remove duplicates and assign to colnames
        mutate("interaction_name2" = interaction_name) %>%
        distinct_at(.vars="interaction_name2", .keep_all = TRUE) %>%
        column_to_rownames("interaction_name2") %>%
        # NOTE: ligand - (subunit_1 + subunit_2)
        mutate(interaction_name_2 = str_glue("{ligand} - {receptor}")) %>%
        mutate(interaction_name_2 = ifelse(str_detect(.data$receptor, "^COMPLEX"),
                                           str_glue("{ligand} -", "{str_split(receptor, pattern='_')}"), interaction_name_2)) %>%
        mutate(interaction_name_2 = ifelse(str_detect(.data$ligand, "^COMPLEX"),
                                           str_glue("{ligand} -", "{str_split(ligand, pattern='_')}"), interaction_name_2)) %>%
        mutate(interaction_name_2 = str_replace(interaction_name_2, "c", "")) %>%
        mutate(interaction_name_2 = str_replace(interaction_name_2, "COMPLEX:", "")) %>%
        mutate(interaction_name_2 = str_replace(interaction_name_2, "-", "- ")) %>%
        mutate(interaction_name_2 = str_replace_all(interaction_name_2, ", ", "+")) %>%
        mutate(interaction_name_2 = str_replace_all(interaction_name_2, '"', ""))

    # Get Omni Complexes
    omni_complexes <- complex_interactions %>%
        filter(str_detect(ligand, "_") |
                   str_detect(receptor, "_")) %>%
        select(ligand, receptor)

    # Convert to CellChat format
    omni_complexes <- union(omni_complexes$ligand,
                            omni_complexes$receptor) %>%
        enframe() %>%
        separate(col=value, sep="_",
                 into = c("subunit_1", "subunit_2",
                          "subunit_3", "subunit_4", "subunit_5"), remove=FALSE) %>%
        mutate_at(vars(everything()), ~ replace(., is.na(.), "")) %>%
        select(-name) %>%
        column_to_rownames("value")

    # Replace Default DB with OmniPath Resource
    ccDB$interaction <- omni_interactions %>%
        filter(!(annotation %in% exclude_anns))

    ccDB$complex <- omni_complexes

    return(ccDB)
}
