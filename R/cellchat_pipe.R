#' Connectome x Omni Pipe function
#' @param op_resource OmniPath Intercell Resource DN
#' @param seurat_object Seurat object as input
#' @param exclude_anns Annotation criteria to be excluded
#' @param nboot number of bootstraps to calculate p-value
#' @param .format bool whether to format output
#' @param .normalize # bool whether to normalize non-normalized data with
#'  internal func
# #' @inheritDotParams CellChat::subsetCommunication
#'
#' @return A DF of intercellular communication network
#'
# #' @importFrom CellChat subsetCommunication createCellChat computeCommunProb
# #' @importFrom CellChat subsetData identifyOverExpressedGenes
# #' @importFrom CellChat identifyOverExpressedInteractions filterCommunication
# #' @importFrom Seurat Idents GetAssayData NormalizeData
#' @importFrom purrr pmap
#' @importFrom magrittr %>%
#' @importFrom dplyr select mutate mutate_at distinct_at filter
#' @importFrom tibble column_to_rownames enframe
#' @importFrom tidyr unite unnest separate
#' @importFrom stringr str_glue
#' @export
call_cellchat <- function(op_resource,
                          seurat_object,
                          .format = TRUE,
                          exclude_anns = c("ecm-receptor",
                                           "cell_surface_ligand-receptor"),
                          nboot = 100,
                          assay = "SCT",
                          .seed = 1004,
                          .normalize = FALSE,
                          .do_parallel = FALSE,
                          ...
                          ){

    stringsAsFactors <- options('stringsAsFactors')[[1]]

    options(stringsAsFactors = FALSE)

    data.input <- as.matrix(GetAssayData(seurat_object,
                                         assay = assay,
                                         slot = "data")) # data matrix
    if(.normalize){
        # function name NormalizeData is not capitalized?
        data.input <- normalizeData(data.input)
    }

    # create a dataframe of the cell labels
    labels <- Idents(seurat_object)

    meta <- data.frame(group = labels, row.names = names(labels))

    cellchat.omni <- createCellChat(object = data.input,
                               meta = meta,
                               group.by = "group")


    if(.do_parallel){
        future::plan("multiprocess") # do parallel
    }

    # load CellChatDB
    CellChatDB.omni <- CellChatDB.human

    if(!is.null(op_resource)){
        # get complexes and interactions from omnipath
        complex_interactions <- op_resource %>%
            select(
                "genesymbol_intercell_source" = target,
                "genesymbol_intercell_target" = source,
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

                log_trace(interaction_name)

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
                                     direction))) %>%
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
        # Here, I filter ECM, as when I don't an error is encountered
        # both with mine and their data set (they too filter for secreted signaling only)
        CellChatDB.omni$interaction <- omni_interactions %>%
            filter(!(annotation %in% exclude_anns))

        CellChatDB.omni$complex <- omni_complexes
    } else{
        CellChatDB.omni$interaction <- CellChatDB.omni$interaction %>%
            filter(!(annotation %in% exclude_anns))
    }

    log_info("Number of interactions:
                     {length(CellChatDB.omni$interaction$interaction_name)}")

    ## set the used database in the object
    cellchat.omni@DB <- CellChatDB.omni

    ## subset the expression data of signaling genes
    cellchat.omni <- subsetData(cellchat.omni)

    # Infer the cell state-specific communications,
    cellchat.omni <- identifyOverExpressedGenes(cellchat.omni)
    cellchat.omni <- identifyOverExpressedInteractions(cellchat.omni)

    ## Compute the communication probability and infer cellular communication network
    cellchat.omni <- projectData(cellchat.omni, PPI.human)
    cellchat.omni <- computeCommunProb(cellchat.omni,
                                       raw.use=FALSE,
                                       seed.use = .seed,
                                       nboot = nboot)

    # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
    cellchat.omni <- filterCommunication(cellchat.omni,
                                         min.cells = 10)


    # Extract the inferred cellular communication network
    df.omni <- subsetCommunication(cellchat.omni, ...)

    if(.format){
        df.omni <- df.omni %>%
            select(source,
                   target,
                   ligand,
                   receptor,
                   prob,
                   pval)
    }

    options(stringsAsFactors = stringsAsFactors)

    return(df.omni)

}