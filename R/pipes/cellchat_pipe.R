#' Connectome x Omni Pipe function
#' @param op_resource OmniPath Intercell Resource DN
#' @param seurat_object Seurat object as input
#' @param use_omni_resource whether to use Omni resources or run the analysis
#' with the default DB
#' @param exclude_anns Annotation criteria to be excluded
#' @inheritDotParams CellChat::subsetCommunication
#' @return A DF of intercellular communication network
call_cellchat <- function(op_resource,
                          seurat_object,
                          use_omni_resource = TRUE,
                          exclude_anns = c("ecm-receptor",
                                           "cell_surface_ligand-receptor"),
                          ...
                          ){
    library(CellChat)
    library(patchwork)
    library(ggalluvial)
    library(igraph)

    # initialize object
    cellchat.omni <- createCellChat(object = seurat_object,
                                    group.by = "ident")
    # future::plan("multiprocess", workers = 4) # do parallel

    # load CellChatDB
    CellChatDB.omni <- CellChatDB.human

    if(use_omni_resource){

        # get complexes and interactions from omnipath
        complex_interactions <- op_resource %>%
            select("ligand" = source_genesymbol,
                   "receptor" = target_genesymbol,
                   "evidence" = sources,
                   category_intercell_source,
                   category_intercell_target,
                   genesymbol_intercell_source,
                   genesymbol_intercell_target,
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

                suppressMessages(message(interaction_name))

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
            mutate(interaction_name_2 = ifelse(str_detect(.data$genesymbol_intercell_target, "^COMPLEX"),
                                               str_glue("{ligand} -", "{str_split(genesymbol_intercell_target, pattern='_')}"), interaction_name_2)) %>%
            mutate(interaction_name_2 = ifelse(str_detect(.data$genesymbol_intercell_source, "^COMPLEX"),
                                               str_glue("{ligand} -", "{str_split(genesymbol_intercell_source, pattern='_')}"), interaction_name_2)) %>%
            mutate(interaction_name_2 = str_replace(interaction_name_2, "c", "")) %>%
            mutate(interaction_name_2 = str_replace(interaction_name_2, "COMPLEX:", "")) %>%
            mutate(interaction_name_2 = str_replace(interaction_name_2, "-", "- ")) %>%
            mutate(interaction_name_2 = str_replace_all(interaction_name_2, ", ", "+")) %>%
            mutate(interaction_name_2 = str_replace_all(interaction_name_2, '"', ""))
        # select(all_of(names(interaction_input)))

        # Get Omni Complexes
        omni_complexes <- complex_interactions %>%
            filter(str_detect(genesymbol_intercell_source, "COMPLEX") |
                       str_detect(genesymbol_intercell_target, "COMPLEX")) %>%
            select(genesymbol_intercell_source,
                   genesymbol_intercell_target)

        # Convert to CellChat format
        omni_complexes <- union(omni_complexes$genesymbol_intercell_source,
                                omni_complexes$genesymbol_intercell_target) %>%
            str_subset(pattern = "COMPLEX") %>%
            str_replace(pattern = "COMPLEX:", "") %>%
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

    message(str_glue("Number of interactions:
                     {length(CellChatDB.omni$interaction$interaction_name)}"))

    ## set the used database in the object
    cellchat.omni@DB <- CellChatDB.omni

    ## subset the expression data of signaling genes
    cellchat.omni <- subsetData(cellchat.omni)

    # Infer the cell state-specific communications,
    cellchat.omni <- identifyOverExpressedGenes(cellchat.omni)
    cellchat.omni <- identifyOverExpressedInteractions(cellchat.omni)

    ## Compute the communication probability and infer cellular communication network
    # NOTE !!!!! here we might be able to extend it further with OmniPath
    # However, this might be too time-consuming and irrelevant.
    # cellchat.omni <- projectData(cellchat.omni, PPI.human)
    cellchat.omni <- computeCommunProb(cellchat.omni,
                                       raw.use=TRUE)

    # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
    cellchat.omni <- filterCommunication(cellchat.omni,
                                         min.cells = 10)


    # Extract the inferred cellular communication network
    df.omni <- subsetCommunication(cellchat.omni, ...)
    return(df.omni)
}
