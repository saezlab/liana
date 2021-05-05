#' Shuffle OmniPath Intercell DB
#' @param op_resource Intrcell DB to shuffle
#' @param .seed Value for set.seed
#' @return A shuffled omnipath-formatted resource
#' @import BiRewire tibble
#' @export
shuffle_omnipath <- function(op_resource,
                             .seed = 1004){

    ### These packages could go to "Suggests" in DESCRIPTION
    ### because not all users want to install all the tools
    ### to run one of them. Functions from these packages
    ### should be referred by :: to avoid warnings
    # library(BiRewire)
    set.seed(.seed)

    # make a vector proportional to the number of consesus directions
    stimul_num <- round(mean(op_resource$consensus_stimulation)/mean(op_resource$consensus_direction) * 100000)
    directed_vector <- append(rep(1, stimul_num), rep(-1,100000 - stimul_num) )

    op_prep <- op_resource %>%
        # filter(entity_type_intercell_source != "complex",
        #        entity_type_intercell_target != "complex") %>%
        select(source_genesymbol, target_genesymbol,
               is_directed, is_stimulation, is_inhibition,
               consensus_direction, consensus_stimulation,
               consensus_inhibition) %>%
        mutate(mor = if_else(consensus_stimulation==1,
                             true=1,
                             if_else(consensus_inhibition==1, -1, 0))) %>%
        mutate(mor = if_else(mor==0, sample(directed_vector, 1), mor)) %>%
        select(source_genesymbol, mor, target_genesymbol) %>%
        distinct()

    # Induced bipartite and SIF builder
    op_dsg = birewire.induced.bipartite(op_prep,
                                        delimitators = list(negative = '-1',
                                                            positive = '1'))
    op_sif = birewire.build.dsg(op_dsg,
                                delimitators = list(negative = '-1',
                                                    positive = '1'))

    # Rewire dsg
    random_dsg = birewire.rewire.dsg(dsg = op_dsg)
    random_sif = birewire.build.dsg(random_dsg,
                                    delimitators = list(negative = '-1',
                                                        positive = '1'))
    # Jacard dsg
    message(str_glue("Jaccard index between random and original resource: ",
                     {birewire.similarity.dsg(op_dsg, random_dsg)}))

    # format to OmniPath
    op_random <- random_sif %>%
        select(
            source,
            target,
            sign) %>%
        mutate(
            source_genesymbol = source,
            target_genesymbol = target,
            is_directed = 1,
            is_stimulation = if_else(sign==1, 1, 0),
            consensus_stimulation = if_else(sign==1, 1, 0),
            is_inhibition = if_else(sign==-1, 1, 0),
            consensus_inhibition = if_else(sign==-1, 1, 0),
            category_intercell_source = "ligand",
            category_intercell_target = "receptor",
            genesymbol_intercell_source = source_genesymbol,
            genesymbol_intercell_target = target_genesymbol,
            entity_type_intercell_target = "protein",
            sources = "RANDOM",
            references = "BiRewire",
            entity_type_intercell_source = "protein",
            entity_type_intercell_target
            ) %>%
        as_tibble()


}
