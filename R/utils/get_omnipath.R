#' Function to call OmniPath resources
#' @return A list with intercellular resources from OmniPath
get_omni_resources <- function(){

    require(tidyverse)
    require(OmnipathR)

    omni_list <- list(
        'CellChatDB',
        'CellPhoneDB',
        'Ramilowski2015',
        'Baccin2019',
        'LRdb',
        'Kirouac2010',
        'ICELLNET',
        'iTALK',
        'EMBRACE',
        'HPMR',
        'Guide2Pharma',
        'connectomeDB2020',
        'talklr',
        'CellTalkDB',
        'OmniPath'
    )


    # Need to come back and expand on categories
    # Replace ligand and receptor with more descriptive ones when available?
    # Is cell-cell contact taken into account?
    # adhesion, ECM, etc
    setClass("OmniCriteria",
             slots=list(interactions_param="list",
                        transmitter_param="list",
                        receiver_param="list"))


    # Get a list with dataframes of omnipath resources
    omni_resources <- omni_list %>%
        map(function(x){
            message(x)

            if(x!="OmniPath"){
                x_obj = methods::new("OmniCriteria",
                                     interactions_param=list(resource=x),
                                     transmitter_param=list(resource=x,
                                                            category="ligand"),
                                     receiver_param=list(resource=x,
                                                         category="receptor"))

                import_intercell_network(
                    interactions_param = x_obj@interactions_param,
                    transmitter_param = x_obj@transmitter_param,
                    receiver_param = x_obj@receiver_param)
            } else{
                import_intercell_network(transmitter_param=list(category="ligand"),
                                         receiver_param=list(category="receptor"))
            }
        }) %>%
        setNames(omni_list)

    return(omni_resources)
}

