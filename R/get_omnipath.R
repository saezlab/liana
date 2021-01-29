omni_list <- list(
    # 'CellChatDB',
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


setClass("OmniCriteria",
         slots=list(interactions_param="list",
                    transmitter_param="list",
                    receiver_param="list"))


# Get a list with dataframes of omnipath resources
omni_resources <- omni_list %>%
    map(function(x){
        print(x)

        if(x!="OmniPath"){
            x_obj = methods::new("OmniCriteria",
                                 interactions_param=list(resource=x),
                                 transmitter_param=list(resource=x),
                                 receiver_param=list(resource=x))

            import_intercell_network(
                interactions_param = x_obj@interactions_param,
                transmitter_param = x_obj@transmitter_param,
                receiver_param = x_obj@receiver_param)
        } else{
            import_intercell_network()
        }
    }) %>%
    setNames(omni_list)
