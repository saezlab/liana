# SUPPORT FUNCTIONS
get_resource_list <- function(){
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
}

make_upset_list <- function(ligrec_list){
  lapply(ligrec_list, function(x){
    if(typeof(x) == "list"){
      unlist(dplyr::select(x,rowid), use.names=FALSE)
    }
    else{
      return(integer(0))
    }
  })
}

