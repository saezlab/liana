# SUPPORT FUNCTIONS
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
