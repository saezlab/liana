tib2list <- function(op_resource){
    if(is_tibble(op_resource)){
        op_resource %<>% list("Placeholder" = op_resource[[1]])
    } else{
        op_resource
    }
}

list2tib <- function(res){
    if(length(res)==1){res %>% pluck(1)} else{res}
}
