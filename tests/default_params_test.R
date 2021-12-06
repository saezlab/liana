# Here we test the way that LIANA should be called for the manuscript
liana_path <- system.file(package = "liana")
seurat_object <- readRDS(file.path(liana_path , "testdata",
                                   "input", "testdata.rds"))


# Function to check the passed default params
liana_def_test <- function(seurat_object,
                           method,
                           resource = c('OmniPath'),
                           external_resource,
                           .simplify = TRUE,
                           ...){

    if(resource!='custom' & length(setdiff(resource, c(show_resources(), "all"))) > 0){
        stop(str_glue("{setdiff(resource, show_resources())} not part of LIANA "))
    }

    if(resource!='custom'){
        resource %<>% select_resource
    } else{
        resource = list('custom_resource'=external_resource)
    }

    ldefs <- liana_defaults(...)

    return(ldefs)
}


# Define and check list
def_list <- liana_def_test(seurat_object,
                           expr_prop=0,
                           squidpy.params=list(threshold = 0.1),
                           cellchat.params=list(nboot=1000),
                           method = c('call_natmi', 'call_connectome', 'logfc',
                                      'cellchat', 'call_sca', 'squidpy'))


#
liana_path <- system.file(package = "liana")
testdata <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))

liana_res <- liana_wrap(testdata,
                        method = c(
                            'call_natmi',
                            'call_connectome',
                            'logfc',
                            'cellchat',
                            'call_sca',
                            'cellphonedb',
                            "cytotalk"
                            ),
                        expr_prop=0.1,
                        cellchat.params = list(nboot=1000,
                                               expr_prop = 0))

saveRDS(liana_res, "data/output/test123.RDS")


