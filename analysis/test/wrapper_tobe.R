# Run LIANA Wrapper
#' @param seurat_object seurat object
#' @param method method(s) to be run via liana
#' @param resource resource(s) to be used by the methods
run_liana <- function(seurat_object, method, resource){
    resource %<>% .select_resource

    .select_method(method) %>%
        map2(names(.),
             function(.method, method_name){
                 if(!(method_name %in% c('squidpy', 'natmi'))){
                     map(resource, function(reso){
                         args <- append(
                             list(seurat_object = seurat_object,
                                  op_resource = reso),
                             options(str_glue('{method_name}.defaults'))[[1]]
                         )
                         exec(.method,  !!!args)
                     })
                 } else{
                 args <- append(
                     list(seurat_object = seurat_object,
                          op_resource = resource),
                     options(str_glue('{method_name}.defaults'))[[1]]
                 )
                 exec(.method,  !!!args)
                 }
        })
}


# Function to Change Assay for all methods

# Handle resource choices
.select_resource <- function(resource){
    omni_resources <- readRDS("input/omni_resources.rds")

    if(tolower(resource)=="all"){
        omni_resources[as.character(get_lr_resources())]
    } else{
        omni_resources[resource]
    }
}


# Select methods to run
# Adapted from decoupleR\https://github.com/saezlab/decoupleR/ (@jvelezmagic)
.select_method <- function(method){
    available_method <-
        list(
            cellchat = expr(call_cellchat),
            connectome = expr(call_connectome),
            italk = expr(call_italk),
            natmi = expr(call_natmi),
            sca = expr(call_sca),
            squidpy = expr(call_squidpyR)
        )

   method %>%
       tolower() %>%
       match.arg(names(available_method), several.ok = TRUE) %>%
       available_method[.]
}


# Helper Function to  Handle list or not list method calls
.list2tib <- function(res){
    if(length(res)==1){res %>% pluck(1)} else{res}
}





# Defaults
#' Default Arguments Helper Function
#' @param cellchat.defaults CellChat Default Parametres \link liana::call_cellchat
#' @details etc
liana_defaults <- function(
    cellchat.defaults = NULL,
    connectome.defaults = NULL,
    italk.defaults = NULL,
    natmi.defaults = NULL,
    sca.defaults = NULL,
    squidpy.defaults = NULL,
    ...){

    mod_args <- list(...)

    default_args <- list(
        "cellchat" = cellchat.defaults %<>%
            `%||%`(list(
                nboot = 100,
                exclude_anns = NULL,
                thresh = 1,
                assay = "RNA",
                .normalize = TRUE,
                .do_parallel = FALSE,
                .raw_use = TRUE
            )),

        'connectome' = connectome.defaults %<>%
            `%||%`(list(
                .spatial = FALSE,
                min.cells.per.ident = 1,
                p.values = TRUE,
                calculate.DOR = FALSE,
                assay = 'RNA',
                .format = TRUE
            )),

        'italk' = italk.defaults %<>%
            `%||%`(list(
                assay = 'RNA',
                .format = TRUE,
                .DE = TRUE
            )),

        'natmi' = natmi.defaults %<>%
            `%||%`(list(
                omnidbs_path = "data/input/omnipath_NATMI",
                natmi_path = "NATMI/",
                em_path = "data/input/test_em.csv",
                ann_path = "data/input/test_metadata.csv",
                output_path = "data/output/NATMI_test",
                assay = "RNA",
                .format = TRUE,
                .write_data = TRUE,
                .seed = 1004,
                .num_cor = 4)),

        'sca' = sca.defaults %<>%
            `%||%`(list(
                assay = 'RNA',
                .format = TRUE,
                s.score = 0,
                logFC = log2(1.5))),

        'squidpy' = squidpy.defaults %<>%
            `%||%`(list(
                python_path = "/home/dbdimitrov/anaconda3/envs/theisverse/bin/python",
                cluster_key="seurat_annotations",
                n_perms=1000,
                threshold=0.01,
                seed=as.integer(1004)
            ))
    )

    default_args %<>%
        map2(names(.), function(method, method_name){
            default_args[[method_name]] %>%
                map2(names(.), function(arg, arg_name){
                    if(!is.null(mod_args[[method_name]][[arg_name]])){
                        mod_args[[method_name]][[arg_name]]
                    } else{
                        arg
                    }
                })
        })

    default_args
}


liana_defaults()

default_args <- liana_defaults()

default_args$squidpy$python_path

liana_defaults(cellchat.defaults = NULL,
               connectome.defaults = NULL,
               italk.defaults = NULL,
               natmi.defaults = NULL,
               sca.defaults = NULL,
               squidpy.defaults = NULL,
               "squidpy" = list("python_path" = ":)",
                                "xd" = "):"))

mod_args <- list("squidpy" = list("python_path" = ":)",
                                  "xd" = "):"))

default_args %>%
    map2(names(.), function(method, method_name){
        mod_args[[method_name]] %>%
            map2(names(.), function(arg, arg_name){
                if((arg_name %in% names(mod_args[[method_name]]))){
                    arg
                } else{
                    .[[arg_name]]
                }
            })
    })
