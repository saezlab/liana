# Handle method choices
run_liana <- function(seurat_object, method, resource){
    resource %<>% .select_resource

    .select_method(method) %>%
        map2(names(.),
             function(.method, method_name){

                 resource %<>% if_else(method_name %in% c("squidpy", "natmi"),
                                       .list2tib(resource))

                 args <- append(
                     list(seurat_object = seurat_object,
                          op_resource = resource),
                     options(str_glue('{method_name}.defaults'))[[1]]
                     )

            exec(.method,  !!!args)
        })
    # Handle Exceptions, if possible with catch
}

run_liana(seurat_object, 'italk', 'Default')


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

# call cellchat
args <- append(
    list(seurat_object = seurat_object,
         op_resource = NULL),
    options('cellchat.defaults')[[1]]
    )
.select_method(c('cellchat')) %>%
    map(function(method) exec(method, !!!args))


# Handle default options
# CellChat
options('cellchat.defaults')[[1]] %>%
    `%||%`(list(
        nboot = 100,
        exclude_anns = NULL,
        thresh = 1,
        assay = "RNA",
        .normalize = TRUE,
        .do_parallel = FALSE,
        .raw_use = TRUE
    )) %T>%
    options(cellchat.defaults = .)
options('cellchat.defaults')[[1]]

options('connectome.defaults')[[1]] %>%
    `%||%`(list(
        .spatial = FALSE,
        min.cells.per.ident = 1,
        p.values = TRUE,
        calculate.DOR = FALSE,
        assay = 'RNA',
        .format = TRUE
    )) %T>%
    options(connectome.defaults = .)
options('connectome.defaults')[[1]]


options('italk.defaults')[[1]] %>%
    `%||%`(list(
        assay = 'RNA',
        .format = TRUE,
        .DE = TRUE
    )) %T>%
    options(italk.defaults = .)
options('italk.defaults')[[1]]

options('natmi.defaults')[[1]] %>%
    `%||%`(list(
        omnidbs_path = "input/omnipath_NATMI",
        natmi_path = "NATMI/",
        em_path = "input/test_em.csv",
        ann_path = "input/test_metadata.csv",
        output_path = "output/NATMI_test",
        .write_data = TRUE,
        assay = "RNA"
    )) %T>%
    options(natmi.defaults = .)
options('natmi.defaults')[[1]]


options('sca.defaults')[[1]] %>%
    `%||%`(list(
        assay = 'RNA',
        .format = TRUE,
        s.score = 0,
        logFC = log2(1.5)
    )) %T>%
    options(sca.defaults = .)
options('sca.defaults')[[1]]


options('squidpy.defaults')[[1]] %>%
    `%||%`(list(
        python_path = "/home/dbdimitrov/anaconda3/envs/theisverse/bin/python",
        cluster_key="seurat_annotations",
        n_perms=1000,
        threshold=0.01,
        seed=as.integer(1004)
    )) %T>%
    options(squidpy.defaults = .)
options('squidpy.defaults')[[1]]



# Handle pre-defined options
.tib2list <- function(op_resource){
    if(is_tibble(op_resource)){
        op_resource %<>% list("Placeholder" = op_resource[[1]])
    } else{
        op_resource
    }
}

# Handle list or not list method calls
.list2tib <- function(res){
    if(length(res)==1){res %>% pluck(1)} else{res}
}
