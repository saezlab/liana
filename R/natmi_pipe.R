#' Call NATMI Pipeline from R with Resources Querried from OmniPath
#' @param op_resource List of OmniPath resources
#' @param seurat_object Seurat object
#' @param omnidbs_dir path of saved omnipath resources
#' @param expr_file expression matrix file name
#' @param meta_file annotations (i.e. clusters) file name
#' @param output_dir NATMI output directory
#' @param .format bool whether to format output
#' @param .write_data bool whether Extract data from Seurat Object
#' @param .default_run bool whether to run default DBs or not
#' @param .natmi_path path of NATMI code and dbs (by default set to liana path)
#' @param assay Seurat assay to be used
#' @return DF with NATMI results
#'
#' @details
#'     This function will save NATMI dbs folder, then it will call the
#'     NATMI Python from the NATMI dir and save the output into a specified
#'     directory in NATMI's path.
#'     It will then load the csvs and format the output to a list of lists.
#'
#'     By default, NATMI's path is set to that of LIANA, but any alternative
#'     path can be passed
#'
#'==============================================================================
#' NATMI Arguments:
#'   --interDB INTERDB
#'                         lrc2p (default) has literature supported ligand-receptor pairs | lrc2a has putative and literature supported ligand-receptor pairs | the user-supplied interaction database can also be used by calling the name of database file without extension
#'   --interSpecies INTERSPECIES
#'                         human (default) | mouse | expandp | expanda
#'   --emFile EMFILE       the path to the normalised expression matrix file with row names (gene identifiers) and column names (cell-type/single-cell identifiers)
#'   --annFile ANNFILE     the path to the metafile in which column one has single-cell identifiers and column two has corresponding cluster IDs (see file 'toy.sc.ann.txt' as an example)
#'   --species SPECIES     human (default) | mouse | rat | zebrafish | fruitfly | chimpanzee | dog | monkey | cattle | chicken | frog | mosquito | nematode | thalecress | rice | riceblastfungus | bakeryeast | neurosporacrassa | fissionyeast | eremotheciumgossypii | kluyveromyceslactis, 21 species are supported
#'   --idType IDTYPE       symbol (default) | entrez(https://www.ncbi.nlm.nih.gov/gene) | ensembl(https://www.ensembl.org/) | uniprot(https://www.uniprot.org/) | hgnc(https://www.genenames.org/) | mgi(http://www.informatics.jax.org/mgihome/nomen/index.shtml) | custom(gene identifier used in the expression matrix)
#'   --coreNum CORENUM     the number of CPU cores used, default is one
#'   --out OUT             the path to save the analysis results
#' (Taken From NATMI's GitHub Page)
#'
#' Stats:
#' 1) The mean-expression edge weights
#' 2) The specificity-based edge weights
#' * a weight of 1 means both the ligand and receptor are only expressed
#'  in one cell type
#'
#' @importFrom reticulate py_set_seed
#' @importFrom stringr str_glue
#' @importFrom Seurat GetAssayData Idents
#' @import dplyr reticulate
#'
#' @export
call_natmi <- function(
    op_resource,
    seurat_object,
    expr_file = "em.csv",
    meta_file = "metadata.csv",
    output_dir = "NATMI_test",
    assay = "RNA",
    num_cor = 4,
    conda_env = NULL,
    .format = TRUE,
    .write_data = TRUE,
    .seed = 1004,
    .natmi_path = NULL){

    if(is_tibble(op_resource)){
        op_resource <- list("placeholder" = op_resource)
    }

    reticulate::use_condaenv(condaenv = conda_env %>% `%||%`("liana_env"),
                             conda = "auto",
                             required = TRUE)
    py$pd <- reticulate::import("pandas")
    python_path <- reticulate::py_discover_config()[[1]]
    reticulate::py_set_seed(.seed)

    .natmi_path %<>% `%||%`(system.file('NATMI/', package = 'liana'))
    .input_path = file.path(.natmi_path, 'data', 'input')
    .output_path = file.path(.natmi_path, 'data', 'output', output_dir)

    if(!dir.exists(file.path(.input_path))){
        log_success(str_glue("Input path created: {.input_path}"))
        dir.create(file.path(.input_path), recursive = TRUE)
    }

    if(!dir.exists(file.path(.output_path))){
        log_success(str_glue("Output path created: {.output_path}"))
        dir.create(file.path(.output_path), recursive = TRUE)
    }

    # append default resources to OmniPath ones
    if("DEFAULT" %in% toupper(names(op_resource))){
        op_resource %<>% purrr::list_modify("Default" = NULL)
        resource_names <- append(as.list(names(op_resource)), "lrc2p")
    } else{
        resource_names <- as.list(names(op_resource))
    }

    if(.write_data){
        log_info(str_glue("Writing EM to {.input_path}/{expr_file}"))
        write.csv(100 * (exp(as.matrix(GetAssayData(object = seurat_object,
                                                    assay = assay,
                                                    slot = "data"))) - 1),
                  file = file.path(.input_path, expr_file),
                  row.names = TRUE)
        log_info(str_glue("Writing Annotations to {.input_path}/{meta_file}"))
        write.csv(Idents(seurat_object)  %>%
                      enframe(name="barcode", value="annotation"),
                  file = file.path(.input_path, meta_file),
                  row.names = FALSE)

        # save OmniPath resources to NATMI dir
        log_info(str_glue("Saving resources to {.natmi_path}/lrdbs"))
        omni_to_NATMI(op_resource, file.path(.natmi_path, "lrdbs"))

    }

    # submit native sys request
    resource_names %>% map(function(resource){

        log_success(str_glue("Now Running: {resource}"))
            system(str_glue("{python_path} {.natmi_path}/ExtractEdges.py ",
                            "--species human ",
                            "--emFile {.input_path}/{expr_file} ",
                            "--annFile {.input_path}/{meta_file} ",
                            "--interDB {.natmi_path}/lrdbs/{resource}.csv ",
                            "--coreNum {num_cor} ",
                            "--out {.output_path}/{resource}",
                            sep = " "))
    })

    # load and format results
    natmi_results <- FormatNatmi(.output_path, resource_names, .format)

    return(natmi_results)
}



#' Reform OmniPath Resource to NATMI format and save to location
#' @param resource_names list of omnipath resources
#' @param omni_path directory in which to save OP resources
omni_to_NATMI <- function(op_resource,
                          omni_path = "input/omnipath_NATMI"){

    names(op_resource) %>%
        map(function(x){
            write.csv(op_resource[[x]]  %>%
                          select("Ligand gene symbol" = source_genesymbol,
                                 "Receptor gene symbol" = target_genesymbol) %>%
                          distinct() %>%
                          as.data.frame(),
                      file = str_glue("{omni_path}/{x}.csv"),
                      row.names = FALSE)
        })
}


#' Load NATMI results from folder and format appropriately
#'
#' @param output_path NATMI output path
#' @param resource_names results for which resources to load
#' @param .format bool whether to format output
#'
#' @return A list of NATMI results per resource loaded from the output
#'     directory.
#'
#' @importFrom tibble enframe deframe
#' @importFrom magrittr %>%
#' @importFrom purrr map
#' @importFrom dplyr mutate select if_else
#' @importFrom readr read_csv
#' @importFrom tidyr separate
FormatNatmi <- function(output_path,
                        resource_names,
                        .format = TRUE){

    list.files(output_path,
               all.files = TRUE,
               recursive = TRUE,
               pattern ="Edges_") %>%
        enframe() %>%
        separate(value,
                 into = c("resource", "file"),
                 remove = FALSE,
                 extra = "drop") %>%
        filter(resource %in% resource_names) %>%
        mutate(value =  value %>% map(function(csv)
            read.csv(str_glue("{output_path}/{csv}")))) %>%
        select(resource, "result" = value) %>%
        mutate(
            result =
                if_else(
                    rep(.format, length(.data$result)),
                    result %>% map(function(df){
                        df %>%
                            select(
                                source = Sending.cluster,
                                target = Target.cluster,
                                ligand = Ligand.symbol,
                                receptor = Receptor.symbol,
                                edge_avg_expr = Edge.average.expression.weight,
                                edge_specificity = Edge.average.expression.derived.specificity
                            ) %>%
                            as_tibble()
                    }),
                    result)) %>%
        deframe() %>%
        plyr::rename(., c("lrc2p" = "Default"), # change this to default
                     warn_missing = FALSE) %>%
        .list2tib()
}
