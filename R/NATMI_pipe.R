#' Call NATMI Pipeline from R with OmniPath
#' @param op_resource List of OmniPath resources
#' @param omnidbs_path path of saved omnipath resources
#' @param natmi_path path of NATMI code and dbs
#' @param em_path expression matrix path
#' @param ann_path annotations (i.e. clusters) path
#' @param output_path NATMI output path
#' @param .format bool whether to format output
#' @param .write_data bool whether Extract data from Seurat Object
#' @param .default_run bool whether to run default DBs or not
#' @param assay Seurat assay to be used
#' @return DF with NATMI results
#'
#' @details
#'     This function will take omnipath resources saved to csvs and copy
#'     them to the NATMI dbs folder, then it will natively call the Python
#'     modules of NATMI in the NATMI dir and save the output into a specified
#'     directory. It will then load and format the output to a DF.
#'
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
#' @importFrom reticulate py_set_seed
#' @importFrom stringr str_glue
#' @importFrom Seurat GetAssayData Idents
#'
#' @export
call_natmi <- function(
    op_resource,
    seurat_object,
    omnidbs_path = "input/omnipath_NATMI",
    natmi_path = "NATMI/",
    em_path = "input/test_em.csv",
    ann_path = "input/test_metadata.csv",
    output_path = "output/NATMI_test",
    assay = "RNA",
    .format = TRUE,
    .write_data = TRUE,
    .seed = 1004,
    .num_cor = 4){

    py_set_seed(.seed)
    .natmi_dir <- system.file('NATMI/', package = 'liana')

    # append default resources to OmniPath ones
    if("DEFAULT" %in% toupper(names(op_resource))){
        op_resource %<>% purrr::list_modify("Default" = NULL)
        resource_names <- append(as.list(names(op_resource)), "lrc2p")
    } else{
        resource_names <- as.list(names(op_resource))
    }

    if(.write_data){
        log_info(str_glue("Writing EM to {em_path}"))
        write.csv(100 * (exp(as.matrix(GetAssayData(object = seurat_object,
                                                    assay = assay,
                                                    slot = "data"))) - 1),
                  file = em_path,
                  row.names = TRUE)
        log_info(str_glue("Writing Annotations to {ann_path}"))
        write.csv(Idents(seurat_object)  %>%
                      enframe(name="barcode", value="annotation"),
                  file = ann_path,
                  row.names = FALSE)

        log_info(str_glue("Saving resources to {omnidbs_path}"))
        omni_to_NATMI(op_resource, omnidbs_path)

        # copy OmniPath resources to NATMI dir
        file.copy(list.files(omnidbs_path, "*.csv$",
                             full.names = TRUE),
                  to=str_glue("{natmi_path}/lrdbs/"),
                  overwrite = TRUE)
    }

    log_success(str_glue("Output to be saved and read from {output_path}"))
    dir.create(file.path(output_path), recursive = TRUE)

    # submit native sys request
    resource_names %>% map(function(resource){

        log_success(str_glue("Now Running: {resource}"))
        system(str_glue("python3 {.natmi_dir}/ExtractEdges.py ",
                        "--species human ",
                        "--emFile {em_path} ",
                        "--annFile {ann_path} ",
                        "--interDB {.natmi_dir}/lrdbs/{resource}.csv ",
                        "--coreNum {.num_cor} ",
                        "--out {output_path}/{resource}",
                        sep = " "))
    })

    # load and format results
    natmi_results <- FormatNatmi(output_path, resource_names, .format)

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
#' @importFrom dplyr mutate select
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
        separate(value, into = c("resource", "file"), remove = FALSE) %>%
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
                                )
                        }),
                    result)) %>%
        deframe() %>%
        plyr::rename(., c("lrc2p" = "Default"), # change this to default
                     warn_missing = FALSE)
}
