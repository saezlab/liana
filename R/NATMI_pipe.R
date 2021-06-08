#' Call NATMI Pipeline from R with OmniPath
#' @param omni_resources List of OmniPath resources
#' @param omnidbs_path path of saved omnipath resources
#' @param natmi_path path of NATMI code and dbs
#' @param em_path expression matrix path
#' @param ann_path annotations (i.e. clusters) path
#' @param output_path NATMI output path
#' @param .format bool whether to format output
#' @param .write_data bool whether Extract data from Seurat Object
#' @param .default_run bool whether to run default DBs or not
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
    omni_resources,
    seurat_object = NULL,
    omnidbs_path = "input/omnipath_NATMI",
    natmi_path = "NATMI/",
    em_path = "input/test_em.csv",
    ann_path = "input/test_metadata.csv",
    output_path = "output/NATMI_test",
    .assay = "RNA",
    .format = TRUE,
    .write_data = TRUE,
    .seed = 1004,
    .num_cor = 4){

    py_set_seed(.seed)
    .wdir <- system.file(package = 'intercell')
    .natmi_dir <- system.file('NATMI/', package = 'intercell')

    # append default resources to OmniPath ones
    if("DEFAULT" %in% toupper(names(omni_resources))){
        omni_resources %<>% purrr::list_modify("Default" = NULL)
        resource_names <- append(as.list(names(omni_resources)),
                                "lrc2p"

        )
    }
    resource_names

    if(.write_data){
        log_info(str_glue("Writing EM to {em_path}"))
        write.csv(100 * (exp(as.matrix(GetAssayData(object = seurat_object,
                                                    assay = .assay,
                                                    slot = "data"))) - 1),
                  file = em_path,
                  row.names = TRUE)
        log_info(str_glue("Writing Annotations to {ann_path}"))
        write.csv(Idents(seurat_object)  %>%
                      enframe(name="barcode", value="annotation"),
                  file = ann_path,
                  row.names = FALSE)

        log_info(str_glue("Saving resources to {omnidbs_path}"))
        omni_to_NATMI(omni_resources, omnidbs_path)
    }

    log_success(str_glue("Output to be saved and read from {output_path}"))
    dir.create(file.path(output_path), recursive = TRUE)

    # copy OmniPath resources to NATMI dir
    file.copy(list.files(omnidbs_path, "*.csv$",
                         full.names = TRUE),
              to=str_glue("{natmi_path}/lrdbs/"),
              overwrite = TRUE)


    # submit native sys requests
    # Check issue with Default
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


    # load results
    natmi_results <- FormatNatmi(output_path, .format)

    return(natmi_results)
}



#' Reform OmniPath Resource to NATMI format and save to location
#' @param resource_names list of omnipath resources
#' @param omni_path directory in which to save OP resources
omni_to_NATMI <- function(omni_resources,
                          omni_path = "input/omnipath_NATMI"){

    omni_resources %<>% purrr::list_modify("Default" = NULL)

    names(omni_resources) %>%
        map(function(x){
            write.csv(omni_resources[[x]]  %>%
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
FormatNatmi <- function(output_path, .format = TRUE){

    list.files(output_path,
               all.files = TRUE,
               recursive = TRUE,
               pattern ="Edges_") %>%
        enframe() %>%
        separate(value, into = c("resource", "file"), remove = FALSE) %>%
        mutate(value =  value %>% map(function(csv)
            read.csv(str_glue("{output_path}/{csv}")))) %>%
        select(resource, "result" = value) %>%
        mutate(result = if_else(rep(.format, length(.data$result)), result %>% map(function(df){
            df %>% select(source = Sending.cluster,
                          target = Target.cluster,
                          ligand = Ligand.symbol,
                          receptor = Receptor.symbol,
                          edge_avg_expr = Edge.average.expression.weight,
                          edge_specificity = Edge.average.expression.derived.specificity)
        }), result)) %>%
        deframe() %>%
        plyr::rename(., c("lrc2p" = "Default")) # change this to default
}


#' Split and format strings for subsampling
#'
#' @param path path to CSV (em/annotations) to split and format to the
#'     current subsampling taken from a seurat object project name
#' @param project_name seurat object project name
#' @return Path to save subsampling EM and Annotations
str_split_helper <- function(path, project_name){
    split_path <- str_split(path, pattern = "\\.", n = 2)[[1]][1]
    str_glue("{split_path}_{project_name}.csv")
}
