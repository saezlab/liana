#' Call NATMI Pipeline from R with Resources Querried from OmniPath [[DEPRECATED]]
#'
#' @param op_resource List of OmniPath resources
#' @param sce Seurat or SingleCellExperiment object
#' @param expr_file expression matrix file name
#' @param meta_file annotations (i.e. clusters) file name
#' @param output_dir NATMI output directory
#' @param .format bool whether to format output
#' @param .overwrite_data bool whether Extract and overwrite csv with
#'    data from Seurat Object
#' @param .natmi_path path of NATMI code and dbs (by default set to liana path)
#' @param assay Seurat assay to be used
#' @param reso_name name of the resource usually in the format list(name = op_resource)
#' @param num_cor number of cores to be used
#' @param conda_env name of the conda environment via which NATMI is called
#' @param assay.type logcounts by default, but it's converted back into counts
#' as suggested by the authors
#' @param .seed random seed
#' @param .delete_input_output logical whether to delete input and
#'  output after run.
#'
#'
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
#' Note that `call_natmi` will write the expression matrix to CSV each time
#' its called, unless .overwrite_data is set to FALSE!
#' This can be an extremely time consuming step when working with large datasets
#'
#' Also, NATMI will sometimes create duplicate files, so please consider
#'  saving each run in a new folder. An easy fix would be to simply delete the
#'  output, but I am reluctant to automatically delete files via an R script.
#'
#'
#' @importFrom reticulate py_set_seed
#' @importFrom stringr str_glue
#' @importFrom SeuratObject GetAssayData Idents
#' @import dplyr reticulate
#'
#' @export
call_natmi <- function(
    sce,
    op_resource,
    expr_file = "em.csv",
    meta_file = "metadata.csv",
    output_dir = "NATMI_test",
    reso_name = "placeholder",
    assay = "RNA",
    num_cor = 4,
    conda_env = NULL,
    assay.type = "logcounts",
    .format = TRUE,
    .overwrite_data = TRUE,
    .seed = 1004,
    .natmi_path = NULL,
    .delete_input_output = FALSE){

    # Convert sce to seurat
    if(class(sce) == "SingleCellExperiment"){
        sce %<>% .liana_convert(., assay=assay)
    }

    # Get Reticulate path
    reticulate::use_condaenv(condaenv = conda_env %>% `%||%`("liana_env"),
                             conda = "auto",
                             required = TRUE)
    py$pd <- reticulate::import("pandas")
    python_path <- reticulate::py_discover_config()[[1]]
    reticulate::py_set_seed(.seed)

    .natmi_path %<>% `%||%`(system.file('NATMI/', package = 'liana'))
    .input_path = file.path(.natmi_path, 'data', 'input')
    .output_path = file.path(.natmi_path, 'data', 'output', output_dir)
    .csv_path = file.path(.input_path, expr_file)

    if(!dir.exists(file.path(.input_path))){
        print(str_glue("Input path created: {.input_path}"))
        dir.create(file.path(.input_path), recursive = TRUE)
    }

        print(str_glue("Output path created: {.output_path}"))
        dir.create(file.path(.output_path), recursive = TRUE)

    if(.overwrite_data || !file.exists(.csv_path)){
        print(str_glue("Writing EM to {.csv_path}"))
        if(assay.type=="counts"){
            write.csv(GetAssayData(object = sce,
                                   assay = "RNA",
                                   slot = "counts"),
                      file = .csv_path,
                      row.names = TRUE)
        } else{
            write.csv(100 * (exp(as.matrix(
                GetAssayData(object = sce,
                             assay = assay,
                             slot = "data"))) - 1),
                file = .csv_path,
                row.names = TRUE)
        }
    }

    print(str_glue("Writing Annotations to {.input_path}/{meta_file}"))
    write.csv(Idents(sce) %>%
                  enframe(name="barcode", value="annotation"),
              file = file.path(.input_path, meta_file),
              row.names = FALSE)

    print(str_glue("Saving resource to {.natmi_path}/lrdbs/{reso_name}"))
    # Deal with Default (i.e. NULL)
    if(is.null(op_resource)){
        reso_name <- "lrc2p"

    } else{
        # save resource to NATMI dir
        omni_to_NATMI(op_resource,
                      reso_name,
                      file.path(.natmi_path, "lrdbs"))
    }

    # submit native sys request
    system(str_glue("{python_path} {.natmi_path}/ExtractEdges.py ",
                    "--species human ",
                    "--emFile {.csv_path} ",
                    "--annFile {.input_path}/{meta_file} ",
                    "--interDB {.natmi_path}/lrdbs/{reso_name}.csv ",
                    "--coreNum {num_cor} ",
                    "--out {.output_path}/{reso_name}",
                    sep = " "))

    # load and format results
    natmi_results <- FormatNatmi(.output_path, reso_name, .format)

    if(.delete_input_output){
        system(str_glue("rm -r {.output_path}/{reso_name}"))
        system(str_glue("rm -r {.input_path}/"))
    }


    return(natmi_results)
}



#' Reform OmniPath Resource to NATMI format and save to location
#'
#' @param op_resource Resource formatted as OmniPath
#' @param reso_name name of the resource
#' @param natmi_db_path directory in which to save the resource
#'
#' @noRd
omni_to_NATMI <- function(op_resource,
                          reso_name = "placeholder",
                          natmi_db_path = "input/omnipath_NATMI"){

    write.csv(op_resource %>%
                  select("Ligand gene symbol" = source_genesymbol,
                         "Receptor gene symbol" = target_genesymbol) %>%
                  distinct() %>%
                  as.data.frame(),
              file = str_glue("{natmi_db_path}/{reso_name}.csv"),
              row.names = FALSE)
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
                 extra = "drop",
                 sep="\\/") %>%
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
