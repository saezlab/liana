#' Call NATMI Pipeline from R with OmniPath
#' @param omni_resources List of OmniPath resources
#' @param omnidbs_path path of saved omnipath resources
#' @param natmi_path path of NATMI code and dbs
#' @param em_path expression matrix path
#' @param ann_path annotations (i.e. clusters) path
#' @param output_path NATMI output path
#' @param .format bool whether to format output
#' @return DF with NATMI results
#' @details This function will take omnipath resources saved to csvs and copy them to the
#' NATMI dbs folder, then it will natively call the python modules of NATMI
#' in the NATMI dir and save the output into a specified directory.
#' It will then load and format the output to a DF.
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
#'
#' Stat details:
#' 1) The mean-expression edge weights are calculated by multiplying the mean-expression level of the ligand in the
#'   sending cell type by the mean expression of the receptor in the target cell type (no discriminatory information)
#' 2) The specificity-based edge weights, help identify the most specific edges in the network
#'  where each specificity is defined as the mean expression of the ligand/receptor in a given cell type
#'  divided by the sum of the mean expression of that ligand/receptor across all cell types

call_natmi <- function(omni_resources,
               omnidbs_path = "~/Repos/ligrec_decoupleR/input/omnipath_NATMI",
               natmi_path = "~/Repos/NATMI",
               em_path = "~/Repos/ligrec_decoupleR/input/test_em.csv",
               ann_path = "~/Repos/ligrec_decoupleR/input/test_metadata.csv",
               output_path = "~/Repos/ligrec_decoupleR/output/NATMI_test",
               .format = TRUE){


    library(rprojroot)
    project_rootdir <- find_rstudio_root_file()

    # copy OmniPath resources to NATMI dir
    file.copy(list.files(omnidbs_path, "*.csv$",
                         full.names = TRUE),
              to=str_glue("{natmi_path}/lrdbs/"),
              overwrite = TRUE)

    # set current dir to NATMI
    setwd(natmi_path)

    # append default resources to OmniPath ones
    omni_list <- append(as.list(names(omni_resources)),
                        list("lrc2p", "lrc2a"))


    # submit native sys requests
    omni_list %>% map(function(resource){
        system(str_glue("python ExtractEdges.py ",
                        "--species human ",
                        "--emFile {em_path} ",
                        "--annFile {ann_path} ",
                        "--interDB {resource} ",
                        "--coreNum 8 ",
                        "--out {output_path}/{resource}",
                        sep = " "))
    })

    # set dir back to project
    setwd(project_rootdir)

    # load results
    natmi_results <- list.files(output_path,
                                all.files = TRUE,
                                recursive = TRUE,
                                pattern ="Edges_") %>%
        enframe() %>%
        separate(value, into = c("resource", "file"), remove = FALSE) %>%
        mutate(value =  value %>% map(function(csv)
            read_csv(str_glue("{output_path}/{csv}")))) %>%
        select(resource, "result" = value) %>%
        mutate(result = ifelse(.format, result %>% map(function(df){
            df %>% select(source = `Sending cluster`,
                          target = `Target cluster`,
                          ligand = `Ligand symbol`,
                          receptor = `Receptor symbol`,
                          edge_avg_expr = `Edge average expression weight`,
                          edge_specificity = `Edge average expression derived specificity`)
        }), result)) %>%
        deframe()

    return(natmi_results)
}



#' Reform OmniPath Resource to NATMI format and save to location
#' @param omni_list list of omnipath resources
#' @param omni_path directory in which to save OP resources
omni_to_NATMI <- function(omni_resources,
                          omni_path = "input/omnipath_NATMI"){

    names(omni_resources) %>%
        map(function(x){
            write.csv(omni_resources[[x]]  %>%
                          # filter complexes
                          filter(entity_type_intercell_source != "complex",
                                 entity_type_intercell_target != "complex") %>%
                          select("Ligand gene symbol" = source_genesymbol,
                                 "Receptor gene symbol" = target_genesymbol) %>%
                          as.data.frame(),
                      file = str_glue("{omni_path}/{x}.csv"),
                      row.names = FALSE)
        })
}
