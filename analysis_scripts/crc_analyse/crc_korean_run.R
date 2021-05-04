require(intercell)
require(tibble)
require(magrittr)

# Load Data and Format data
# crc_korean <- readRDS("input/crc_data/crc_korean.rds") %>%
#     format_crc_meta()
# saveRDS(crc_korean, "input/crc_data/crc_korean_form.rds")
crc_korean <- readRDS("input/crc_data/crc_korean_form.rds")

# Get Full Omni Resources
# omni_resources <- compile_ligrec(lr_pipeline = TRUE)
# saveRDS(omni_resources, "input/omni_resources.rds")
omni_resources <- readRDS("input/omni_resources.rds")

# 1. Squidpy -------------------------------------------------------------------
# squidpy_results <- call_squidpyR(seurat_object = crc_korean,
#                                  omni_resources = omni_resources,
#                                  python_path = "/net/data.isilon/ag-saez/bq_ddimitrov/SOFTWARE/miniconda3/envs/ligrec/bin/python3.8",
#                                  .ident = "Cell_subtype")
# saveRDS(squidpy_results, "output/crc_res/squidpy_results_local.rds")


# 2. NATMI --------------------------------------------------------------------
# save OmniPath Resource to NATMI format
natmi_results <- call_natmi(omni_resources = omni_resources,
                            seurat_object = crc_korean,
                            wd_path = "/net/data.isilon/ag-saez/bq_ddimitrov/Repos/Cell_Cell_Investigation/",
                            omnidbs_path = "/net/data.isilon/ag-saez/bq_ddimitrov/Repos/Cell_Cell_Investigation/input/omni_natmi",
                            natmi_path = "/net/data.isilon/ag-saez/bq_ddimitrov/Repos/NATMI/",
                            em_path = "/net/data.isilon/ag-saez/bq_ddimitrov/Repos/Cell_Cell_Investigation/input/crc_korean_counts.csv",
                            ann_path = "/net/data.isilon/ag-saez/bq_ddimitrov/Repos/Cell_Cell_Investigation/input/crc_korean_ann.csv",
                            output_path = "/net/data.isilon/ag-saez/bq_ddimitrov/Repos/Cell_Cell_Investigation/output/crc_natmi/",
                            .write_data = TRUE,
                            .subsampling_pipe = FALSE,
                            .assay = "RNA",
                            .num_cor = 64
                            )
saveRDS(natmi_results, "output/crc_res/natmi_results.rds")


# 3. SCA ----------------------------------------------------------------------
sca_results <- omni_resources %>%
    map(function(db)
        call_sca(op_resource = db,
                 crc_korean,
                 assay = 'RNA',
                 .format = TRUE,
                 s.score = 0,
                 logFC = log2(1.5)
        ))
saveRDS(sca_results, "output/crc_res/sca_results.rds")


# 4. Connectome ----------------------------------------------------------------
conn_results <- omni_resources %>%
    map(function(db)
        call_connectome(seurat_object = crc_korean,
                        .spatial = FALSE,
                        op_resource = db,
                        min.cells.per.ident = 1,
                        p.values = TRUE,
                        calculate.DOR = FALSE,
                        .format = TRUE,
                        assay = 'RNA'))
saveRDS(conn_results, "output/crc_res/conn_results.rds")


# 5. iTALK ---------------------------------------------------------------------
italk_results <- omni_resources %>%
    map(function(db)
        call_italk(op_resource = db,
                   seurat_object = crc_korean,
                   assay = 'RNA',
                   .format = TRUE,
                   .DE = TRUE
        ))
saveRDS(italk_results, "output/crc_res/italk_results.rds")



# 6. CellChat -----------------------------------------------------------------
# cellchat_results <- omni_resources %>%
#     map(function(db) call_cellchat(op_resource = db,
#                                    seurat_object = crc_korean,
#                                    nboot = 1000,
#                                    exclude_anns = c(),
#                                    thresh = 1,
#                                    assay = "RNA",
#                                    .normalize = FALSE,
#                                    .do_parallel = FALSE,
#                                    .raw_use = TRUE
#                                    )) %>%
#     setNames(names(omni_resources))
# saveRDS(cellchat_results, "output/crc_res/cellchat_results.rds")
