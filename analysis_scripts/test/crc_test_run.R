# Load Data
crc_korean <- readRDS("input/crc_data/crc_korean.rds") %>%
    format_crc_meta()
# saveRDS(crc_korean, "input/crc_data/crc_korean_mod.rds")

crc_korean <- subset(crc_korean, cells = rownames(crc_korean@meta.data)[5000:7500])
crc_korean <- subset(crc_korean, cells = rownames(crc_korean@meta.data))

# Get Full Omni Resources
# omni_resources <- compile_ligrec()
# saveRDS(omni_resources, "input/omni_resources.rds")
omni_resources <- readRDS("input/omni_resources.rds")
omni_resources <- list("Kirouac2010" = omni_resources$Kirouac2010,
                       "ICELLNET" = omni_resources$ICELLNET)


# 1. Squidpy -------------------------------------------------------------------
squidpy_results <- call_squidpyR(seurat_object = readRDS("input/crc_data/crc_korean_mod.rds"),
                                 omni_resources = omni_resources,
                                 python_path = "/home/dbdimitrov/anaconda3/bin/python",
                                 .ident = "Cell_subtype")
saveRDS(squidpy_results, "output/crc_res/squidpy_results.rds")

# 2. NATMI --------------------------------------------------------------------
# save OmniPath Resource to NATMI format
natmi_results <- call_natmi(omni_resources = omni_resources,
                            seurat_object = crc_korean,
                            wd_path = "/home/dbdimitrov/Repos/ligrec_decoupleR",
                            omnidbs_path = "~/Repos/ligrec_decoupleR/input/omnipath_NATMI",
                            natmi_path = "~/Repos/NATMI",
                            em_path = "~/Repos/ligrec_decoupleR/input/crc_korean_em.csv",
                            ann_path = "~/Repos/ligrec_decoupleR/input/crc_korean_ann.csv",
                            output_path = "~/Repos/ligrec_decoupleR/output/crc_res/natmi_results.rds",
                            .write_data = FALSE,
                            .subsampling_pipe = FALSE,
                            .assay = "RNA",
                            .num_cor = 11
                            )
saveRDS(natmi_results, "output/crc_res/natmi_results.rds")


# 3. CellChat -----------------------------------------------------------------
cellchat_results <- omni_resources %>%
    map(function(db) call_cellchat(op_resource = db,
                                   seurat_object = crc_korean,
                                   nboot = 10,
                                   exclude_anns = c(),
                                   thresh = 1,
                                   assay = "RNA",
                                   .normalize = TRUE,
                                   .do_parallel = TRUE)) %>%
    setNames(names(omni_resources))
saveRDS(cellchat_results, "output/crc_res/cellchat_results.rds")

# 4. SCA ----------------------------------------------------------------------
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

# 5. Connectome ----------------------------------------------------------------
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

# 6. iTALK
italk_results <- omni_resources %>%
    map(function(db)
        call_italk(op_resource = db,
                   crc_korean,
                   assay = 'RNA',
                   .format = TRUE,
                   .DE = TRUE
        ))
saveRDS(italk_results, "output/crc_res/italk_results.rds")
