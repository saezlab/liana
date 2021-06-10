seurat_object <- readRDS("input/visium_converted/cere_visium.rds")

# Get Full Omni Resources
omni_resources <- compile_ligrec()
# saveRDS(omni_resources, "input/omni_resources.rds")
omni_resources <- readRDS("input/omni_resources.rds")
omni_resources <- omni_resources[c("OmniPath", "Reshuffled")]

# 1. Squidpy -------------------------------------------------------------------
squidpy_results <- call_squidpyR(seurat_object = seurat_object,
                                 omni_resources = omni_resources,
                                 python_path = "/home/dbdimitrov/anaconda3/bin/python",
                                 n_perms=10000,
                                 threshold=0.1,
                                 seed=as.integer(1004),
                                 n_jobs=as.integer(10),
                                 cluster_key="seurat_clusters")
saveRDS(squidpy_results, "output/visium_runs/cere/squidpy_full.rds")


# 2. NATMI --------------------------------------------------------------------
# save OmniPath Resource to NATMI format
natmi_results <- call_natmi(omni_resources = omni_resources,
                            seurat_object = seurat_object,
                            omnidbs_path = "~/Repos/ligrec_decoupleR/input/omnipath_NATMI",
                            natmi_path = "~/Repos/NATMI",
                            em_path = "~/Repos/ligrec_decoupleR/input/cere_visium_em.csv",
                            ann_path = "~/Repos/ligrec_decoupleR/input/cere_visium_ann.csv",
                            output_path = "~/Repos/ligrec_decoupleR/output/visium_runs/cere/natmi_run",
                            .write_data = TRUE,
                            .subsampling_pipe = FALSE,
                            .assay = "SCT"
)
saveRDS(natmi_results, "output/visium_runs/cere/natmi_full.rds")

# 3. CellChat -----------------------------------------------------------------
cellchat_results <- omni_resources %>%
    map(function(db) call_cellchat(op_resource = db,
                                   seurat_object = seurat_object,
                                   nboot = 1000,
                                   exclude_anns = c(),
                                   thresh = 1,
                                   assay = "SCT",
                                   .normalize = FALSE,
                                   .do_parallel = FALSE,
                                   .raw_use = TRUE)) %>%
    setNames(names(omni_resources))
saveRDS(cellchat_results, "output/visium_runs/cere/cellchat_full.rds")


# 4. SCA ----------------------------------------------------------------------
sca_results <- omni_resources %>%
    map(function(db)
        call_sca(op_resource = db,
                 seurat_object,
                 assay = 'SCT',
                 .format = TRUE,
                 s.score = 0,
                 logFC = log2(1.5)
        ))
saveRDS(sca_results, "output/visium_runs/cere/sca_full.rds")


# 5. Connectome ----------------------------------------------------------------
conn_results <- omni_resources %>%
    map(function(db)
        call_connectome(seurat_object = seurat_object,
                        op_resource = db,
                        min.cells.per.ident = 10,
                        p.values = TRUE,
                        calculate.DOR = FALSE,
                        .format = TRUE,
                        assay = 'SCT'))
saveRDS(conn_results, "output/visium_runs/cere/conn_full.rds")


# 6. iTALK
italk_results <- omni_resources %>%
    map(function(db)
        call_italk(op_resource = db,
                   seurat_object,
                   assay = 'SCT',
                   .format = TRUE,
                   .DE = TRUE
        ))
saveRDS(italk_results, "output/visium_runs/cere/italk_full.rds")
