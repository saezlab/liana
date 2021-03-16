# Load Data
breast_cancer <- readRDS("input/sc_bc/breast_cancer_seurat323.rds")

# Get Full Omni Resources
# omni_resources <- compile_ligrec()
# saveRDS(omni_resources, "input/omni_resources.rds")
# omni_resources <- readRDS("input/omni_resources.rds")



# 1. Squidpy -------------------------------------------------------------------
squidpy_results <- call_squidpyR(seurat_object = breast_cancer,
                                 omni_resources = omni_resources,
                                 python_path = "/home/dbdimitrov/anaconda3/bin/python",
                                 .ident = "seurat_clusters")
saveRDS(squidpy_results, "output/benchmark/main_run/squidpy_full.rds")


# 2. NATMI --------------------------------------------------------------------
# save OmniPath Resource to NATMI format
natmi_results <- call_natmi(omni_resources = omni_resources,
                            seurat_object = breast_cancer,
                            omnidbs_path = "~/Repos/ligrec_decoupleR/input/omnipath_NATMI",
                            natmi_path = "~/Repos/NATMI",
                            em_path = "~/Repos/ligrec_decoupleR/input/bc_em.csv",
                            ann_path = "~/Repos/ligrec_decoupleR/input/bc_ann.csv",
                            output_path = "~/Repos/ligrec_decoupleR/output/benchmark/natmi_full",
                            .write_data = TRUE,
                            .subsampling_pipe = FALSE
                            )
saveRDS(natmi_results, "output/benchmark/main_run/natmi_full.rds")

# 3. CellChat -----------------------------------------------------------------
cellchat_results <- omni_resources %>%
    map(function(db) call_cellchat(op_resource = db,
                                  seurat_object = breast_cancer,
                                  nboot = 100,
                                  exclude_anns = c(),
                                  thresh = 1,
                                  assay = "SCT",
                                  .normalize = FALSE)) %>%
    setNames(names(omni_resources))
saveRDS(cellchat_results, "output/benchmark/main_run/cellchat_full.rds")


call_cellchat(op_resource = omni_resources$Random,
              seurat_object = breast_cancer,
              nboot = 100,
              exclude_anns = c(),
              thresh = 1,
              assay = "SCT",
              .normalize = FALSE)



# 4. SCA ----------------------------------------------------------------------
sca_results <- omni_resources_plus %>%
    map(function(db)
        call_sca(op_resource = db,
                 breast_cancer,
                 assay = 'SCT',
                 .format = TRUE,
                 s.score = 0,
                 logFC = 0.25
        ))
# saveRDS(sca_results, "output/benchmark/main_run/sca_full.rds")


# 5. Connectome ----------------------------------------------------------------
conn_results <- omni_resources_plus %>%
    map(function(db)
    call_connectome(seurat_object = breast_cancer,
                    op_resource = db,
                    min.cells.per.ident = 10,
                    p.values = TRUE,
                    calculate.DOR = FALSE,
                    .format = FALSE,
                    assay = 'SCT'))
# saveRDS(conn_results, "output/benchmark/main_run/conn_full.rds")


# 6. iTALK
italk_results <- omni_resources_plus %>%
    map(function(db)
        call_italk(op_resource = db,
                   breast_cancer,
                   assay = 'SCT',
                   .format = TRUE
        ))
# saveRDS(italk_results, "output/benchmark/main_run/italk_full.rds")
