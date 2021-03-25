# Load Data
crc_belgian <- readRDS("input/crc_data/crc_belgian.rds")
crc_belgian@meta.data <-  crc_belgian@meta.data %>%
    unite(Cell_type, Cell_subtype, col = "Cell_clusters", remove = FALSE, sep="_") %>%
    # Filter Unknown and Unspecified
    filter(!(str_detect("Unknown", Cell_subtype))) %>%
    filter(!(str_starts("Unspecified Plasma", Cell_subtype))) %>%
    # Myoloids_SPP1+A/B into Myoloids_SPP1+ (+ fix typo)
    mutate(Cell_clusters = if_else(str_detect(Cell_clusters, "SPP1"),
                        "Myoloids_SPP1+",
                        Cell_clusters)
    ) %>%
    # Group Healthy Epithelial Cells
    mutate(Cell_clusters = if_else((str_detect(Cell_clusters, "Epithelial") & !str_detect(Cell_clusters, "CMS")),
                        "Epithelial_Healthy",
                        Cell_clusters)
           ) %>%
    mutate(Cell_clusters = factor(Cell_clusters)) %>%
    mutate(Cell_type = factor(Cell_type))
crc_belgian <- subset(crc_belgian, cells = rownames(crc_belgian@meta.data))
crc_belgian <- SetIdent(crc_belgian, value = crc_belgian@meta.data$Cell_type)
crc_belgian <- subset(crc_belgian, cells = rownames(crc_belgian@meta.data)[5000:10000])

# Get Full Omni Resources
# omni_resources <- compile_ligrec()
# saveRDS(omni_resources, "input/omni_resources.rds")
omni_resources <- readRDS("input/omni_resources.rds")
omni_resources <- list("Kirouac2010" = omni_resources$Kirouac2010,
                       "ICELLNET" = omni_resources$ICELLNET)


# 1. Squidpy -------------------------------------------------------------------
squidpy_results <- call_squidpyR(seurat_object = crc_belgian,
                                 omni_resources = omni_resources,
                                 python_path = "/home/dbdimitrov/anaconda3/bin/python",
                                 .ident = "Cell_type")

# 2. NATMI --------------------------------------------------------------------
# save OmniPath Resource to NATMI format
natmi_results <- call_natmi(omni_resources = omni_resources,
                            seurat_object = crc_belgian,
                            wd_path = "/home/dbdimitrov/Repos/ligrec_decoupleR",
                            omnidbs_path = "~/Repos/ligrec_decoupleR/input/omnipath_NATMI",
                            natmi_path = "~/Repos/NATMI",
                            em_path = "~/Repos/ligrec_decoupleR/input/crc_em.csv",
                            ann_path = "~/Repos/ligrec_decoupleR/input/crc_ann.csv",
                            output_path = "~/Repos/ligrec_decoupleR/output/benchmark/natmi_crc",
                            .write_data = TRUE,
                            .subsampling_pipe = FALSE,
                            .assay = "RNA"
                            )

# 3. CellChat -----------------------------------------------------------------
cellchat_results <- omni_resources %>%
    map(function(db) call_cellchat(op_resource = db,
                                   seurat_object = crc_belgian,
                                   nboot = 10,
                                   exclude_anns = c(),
                                   thresh = 1,
                                   assay = "RNA",
                                   .normalize = TRUE,
                                   .do_parallel = TRUE)) %>%
    setNames(names(omni_resources))


# 4. SCA ----------------------------------------------------------------------
sca_results <- omni_resources %>%
    map(function(db)
        call_sca(op_resource = db,
                 crc_belgian,
                 assay = 'RNA',
                 .format = TRUE,
                 s.score = 0,
                 logFC = log2(1.5)
        ))


# 5. Connectome ----------------------------------------------------------------
conn_results <- omni_resources %>%
    map(function(db)
        call_connectome(seurat_object = crc_belgian,
                        .spatial = FALSE,
                        op_resource = db,
                        min.cells.per.ident = 1,
                        p.values = TRUE,
                        calculate.DOR = FALSE,
                        .format = TRUE,
                        assay = 'RNA'))


# 6. iTALK
italk_results <- omni_resources %>%
    map(function(db)
        call_italk(op_resource = db,
                   crc_belgian,
                   assay = 'RNA',
                   .format = TRUE,
                   .DE = TRUE
        ))
