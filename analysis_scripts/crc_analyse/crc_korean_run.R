library(intercell)
sapply(list.files("/net/data.isilon/ag-saez/bq_ddimitrov/Repos/Cell_Cell_Investigation/R/", pattern = ".R", full.names = TRUE), source)
setwd("/net/data.isilon/ag-saez/bq_ddimitrov/Repos/Cell_Cell_Investigation/")


# Load Data
crc_korean <- readRDS("input/crc_data/crc_korean.rds") %>%
    format_crc_meta()

# Get Full Omni Resources
# omni_resources <- compile_ligrec()
# saveRDS(omni_resources, "input/omni_resources.rds")
omni_resources <- readRDS("input/omni_resources.rds")


# 1. Squidpy -------------------------------------------------------------------
squidpy_results <- call_squidpyR(seurat_object = crc_korean,
                                 omni_resources = omni_resources,
                                 python_path = "/net/data.isilon/ag-saez/bq_ddimitrov/SOFTWARE/miniconda3/envs/ligrec/bin/python3.8",
                                 .ident = "Cell_subtype")
saveRDS(squidpy_results, "output/crc_res/squidpy_results.rds")

# 2. NATMI --------------------------------------------------------------------
# save OmniPath Resource to NATMI format
natmi_results <- call_natmi(omni_resources = omni_resources,
                            seurat_object = crc_korean,
                            wd_path = "/net/data.isilon/ag-saez/bq_ddimitrov/Repos/Cell_Cell_Investigation/",
                            omnidbs_path = "/net/data.isilon/ag-saez/bq_ddimitrov/Repos/Cell_Cell_Investigation/input/omni_natmi",
                            natmi_path = "/net/data.isilon/ag-saez/bq_ddimitrov/Repos/NATMI/",
                            em_path = "/net/data.isilon/ag-saez/bq_ddimitrov/Repos/Cell_Cell_Investigation/input/crc_korean_counts.csv",
                            ann_path = "/net/data.isilon/ag-saez/bq_ddimitrov/Repos/Cell_Cell_Investigation/input/crc_korean_ann.csv",
                            output_path = "/net/data.isilon/ag-saez/bq_ddimitrov/Repos/Cell_Cell_Investigation/output/crc_natmi_test/",
                            .write_data = TRUE,
                            .subsampling_pipe = FALSE,
                            .assay = "RNA"
)
saveRDS(natmi_results, "output/crc_res/natmi_results.rds")

