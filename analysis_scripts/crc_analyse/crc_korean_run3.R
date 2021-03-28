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


# 5. iTALK
italk_results <- omni_resources %>%
    map(function(db)
        call_italk(op_resource = db,
                   crc_korean,
                   assay = 'RNA',
                   .format = TRUE,
                   .DE = TRUE
        ))
saveRDS(italk_results, "output/crc_res/italk_results.rds")


