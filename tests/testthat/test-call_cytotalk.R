liana_path <- system.file(package = "liana")
testdata <- Seurat::NormalizeData(readRDS(file.path(liana_path , "testdata", "input", "testdata.rds")))
op_resource <- decomplexify(select_resource("OmniPath")[[1]])

call_cytotalk(SeuratObject = testdata, op_resouce = op_resource)

