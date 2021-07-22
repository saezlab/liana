# input
liana_path <- system.file(package = "liana")
seurat_object <-
    readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))

require(SingleCellExperiment)

# Run /w OmniPath
lr_res <- liana_pipe(seurat_object,
                     select_resource("OmniPath")[[1]],
                     decomplexify = TRUE)
lr_res

# Run /w CellPhoneDB
lr_cdbd <- liana_pipe(seurat_object,
                      select_resource("CellPhoneDB")[[1]])

lr_cdbd

# call_connectome
liana_scores(.score_specs()[["connectome"]],
             lr_res = lr_res)

get_connectome(seurat_object,
               select_resource("OmniPath")[[1]],
               lr_res = lr_res,
               protein = 'subunit')


get_natmi(seurat_object,
          select_resource("OmniPath")[[1]],
          lr_res = lr_res)

get_logfc(seurat_object,
          select_resource("OmniPath")[[1]],
          lr_res = lr_res)
