
library(RSoptSC)

## Load and Preprocess Data
# path to the file containing counts
df <- system.file("extdata", "GSE67602_JoostData.csv.bz2", package = "RSoptSC")
# path the file containing unique gene ids
gf <- system.file("extdata", "GSE67602_JoostGenes.csv.bz2", package = "RSoptSC")

# path the file containing unique cell ids
cf <- system.file("extdata", "GSE67602_JoostCells.csv.bz2", package = "RSoptSC")

# path to the file containing labels, e.g. tissue type of the cells
af <- system.file("extdata", "GSE67602_JoostAnnotation.csv.bz2", package = "RSoptSC")

GSE67602_Joost <- LoadData(df, gf, cf, af)


# Remove Spike-in RNA
# Apply number of features and exclusion threshold
logdata <- log10(data + 1)
gene_expression_threshold <- 0.03
n_features <- 3000
filtered_data <- SelectData(logdata, gene_expression_threshold, n_features)

?SelectData


#
library(knitr)
library(kableExtra)

assayData(seurat_object)

xd <- as.matrix(Seurat::GetAssayData(seurat_object))

lig_rec_path <- system.file("extdata", "tgfb_lig_rec.tsv", package = "RSoptSC")
rec_target_path <- system.file("extdata", "tgfb_rec_target_both.tsv", package = "RSoptSC")

xd <- read.table(lig_rec_path, header = T)
xd2 <- read.table(rec_target_path, header = T)

xd2


?ImportPathway
pathway <- ImportPathway(lig_table_path = lig_rec_path,
                         rec_table_path = rec_target_path,
                         data = xd,
                         gene_names = gene_names)



pathway$pathway %>% kable() %>% kable_styling()
