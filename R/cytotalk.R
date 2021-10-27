#### PEM ####

#' Compute Preferential Expression Measure scores from SingleCellExperiment object
#'
#' Uses the information contained in a SingleCellExperiment object to stratify the
#' expression matrix by cell type, and then computes the PEM scores from the exponential
#' of the normalized expression values. IMPORTANT: This function expects an object
#' containing log-transformed values (not raw counts).
#'
#' @param sce A SingleCellExperiment object containing log-normalized expression
#' values and cell type annotations as colLabels
#' @param assay.type The name of the data slot containing the log-normalized expression
#' values in the SingleCellExperiment object
#'
#' @return A matrix containing the computed PEM values
#'
#' @import SingleCellExperiment SeuratObject
#'
compute_pem_scores <- function(sce, assay.type = "logcounts") {

  # extract normalized data
  norm_data <- t(as.matrix(sce@assays@data[[assay.type]]))

  # expm1(x) computes exp(x) - 1 accurately also for |x| << 1
  exp_data <- expm1(norm_data)

  # calculate mean per cell using expm1 data
  means_per_cell <- exp_data %>%
    aggregate(., list(colLabels(sce)), FUN = mean) %>%
    tibble::as_tibble() %>%
    dplyr::rename(celltype = Group.1) %>%
    pivot_longer(-celltype, names_to = "gene") %>%
    tidyr::pivot_wider(names_from = celltype,
                       id_cols = gene,
                       values_from = value) %>%
    column_to_rownames("gene")

  # calculate PEMs from the mean matrix
  pem_out <- pems_from_means(means_per_cell)

  return(pem_out)

}

#' PEM from mean matrix
#'
#' Computes the Preferential Expression Measure (PEM) scores from a matrix containing
#' genes as rows and cell types as samples. The values in the matrix are supposed
#' to contain the average expression values of the count matrix per cell, obtained
#' from the exponential of the log-transformed count matrix.
#'
#' @param means_per_cell A matrix containing the mean expression per cell type
#'
#' @return A matrix with the computed PEM scores
#'
pems_from_means <- function(means_per_cell) {

  # calculate sum of means for rows and columns
  sums_per_cell <- colSums(means_per_cell)
  sums_per_gene <- rowSums(means_per_cell)
  total_sums <- sum(sums_per_cell)

  # will lead to bugs later, could check downstream
  if (total_sums == Inf) {
    message("Inf introduced into PEM score, is the data log1p tranformed?")
    return(1)
  }

  # what proportion of this cell type's rowmean sum accounts for the whole?
  cell_type_prop <- sums_per_cell / total_sums

  # gene_proportions
  gene_prop <- sapply(cell_type_prop, function(x) x * sums_per_gene)

  # calculate PEM scores
  pem_out <- log10(means_per_cell / gene_prop)

  # set negative and NAs PEMs to 0
  pem_out[pem_out < 0 | is.na(pem_out)] <- 0

  return(pem_out)

}

#### NST ####

#' Compute non-self talk scores from SingleCellExperiment object
#'
#' This function computes the non-self talk scores for all the cell types contained
#' in a SingleCellExperiment object. It iterates through the ligand-receptor pairs
#' that are provided as input and calculates these entropy-based measurements for each
#' cell type using its normalized expression matrix, which should contain log-transformed
#' values.
#'
#' @param sce A SingleCellExperiment object containing log-normalized expression
#' values and cell type annotations as colLabels
#' @param ligand_receptor_df A data frame with, at least, two columns named 'ligand'
#' and 'receptor' containing the ligand-receptor pairs to evaluate
#' @param assay.type The name of the data slot containing the log-normalized expression
#' values in the SingleCellExperiment object
#'
#' @return A data frame with the computed non-self-talk score for each ligand-receptor
#' pair on each cell-type.
#'
compute_nst_scores <- function(sce, ligand_receptor_df, assay.type = "logcounts") {

  # extract normalized data
  norm_data <- as.matrix(sce@assays@data[[assay.type]])

  # iterate through cell types, calculating NST for each ligand-receptor
  cell_types <- levels(colLabels(sce))
  nst_scores <- lapply(cell_types, function(cell) {

    mat <- norm_data[, colLabels(sce) == cell]
    nst_tibble <- compute_nst_from_matrix(mat = mat, ligand_receptor_df = ligand_receptor_df) %>%
      mutate(cell_type = cell)
    return(nst_tibble)

  }) %>%
    bind_rows()

  return(nst_scores)

}

#' Compute non-self talk scores from matrix
#'
#' Calculates the non-self talk scores (also refered to as mutual information distances)
#' for each ligand-receptor pair using the normalized expression matrix for a given
#' cell type. It expects log-transformed expression values.
#'
#' @param mat A matrix containing the log-transformed normalized expression values.
#' @param ligand_receptor_df A data frame with, at least, two columns named 'ligand'
#' and 'receptor' containing the ligand-receptor pairs to evaluate
#'
#' @return A data frame with the non-self-talk score for each pair of genes (ligand-receptors).
#'
compute_nst_from_matrix <- function(mat, ligand_receptor_df) {

  # extract number of columns and bins from the gene expression matrix
  n_cols <- ncol(mat)
  n_bins <- sqrt(n_cols)

  # add noise
  is_zero <- rowSums(mat) == 0
  mat[is_zero,] <- do.call(rbind, lapply(rep(n_cols, sum(is_zero)), noise))

  # filter ligand receptor data frame
  keep <- pmap(ligand_receptor_df,
               function(ligand, receptor)
                 ligand %in% rownames(mat) &
                 receptor %in% rownames(mat)) %>%
    as.logical()

  # calculate mi distances
  nst_score <- pmap(ligand_receptor_df[keep,],
                    function(ligand, receptor) {
                      mi_dist <- compute_mi_dist(exp1 = mat[ligand, ],
                                                 exp2 = mat[receptor, ],
                                                 nbins = n_bins)
    return(mi_dist)

  }) %>%
    as.numeric()

  # return df
  out_df <- mutate(ligand_receptor_df[keep,], nst_score = nst_score)

  return(out_df)

}

#' Compute mutual information distance from expression vectors
#'
#' Given two normalized gene expression vectors, and a given number of bins, this
#' function uses the entropy package to compute the mutual information between the
#' two vectors. Values passed to this function should be log-transformed.
#'
#' @param exp1 A vector containing the normalized gene expression for the first
#' gene (ligand)
#' @param exp2 A vector containing the normalized gene expression for the first
#' gene (receptor)
#' @param nbins Number of bins to discretize the expression values.
#'
#' @return The mututal information distance, also refered to as non-self atalk score
#'
compute_mi_dist <- function(exp1, exp2, nbins) {

  y2d <- entropy::discretize2d(exp1, exp2, nbins, nbins)
  H1 <- entropy::entropy(rowSums(y2d), method = "MM")
  H2 <- entropy::entropy(colSums(y2d), method = "MM")
  H12 <- entropy::entropy(y2d, method = "MM")
  mi <- H1 + H2 - H12

  norm <- min(c(H1, H2))
  mi_dist <- -log10(mi / norm)

  return(mi_dist)

}

#' Generate random noise
#'
#' This function generates a vector of length n with a noise value inserted in a
#' random position of the vector. This is needed to compute the entropy-related
#' values.
#'
#' @param n Length of the vector to generate
#'
#' @return A vector of length n with noise in a random position
#'
noise <- function(n, noise_value = 1e-20) {
  rng <- seq_len(n)
  (sample(rng, 1) == rng) * noise_value
}

#### CROSS-TALK SCORE ####

#' Compute cross-talk score from a Seurat Object
#'
#' This function executes all the required functions to extract and compute the
#' cross-talk scores as defined by the Cytotalk authors. It firsts transforms the SeuratObject
#' into a SingleCellExperiment object and then joins the results to the output of
#' the liana_pipe function.
#'
#' @param SeuratObject Seurat object as input
#' @param op_resouce The resource of ligand-receptor interactions
#'
#' @return A data frame containing the resulting scores
#'
call_cytotalk <- function(SeuratObject, op_resouce) {

  # Seurat to SCE
  ligand_receptor_df <- op_resource[, c("source_genesymbol", "target_genesymbol")]
  colnames(ligand_receptor_df) <- c("ligand", "receptor")
  entity_genes <- union(op_resource$source_genesymbol, op_resource$target_genesymbol)
  sce <- seurat_to_sce(testdata, entity_genes = entity_genes, assay = "RNA")
  sce <- sce[, colSums(SingleCellExperiment::counts(sce)) > 0]

  # apply liana_pipe to obtain intermediate results
  liana_pipe_results <- liana_pipe(seurat_object = testdata, op_resource = op_resource)

  # apply cytotalk functions
  pem_scores <- compute_pem_scores(sce = sce, assay.type = "logcounts")
  nst_scores <- compute_nst_scores(sce = sce, ligand_receptor_df = ligand_receptor_df, assay.type = "logcounts")

  # filter nst scores to remove NAs
  nst_scores <- tidyr::drop_na(nst_scores)

  # join with liana pipe results results
  liana_pipe_cytotalk <- liana_pipe_results %>%
    join_means(means = pem_scores,
               source_target = "source",
               entity = "ligand",
               type = "pem") %>%
    join_means(means = pem_scores,
               source_target = "target",
               entity = "receptor",
               type = "pem") %>%
    left_join(nst_scores, by = c("source" = "cell_type", "ligand", "receptor")) %>%
    dplyr::rename(nst_score_source = nst_score) %>%
    left_join(nst_scores, by = c("target" = "cell_type", "ligand", "receptor")) %>%
    dplyr::rename(nst_score_target = nst_score)

  # calculate final scores and retrieve output
  liana_pipe_cytotalk %<>%
    mutate(expression_score = (ligand.pem + receptor.pem) / 2,
           nst = (nst_score_source + nst_score_target) / 2)

  # normalize using minmax
  es <- liana_pipe_cytotalk$expression_score
  nst <- liana_pipe_cytotalk$nst
  Nes <- (es - min(es, na.rm = TRUE)) / (max(es, na.rm = TRUE) - min(es, na.rm = TRUE))
  Nnst <- (nst - min(nst, na.rm = TRUE)) / (max(nst, na.rm = TRUE) - min(nst, na.rm = TRUE))
  cross_talk_score <- Nes * Nnst

  # prepare out df
  out_df <-  liana_pipe_cytotalk %>%
    dplyr::select(ligand, receptor, source, target, ends_with("pem"), starts_with("nst"),
                  expression_score, nst) %>%
    mutate(normalized_expression_score = Nes, normalized_nst = Nnst, cross_talk_score = cross_talk_score)

  return(out_df)

}

