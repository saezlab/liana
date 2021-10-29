# CytoTalk Test -----
### LIANA
liana_path <- system.file(package = "liana")
seurat_object <- readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))
assay.type = "logcounts"
op_resource <- select_resource("OmniPath")[[1]] %>% decomplexify()

# Get cytotalk scores
martin_crosstalk <- call_cytotalk(seurat_object = seurat_object,
                                 op_resouce = op_resource)
martin_crosstalk

# LIANA pipe out
liana_out <- liana_pipe(seurat_object = seurat_object,
                        op_resource = op_resource,
                        expr_prop = 0)

# Basic CCLR scaffold
liana_scaffold <- liana_out %>%
    select(source, ligand, receptor, target, ligand.expr, receptor.expr)

# Seurat to SCE
entity_genes <- union(op_resource$source_genesymbol, op_resource$target_genesymbol)
sce <- seurat_to_sce(seurat_object, entity_genes = entity_genes, assay = "RNA")



# LIANA PIPE
pem_scores <- compute_pem_scores(sce = sce,
                                 assay.type = assay.type)

# ADD PEM to LIANA PIPE scaffold
lr_res <- liana_scaffold %>%
    join_means(means = pem_scores,
               source_target = "target",
               entity = "receptor",
               type = "pem") %>%
    join_means(means = pem_scores,
               source_target = "source",
               entity = "ligand",
               type = "pem")


# NST - in LIANA wrap
nst_scores <- compute_nst_scores(sce = sce,
                                 ligand_receptor_df = lr_res %>%
                                     select(ligand,receptor) %>%
                                     distinct(),
                                 assay.type = assay.type,
                                 seed = 1234)

liana_cts <- lr_res %>%
    left_join(nst_scores, by = c("target" = "celltype", "ligand", "receptor")) %>%
    dplyr::rename(target.nst = nst_score) %>%
    left_join(nst_scores, by = c("source" = "celltype", "ligand", "receptor")) %>%
    dplyr::rename(source.nst = nst_score) %>%
    group_by(source, target) %>%
    # compute cross-talk
    mutate(es = (ligand.pem + receptor.pem)/ 2, # calculate Expression Score
           nst = (source.nst + target.nst)/ 2 # combine -log10 NST scores
    ) %>%
    # normalize ES and NST
    mutate(Nes = minmax(es), Nnst = minmax(nst)) %>%
    mutate(cross_talk_score = if_else(source!=target,
                                      Nes * Nnst, # Para as CytoTalk
                                      Nes * (1-Nnst) # Autocrine we inverse NST
                                      )
           ) %>%
    ungroup() # %>%
    # filter(!(ligand.expr == 0)) %>%
    # filter(!(receptor.expr == 0))


