# Cell-Cell Communication Investigation

## Goal

The continuous developments of single-cell RNA-Seq (scRNA-Seq) have sparked an immense interest in understanding multi-cellular crosstalk. Multiple methods and resources that aid the investigation of cell-cell communication have been recently published.
However, these methods and resources are usually in a fixed combination of a method and its corresponding resource, but in principle any resource could be combined with any statistical method. Yet, it is largely unclear the difference that the choice of resource and method can have on the predicted CCC events. Thus, we attempt to shed some light on this topic via a systematic overview of how different combinations might influence CCC inference, by decoupling the methods from their corresponding resources.

As such, we compare all combinations between 15 resources and 7 methods and explore the effect on downstream results.

We provide the resources and methods as a pipeline for further use in this repository.


## Methods

The methods implemented in this repository are:

CellPhoneDB algorithm (via Squidpy);
CellChat;
NATMI;
Connectome;
SingleCellSignalR (SCAomni);
iTALK;
scTalk - to be implemented;
cell2cell - to be implemented.

## Resources

The intercellular signalling resources were queried from OmniPath and are the following:
CellChatDB;
CellPhoneDB;
Ramilowski2015;
Baccin2019;
LRdb;
Kiroauc2010;
ICELLNET;
iTALK;
EMBRACE;
HPMR;
Guide2Pharma;
connectomeDB2020;
talklr;
CellTalkDB;
OmniPathDB - a composite resource from all resources.

Moreover, a Randomized resource can be generated via reshuffling any of the abovementioned using BiRewire, and each tool can be run with its 'Default' inbuilt resource.



## Dependancies
Please check the .yml file (Recommended to use it to set up a conda environment) and the also install the following in R:

library(devtools)
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE) # iTALK throws a warning...
devtools::install_github("sqjin/CellChat")
devtools::install_github('msraredon/Connectome', ref = 'master')
devtools::install_github('Coolgenome/iTALK', ref = 'biorxiv')
devtools::install_github('saezlab/Cell_Cell_Investigation', ref = 'biorxiv')

#### Modified version of SingleCellSignalR (SCA)
devtools::install_github(repo = "https://github.com/CostaLab/SingleCellSignalR_v1", subdir = "SingleCellSignalR")    

#### Clone NATMI, as it is ran natively (i.e. in cmd)
git clone https://github.com/asrhou/NATMI.git

### R sessionInfo()
R version 4.0.3 (2020-10-10)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.2 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=de_DE.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=de_DE.UTF-8   
 [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=de_DE.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] intercell_0.0.1     ggfortify_0.4.11    RColorBrewer_1.1-2  UpSetR_1.4.0        magrittr_2.0.1      yardstick_0.0.7     BiRewire_3.22.0    
 [8] Matrix_1.2-18       tsne_0.1-3          slam_0.1-48         reticulate_1.18     rprojroot_2.0.2     SCAomni_0.0.1.8     Seurat_3.2.3       
[15] iTALK_0.1.0         CellChat_0.5.5      Connectome_1.0.1    OmnipathR_2.0.0     jsonlite_1.7.2      igraph_1.2.6        forcats_0.5.1      
[22] stringr_1.4.0       dplyr_1.0.4         purrr_0.3.4         readr_1.4.0         tidyr_1.1.2         tibble_3.0.6        ggplot2_3.3.3      
[29] tidyverse_1.3.0     Biobase_2.50.0      BiocGenerics_0.36.0

loaded via a namespace (and not attached):
  [1] statnet.common_4.4.1        rsvd_1.0.3                  ica_1.0-2                   svglite_2.0.0               ps_1.5.0                   
  [6] foreach_1.5.1               lmtest_0.9-38               crayon_1.4.1                V8_3.4.0                    MASS_7.3-53                
 [11] MAST_1.16.0                 nlme_3.1-149                backports_1.2.1             qlcMatrix_0.9.7             reprex_1.0.0               
 [16] rlang_0.4.10                XVector_0.30.0              ROCR_1.0-11                 readxl_1.3.1                irlba_2.3.3                
 [21] SparseM_1.81                callr_3.5.1                 limma_3.46.0                BiocParallel_1.24.1         rjson_0.2.20               
 [26] bit64_4.0.5                 glue_1.4.2                  pheatmap_1.0.12             rngtools_1.5                sctransform_0.3.2          
 [31] processx_3.4.5              AnnotationDbi_1.52.0        VGAM_1.1-5                  haven_2.3.1                 tidyselect_1.1.0           
 [36] SummarizedExperiment_1.20.0 usethis_2.0.1               fitdistrplus_1.1-3          XML_3.99-0.5                DEsingle_1.10.0            
 [41] zoo_1.8-8                   gg.gap_1.3                  xtable_1.8-4                MatrixModels_0.4-1          cli_2.3.0                  
 [46] zlibbioc_1.36.0             rstudioapi_0.13             miniUI_0.1.1.1              rpart_4.1-15                shiny_1.6.0                
 [51] clue_0.3-58                 multtest_2.46.0             pkgbuild_1.2.0              cluster_2.1.0               caTools_1.18.1             
 [56] pcaMethods_1.82.0           quantreg_5.83               ggrepel_0.9.1               listenv_0.8.0               png_0.1-7                  
 [61] future_1.21.0               withr_2.4.1                 rle_0.9.2                   bitops_1.0-6                plyr_1.8.6                 
 [66] cellranger_1.1.0            sparsesvd_0.2               pROC_1.17.0.1               pracma_2.3.3                coda_0.19-4                
 [71] pillar_1.4.7                gplots_3.1.1                GlobalOptions_0.1.2         cachem_1.0.4                pscl_1.5.5                 
 [76] fs_1.5.0                    flexmix_2.3-17              GetoptLong_1.0.5            gamlss.data_5.1-4           vctrs_0.3.6                
 [81] ellipsis_0.3.1              generics_0.1.0              randomcoloR_1.1.0.1         devtools_2.3.2              NMF_0.23.0                 
 [86] tools_4.0.3                 munsell_0.5.0               distillery_1.2              DelayedArray_0.16.1         fastmap_1.1.0              
 [91] compiler_4.0.3              HSMMSingleCell_1.10.0       pkgload_1.1.0               abind_1.4-5                 httpuv_1.5.5               
 [96] extRemes_2.1                sessioninfo_1.1.1           pkgmaker_0.32.2             plotly_4.9.3                GenomeInfoDbData_1.2.4     
[101] gridExtra_2.3               edgeR_3.32.1                lattice_0.20-41             deldir_0.2-10               later_1.1.0.1              
[106] scales_1.1.1                docopt_0.7.1                pbapply_1.4-3               genefilter_1.72.1           lazyeval_0.2.2             
[111] promises_1.2.0.1            spatstat_1.64-1             doParallel_1.0.16           goftest_1.2-2               spatstat.utils_2.1-0       
[116] brew_1.0-6                  sna_2.6                     sandwich_3.0-0              cowplot_1.1.1               Rtsne_0.15                 
[121] uwot_0.1.10                 Rook_1.1-1                  survival_3.2-7              numDeriv_2016.8-1.1         systemfonts_1.0.1          
[126] DDRTree_0.1.5               htmltools_0.5.1.1           memoise_2.0.0               modeltools_0.2-23           locfit_1.5-9.4             
[131] IRanges_2.24.1              viridisLite_0.3.0           digest_0.6.27               assertthat_0.2.1            mime_0.10                  
[136] densityClust_0.3            registry_0.5-1              SIMLR_1.16.0                RSQLite_2.2.3               future.apply_1.7.0         
[141] remotes_2.2.0               RcppArmadillo_0.10.2.1.0    data.table_1.14.0           blob_1.2.1                  S4Vectors_0.28.1           
[146] fastICA_1.2-2               splines_4.0.3               Cairo_1.5-12.2              RCurl_1.98-1.2              broom_0.7.5                
[151] monocle_2.18.0              hms_1.0.0                   modelr_0.1.8                colorspace_2.0-0            GenomicRanges_1.42.0       
[156] shape_1.4.5                 maxLik_1.4-6                nnet_7.3-14                 Rcpp_1.0.6                  RANN_2.6.1                 
[161] mvtnorm_1.1-1               circlize_0.4.12             conquer_1.0.2               parallelly_1.23.0           R6_2.5.0                   
[166] grid_4.0.3                  ggridges_0.5.3              lifecycle_1.0.0             miscTools_0.6-26            curl_4.3                   
[171] leiden_0.3.7                testthat_3.0.2              desc_1.2.0                  scde_2.18.0                 RcppAnnoy_0.0.18           
[176] iterators_1.0.13            htmlwidgets_1.5.3           polyclip_1.10-0             RMTstat_0.3                 network_1.16.1             
[181] gamlss_5.2-0                rvest_0.3.6                 ComplexHeatmap_2.6.2        mgcv_1.8-33                 globals_0.14.0             
[186] patchwork_1.1.1             bdsmatrix_1.3-4             codetools_0.2-16            Lmoments_1.3-1              matrixStats_0.58.0         
[191] lubridate_1.7.9.2           FNN_1.1.3                   gtools_3.8.2                prettyunits_1.1.1           SingleCellExperiment_1.12.0
[196] dbplyr_2.1.0                gridBase_0.4-7              RSpectra_0.16-0             GenomeInfoDb_1.26.2         gtable_0.3.0               
[201] DBI_1.1.1                   ggalluvial_0.12.3           stats4_4.0.3                tensor_1.5                  httr_1.4.2                 
[206] KernSmooth_2.23-17          stringi_1.5.3               progress_1.2.2              reshape2_1.4.4              annotate_1.68.0            
[211] viridis_0.5.1               xml2_1.3.2                  combinat_0.0-8              bbmle_1.0.23.1              geneplotter_1.68.0         
[216] scattermore_0.7             DESeq2_1.30.1               bit_4.0.4                   MatrixGenerics_1.2.1        spatstat.data_2.1-0        
[221] pkgconfig_2.0.3             gamlss.dist_5.1-7  
