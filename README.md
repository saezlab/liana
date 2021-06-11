# LIANA: a LIgand-receptor ANalysis frAmework
  
![landingpage](ligrec_pipe.png)

## Objectives
  
The continuous developments of single-cell RNA-Seq (scRNA-Seq) have sparked
an immense interest in understanding intercellular crosstalk. Multiple
tools and resources that aid the investigation of cell-cell communication (CCC)
were published recently.
However, these methods and resources are usually in a fixed combination of a
tool and its corresponding resource, but in principle any resource could be
combined with any method.

## LIANA Framework
  
To this end we built a framework to decouple the tools from their corresponding resources.
In turn, we used the pipeline for the systematic comparison of all combinations between 15 resources and 6 tools.
  
LIANA is open access and available in this repository.
  
  

## Tools

The tools implemented in this repository are:

- CellPhoneDB algorithm (via [Squidpy](https://squidpy.readthedocs.io/en/latest/))
- CellChat
- NATMI
- Connectome
- SingleCellSignalR (SCA)
- iTALK
  
  
  
## Resources

### Cell-cell Communication resources

The following CCC resources are accessible via this pipeline:

- CellChatDB
- CellPhoneDB
- Ramilowski2015
- Baccin2019
- LRdb
- Kiroauc2010
- ICELLNET
- iTALK
- EMBRACE
- HPMR
- Guide2Pharma
- connectomeDB2020
- talklr
- CellTalkDB
- OmniPath
  
  
### OmniPath
  
All the resources above are retrieved from OmniPath (https://omnipathdb.org/).
OmniPath itself is also a composite resource combining all the ones listed
above. However the cell-cell interactions in OmniPath are more than simply
the union of the ligand-receptor resources. OmniPath uses several further
databases to collect information about the roles of proteins in intercellular
communication and again other databases to find connections between them. At
the same time, OmniPath blacklists certain wrong annotations, removing some
of the contents of the original resources. However the data of individual
resources retrieved from the OmniPath web service is not affected by this,
each resource supposed to be identical to its original form, apart from minor
processing imperfections.
  
  
### Reshuffled and Default Resources
  
Moreover, a Reshuffled resource can be generated via reshuffling any of the
aforementioned using the `BiRewire` package, and each tool can be run with
its 'Default' resource, the dataset used in its original publication.
  
  
  
## Dependancies
Please check the [.yml](https://github.com/saezlab/Cell_Cell_Investigation/blob/biorxiv/intercell.yml) file
(Recommended to use it to set up a conda environment) and then also install the following in R:
   
   
```{r}
library(devtools)  
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE) # ignore warning from iTALK 
install_github('saezlab/OmnipathR@ff3ad88e3915747e1b557bf44ac5396f9525dd7e') # install 4.0 version of OmnipathR
devtools::install_github("sqjin/CellChat")  
devtools::install_github('msraredon/Connectome', ref = 'master')   
devtools::install_github('Coolgenome/iTALK', ref = 'biorxiv')   
# A modified version of SingleCellSignalR (SCA) that enables external resources
devtools::install_github(repo = "CostaLab/SingleCellSignalR_v1", subdir = "SingleCellSignalR")
# Install this repos
devtools::install_github('saezlab/liana')
```
  
  
#### Clone modified NATMI into liana path, as it is run via the command line 
##### NATMI was forked from its [github repo](https://github.com/asrhou/NATMI.git) and changed to be agnostic of the working directory when loading the resources:
```{sh}
cd *insert fullpath*/liana
git clone https://github.com/saezlab/NATMI
```
  
### Citing us:
Please cite our recent [BioRxiv](https://www.biorxiv.org/content/10.1101/2021.05.21.445160v1) manuscript, if you use LIANA.


## References
Baccin, C., Al-Sabah, J., Velten, L., Helbling, P.M., Grünschläger, F., Hernández-Malmierca, P., Nombela-Arrieta, C., Steinmetz, L.M., Trumpp, A., and Haas, S. (2020). Combined single-cell and spatial transcriptomics reveal the molecular, cellular and spatial bone marrow niche organization. Nat. Cell Biol. 22, 38–48.

Ben-Shlomo, I., Yu Hsu, S., Rauch, R., Kowalski, H.W., and Hsueh, A.J.W. (2003). Signaling receptome: a genomic and evolutionary perspective of plasma membrane receptors involved in signal transduction. Sci. STKE 2003, RE9.

Cabello-Aguilar, S., Alame, M., Kon-Sun-Tack, F., Fau, C., Lacroix, M., and Colinge, J. (2020). SingleCellSignalR: inference of intercellular networks from single-cell transcriptomics. Nucleic Acids Res. 48, e55.

Efremova, M., Vento-Tormo, M., Teichmann, S.A., and Vento-Tormo, R. (2020). CellPhoneDB: inferring cell-cell communication from combined expression of multi-subunit ligand-receptor complexes. Nat. Protoc. 15, 1484–1506.

Harding, S.D., Sharman, J.L., Faccenda, E., Southan, C., Pawson, A.J., Ireland, S., Gray, A.J.G., Bruce, L., Alexander, S.P.H., Anderton, S., et al. (2018). The IUPHAR/BPS Guide to PHARMACOLOGY in 2018: updates and expansion to encompass the new guide to IMMUNOPHARMACOLOGY. Nucleic Acids Res. 46, D1091–D1106.

Hou, R., Denisenko, E., Ong, H.T., Ramilowski, J.A., and Forrest, A.R.R. (2020). Predicting cell-to-cell communication networks using NATMI. Nat. Commun. 11, 5011.

Jin, S., Guerrero-Juarez, C.F., Zhang, L., Chang, I., Ramos, R., Kuan, C.-H., Myung, P., Plikus, M.V., and Nie, Q. (2021). Inference and analysis of cell-cell communication using CellChat. Nat. Commun. 12, 1088.


Noël, F., Massenet-Regad, L., Carmi-Levy, I., Cappuccio, A., Grandclaudon, M., Trichot, C., Kieffer, Y., Mechta-Grigoriou, F., and Soumelis, V. (2021). Dissection of intercellular communication using the transcriptome-based framework ICELLNET. Nat. Commun. 12, 1089.

Palla, G., Spitzer, H., Klein, M., Fischer, D.S., Schaar, A.C., Kuemmerle, L.B., Rybakov, S., Ibarra, I.L., Holmberg, O., Virshup, I., et al. (2021). Squidpy: a scalable framework for spatial single cell analysis. BioRxiv.

Ramilowski, J.A., Goldberg, T., Harshbarger, J., Kloppmann, E., Lizio, M., Satagopam, V.P., Itoh, M., Kawaji, H., Carninci, P., Rost, B., et al. (2015). A draft network of ligand-receptor-mediated multicellular signalling in human. Nat. Commun. 6, 7866.

Raredon, M.S.B., Yang, J., Garritano, J., Wang, M., Kushnir, D., Schupp, J.C., Adams, T.S., Greaney, A.M., Leiby, K.L., Kaminski, N., et al. (2021). Connectome: computation and visualization of cell-cell signaling topologies in single-cell systems data. BioRxiv.

Shao, X., Liao, J., Li, C., Lu, X., Cheng, J., and Fan, X. (2020). CellTalkDB: a manually curated database of ligand-receptor interactions in humans and mice. Brief. Bioinformatics.

Sheikh, B.N., Bondareva, O., Guhathakurta, S., Tsang, T.H., Sikora, K., Aizarani, N., Sagar, Holz, H., Grün, D., Hein, L., et al. (2019). Systematic Identification of Cell-Cell Communication Networks in the Developing Brain. IScience 21, 273–287.

Türei, D., Valdeolivas, A., Gul, L., Palacio-Escat, N., Klein, M., Ivanova, O., Ölbei, M., Gábor, A., Theis, F., Módos, D., et al. (2021). Integrated intra- and intercellular signaling knowledge for multicellular omics analysis. Mol. Syst. Biol. 17, e9923.

Vento-Tormo, R., Efremova, M., Botting, R.A., Turco, M.Y., Vento-Tormo, M., Meyer, K.B., Park, J.-E., Stephenson, E., Polański, K., Goncalves, A., et al. (2018). Single-cell reconstruction of the early maternal-fetal interface in humans. Nature 563, 347–353.

Wang, Y. (2020). talklr uncovers ligand-receptor mediated intercellular crosstalk. BioRxiv.

Wang, Y., Wang, R., Zhang, S., Song, S., Jiang, C., Han, G., Wang, M., Ajani, J., Futreal, A., and Wang, L. (2019). iTALK: an R Package to Characterize and Illustrate Intercellular Communication. BioRxiv.


  
### R sessionInfo()
```r
─ Session info ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
 setting  value                       
 version  R version 4.0.3 (2020-10-10)
 os       Ubuntu 20.04.2 LTS          
 system   x86_64, linux-gnu           
 ui       RStudio                     
 language (EN)                        
 collate  en_US.UTF-8                 
 ctype    en_US.UTF-8                 
 tz       Europe/Berlin               
 date     2021-06-11                  

─ Packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
 package              * version    date       lib source                               
 abind                  1.4-5      2016-07-21 [1] CRAN (R 4.0.2)                       
 annotate               1.68.0     2020-10-27 [1] Bioconductor                         
 AnnotationDbi          1.52.0     2020-10-27 [1] Bioconductor                         
 assertthat             0.2.1      2019-03-21 [1] CRAN (R 4.0.2)                       
 backports              1.2.1      2020-12-09 [1] CRAN (R 4.0.2)                       
 bbmle                  1.0.23.1   2020-02-03 [1] CRAN (R 4.0.3)                       
 bdsmatrix              1.3-4      2020-01-13 [1] CRAN (R 4.0.3)                       
 Biobase              * 2.50.0     2020-10-27 [1] Bioconductor                         
 BiocGenerics         * 0.36.0     2020-10-27 [1] Bioconductor                         
 BiocParallel           1.24.1     2020-11-06 [1] Bioconductor                         
 bit                    4.0.4      2020-08-04 [1] CRAN (R 4.0.2)                       
 bit64                  4.0.5      2020-08-30 [1] CRAN (R 4.0.2)                       
 bitops                 1.0-6      2013-08-17 [1] CRAN (R 4.0.2)                       
 blob                   1.2.1      2020-01-20 [1] CRAN (R 4.0.2)                       
 brew                   1.0-6      2011-04-13 [1] CRAN (R 4.0.2)                       
 cachem                 1.0.4      2021-02-13 [1] CRAN (R 4.0.3)                       
 Cairo                  1.5-12.2   2020-07-07 [1] CRAN (R 4.0.3)                       
 caTools                1.18.1     2021-01-11 [1] CRAN (R 4.0.3)                       
 CellChat               1.1.0      2021-05-02 [1] Github (sqjin/CellChat@b3ccf96)      
 cellranger             1.1.0      2016-07-27 [1] CRAN (R 4.0.2)                       
 checkmate              2.0.0      2020-02-06 [1] CRAN (R 4.0.3)                       
 circlize               0.4.12     2021-01-08 [1] CRAN (R 4.0.3)                       
 cli                    2.3.0      2021-01-31 [1] CRAN (R 4.0.3)                       
 clipr                  0.7.1      2020-10-08 [1] CRAN (R 4.0.2)                       
 clue                   0.3-58     2020-12-03 [1] CRAN (R 4.0.3)                       
 cluster                2.1.0      2019-06-19 [4] CRAN (R 4.0.0)                       
 coda                   0.19-4     2020-09-30 [1] CRAN (R 4.0.3)                       
 codetools              0.2-16     2018-12-24 [4] CRAN (R 4.0.0)                       
 colorspace             2.0-0      2020-11-11 [1] CRAN (R 4.0.2)                       
 combinat               0.0-8      2012-10-29 [1] CRAN (R 4.0.2)                       
 ComplexHeatmap         2.6.2      2020-11-12 [1] Bioconductor                         
 Connectome             1.0.1      2021-02-21 [1] Github (msraredon/Connectome@4531dd9)
 conquer                1.0.2      2020-08-27 [1] CRAN (R 4.0.2)                       
 cowplot                1.1.1      2020-12-30 [1] CRAN (R 4.0.3)                       
 crayon                 1.4.1      2021-02-08 [1] CRAN (R 4.0.3)                       
 curl                   4.3        2019-12-02 [1] CRAN (R 4.0.2)                       
 data.table             1.14.0     2021-02-21 [1] CRAN (R 4.0.3)                       
 DBI                    1.1.1      2021-01-15 [1] CRAN (R 4.0.3)                       
 DDRTree                0.1.5      2017-04-30 [1] CRAN (R 4.0.2)                       
 DelayedArray           0.16.1     2021-01-22 [1] Bioconductor                         
 deldir                 0.2-10     2021-02-16 [1] CRAN (R 4.0.3)                       
 densityClust           0.3        2017-10-24 [1] CRAN (R 4.0.2)                       
 desc                   1.2.0      2018-05-01 [1] CRAN (R 4.0.2)                       
 DESeq2                 1.30.1     2021-02-19 [1] Bioconductor                         
 DEsingle               1.10.0     2020-10-27 [1] Bioconductor                         
 digest                 0.6.27     2020-10-24 [1] CRAN (R 4.0.2)                       
 distillery             1.2        2020-11-21 [1] CRAN (R 4.0.3)                       
 docopt                 0.7.1      2020-06-24 [1] CRAN (R 4.0.2)                       
 doParallel             1.0.16     2020-10-16 [1] CRAN (R 4.0.2)                       
 dplyr                  1.0.4      2021-02-02 [1] CRAN (R 4.0.3)                       
 edgeR                  3.32.1     2021-01-14 [1] Bioconductor                         
 ellipsis               0.3.1      2020-05-15 [1] CRAN (R 4.0.2)                       
 evaluate               0.14       2019-05-28 [1] CRAN (R 4.0.2)                       
 extRemes               2.1        2020-11-24 [1] CRAN (R 4.0.3)                       
 fastICA                1.2-2      2019-07-08 [1] CRAN (R 4.0.2)                       
 fastmap                1.1.0      2021-01-25 [1] CRAN (R 4.0.3)                       
 fitdistrplus           1.1-3      2020-12-05 [1] CRAN (R 4.0.2)                       
 flexmix                2.3-17     2020-10-12 [1] CRAN (R 4.0.3)                       
 FNN                    1.1.3      2019-02-15 [1] CRAN (R 4.0.2)                       
 foreach                1.5.1      2020-10-15 [1] CRAN (R 4.0.2)                       
 future                 1.21.0     2020-12-10 [1] CRAN (R 4.0.2)                       
 future.apply           1.7.0      2021-01-04 [1] CRAN (R 4.0.3)                       
 gamlss                 5.2-0      2020-09-12 [1] CRAN (R 4.0.3)                       
 gamlss.data            5.1-4      2019-05-15 [1] CRAN (R 4.0.3)                       
 gamlss.dist            5.1-7      2020-07-13 [1] CRAN (R 4.0.3)                       
 genefilter             1.72.1     2021-01-21 [1] Bioconductor                         
 geneplotter            1.68.0     2020-10-27 [1] Bioconductor                         
 generics               0.1.0      2020-10-31 [1] CRAN (R 4.0.2)                       
 GenomeInfoDb           1.26.2     2020-12-08 [1] Bioconductor                         
 GenomeInfoDbData       1.2.4      2020-12-14 [1] Bioconductor                         
 GenomicRanges          1.42.0     2020-10-27 [1] Bioconductor                         
 GetoptLong             1.0.5      2020-12-15 [1] CRAN (R 4.0.3)                       
 gg.gap                 1.3        2019-09-30 [1] CRAN (R 4.0.3)                       
 ggalluvial             0.12.3     2020-12-05 [1] CRAN (R 4.0.3)                       
 ggplot2                3.3.3      2020-12-30 [1] CRAN (R 4.0.3)                       
 ggrepel                0.9.1      2021-01-15 [1] CRAN (R 4.0.3)                       
 ggridges               0.5.3      2021-01-08 [1] CRAN (R 4.0.3)                       
 GlobalOptions          0.1.2      2020-06-10 [1] CRAN (R 4.0.3)                       
 globals                0.14.0     2020-11-22 [1] CRAN (R 4.0.2)                       
 glue                   1.4.2      2020-08-27 [1] CRAN (R 4.0.2)                       
 goftest                1.2-2      2019-12-02 [1] CRAN (R 4.0.2)                       
 gplots                 3.1.1      2020-11-28 [1] CRAN (R 4.0.2)                       
 gridBase               0.4-7      2014-02-24 [1] CRAN (R 4.0.3)                       
 gridExtra              2.3        2017-09-09 [1] CRAN (R 4.0.2)                       
 gtable                 0.3.0      2019-03-25 [1] CRAN (R 4.0.2)                       
 gtools                 3.8.2      2020-03-31 [1] CRAN (R 4.0.2)                       
 hms                    1.0.0      2021-01-13 [1] CRAN (R 4.0.3)                       
 HSMMSingleCell         1.10.0     2020-10-29 [1] Bioconductor                         
 htmltools              0.5.1.1    2021-01-22 [1] CRAN (R 4.0.3)                       
 htmlwidgets            1.5.3      2020-12-10 [1] CRAN (R 4.0.2)                       
 httpuv                 1.5.5      2021-01-13 [1] CRAN (R 4.0.3)                       
 httr                   1.4.2      2020-07-20 [1] CRAN (R 4.0.2)                       
 ica                    1.0-2      2018-05-24 [1] CRAN (R 4.0.2)                       
 igraph                 1.2.6      2020-10-06 [1] CRAN (R 4.0.2)                       
 IRanges                2.24.1     2020-12-12 [1] Bioconductor                         
 irlba                  2.3.3      2019-02-05 [1] CRAN (R 4.0.2)                       
 iTALK                  0.1.0      2021-02-21 [1] Github (Coolgenome/iTALK@6d9b390)    
 iterators              1.0.13     2020-10-15 [1] CRAN (R 4.0.2)                       
 jsonlite               1.7.2      2020-12-09 [1] CRAN (R 4.0.2)                       
 KernSmooth             2.23-17    2020-04-26 [4] CRAN (R 4.0.0)                       
 knitr                  1.31       2021-01-27 [1] CRAN (R 4.0.3)                       
 later                  1.1.0.1    2020-06-05 [1] CRAN (R 4.0.2)                       
 lattice                0.20-41    2020-04-02 [4] CRAN (R 4.0.0)                       
 lazyeval               0.2.2      2019-03-15 [1] CRAN (R 4.0.2)                       
 leiden                 0.3.7      2021-01-26 [1] CRAN (R 4.0.3)                       
 liana                * 0.0.1      2021-06-11 [1] local                                
 lifecycle              1.0.0      2021-02-15 [1] CRAN (R 4.0.3)                       
 limma                  3.46.0     2020-10-27 [1] Bioconductor                         
 listenv                0.8.0      2019-12-05 [1] CRAN (R 4.0.2)                       
 Lmoments               1.3-1      2019-03-15 [1] CRAN (R 4.0.3)                       
 lmtest                 0.9-38     2020-09-09 [1] CRAN (R 4.0.2)                       
 locfit                 1.5-9.4    2020-03-25 [1] CRAN (R 4.0.2)                       
 logger                 0.2.0      2021-03-04 [1] CRAN (R 4.0.3)                       
 magrittr               2.0.1      2020-11-17 [1] CRAN (R 4.0.2)                       
 MASS                   7.3-53     2020-09-09 [4] CRAN (R 4.0.2)                       
 MAST                   1.16.0     2020-10-27 [1] Bioconductor                         
 Matrix                 1.2-18     2019-11-27 [4] CRAN (R 4.0.0)                       
 MatrixGenerics         1.2.1      2021-01-30 [1] Bioconductor                         
 MatrixModels           0.4-1      2015-08-22 [1] CRAN (R 4.0.2)                       
 matrixStats            0.58.0     2021-01-29 [1] CRAN (R 4.0.3)                       
 maxLik                 1.4-6      2020-11-24 [1] CRAN (R 4.0.3)                       
 memoise                2.0.0      2021-01-26 [1] CRAN (R 4.0.3)                       
 mgcv                   1.8-33     2020-08-27 [4] CRAN (R 4.0.2)                       
 mime                   0.10       2021-02-13 [1] CRAN (R 4.0.3)                       
 miniUI                 0.1.1.1    2018-05-18 [1] CRAN (R 4.0.2)                       
 miscTools              0.6-26     2019-12-08 [1] CRAN (R 4.0.3)                       
 modeltools             0.2-23     2020-03-05 [1] CRAN (R 4.0.3)                       
 monocle                2.18.0     2020-10-27 [1] Bioconductor                         
 multtest               2.46.0     2020-10-27 [1] Bioconductor                         
 munsell                0.5.0      2018-06-12 [1] CRAN (R 4.0.2)                       
 mvtnorm                1.1-1      2020-06-09 [1] CRAN (R 4.0.2)                       
 network                1.16.1     2020-10-07 [1] CRAN (R 4.0.2)                       
 nlme                   3.1-149    2020-08-23 [4] CRAN (R 4.0.2)                       
 NMF                    0.23.0     2020-08-01 [1] CRAN (R 4.0.3)                       
 nnet                   7.3-14     2020-04-26 [4] CRAN (R 4.0.0)                       
 numDeriv               2016.8-1.1 2019-06-06 [1] CRAN (R 4.0.3)                       
 OmnipathR              2.99.18    2021-06-09 [1] Github (saezlab/OmnipathR@ff3ad88)   
 parallelly             1.24.0     2021-03-14 [1] CRAN (R 4.0.3)                       
 patchwork              1.1.1      2020-12-17 [1] CRAN (R 4.0.3)                       
 pbapply                1.4-3      2020-08-18 [1] CRAN (R 4.0.2)                       
 pcaMethods             1.82.0     2020-10-27 [1] Bioconductor                         
 pheatmap               1.0.12     2019-01-04 [1] CRAN (R 4.0.2)                       
 pillar                 1.4.7      2020-11-20 [1] CRAN (R 4.0.2)                       
 pkgconfig              2.0.3      2019-09-22 [1] CRAN (R 4.0.2)                       
 pkgload                1.1.0      2020-05-29 [1] CRAN (R 4.0.2)                       
 pkgmaker               0.32.2     2020-10-20 [1] CRAN (R 4.0.3)                       
 plotly                 4.9.3      2021-01-10 [1] CRAN (R 4.0.3)                       
 plyr                   1.8.6      2020-03-03 [1] CRAN (R 4.0.2)                       
 png                    0.1-7      2013-12-03 [1] CRAN (R 4.0.2)                       
 polyclip               1.10-0     2019-03-14 [1] CRAN (R 4.0.2)                       
 pracma                 2.3.3      2021-01-23 [1] CRAN (R 4.0.3)                       
 prettyunits            1.1.1      2020-01-24 [1] CRAN (R 4.0.2)                       
 progress               1.2.2      2019-05-16 [1] CRAN (R 4.0.2)                       
 promises               1.2.0.1    2021-02-11 [1] CRAN (R 4.0.3)                       
 pscl                   1.5.5      2020-03-07 [1] CRAN (R 4.0.3)                       
 purrr                * 0.3.4      2020-04-17 [1] CRAN (R 4.0.2)                       
 qlcMatrix              0.9.7      2018-04-20 [1] CRAN (R 4.0.2)                       
 quantreg               5.83       2021-01-22 [1] CRAN (R 4.0.3)                       
 R6                     2.5.0      2020-10-28 [1] CRAN (R 4.0.2)                       
 randomcoloR            1.1.0.1    2019-11-24 [1] CRAN (R 4.0.3)                       
 RANN                   2.6.1      2019-01-08 [1] CRAN (R 4.0.2)                       
 rappdirs               0.3.3      2021-01-31 [1] CRAN (R 4.0.3)                       
 RColorBrewer           1.1-2      2014-12-07 [1] CRAN (R 4.0.2)                       
 Rcpp                   1.0.6      2021-01-15 [1] CRAN (R 4.0.3)                       
 RcppAnnoy              0.0.18     2020-12-15 [1] CRAN (R 4.0.3)                       
 RcppArmadillo          0.10.2.1.0 2021-02-09 [1] CRAN (R 4.0.3)                       
 RCurl                  1.98-1.3   2021-03-16 [1] CRAN (R 4.0.3)                       
 readr                  1.4.0      2020-10-05 [1] CRAN (R 4.0.2)                       
 readxl                 1.3.1      2019-03-13 [1] CRAN (R 4.0.2)                       
 registry               0.5-1      2019-03-05 [1] CRAN (R 4.0.3)                       
 reshape2               1.4.4      2020-04-09 [1] CRAN (R 4.0.2)                       
 reticulate             1.18       2020-10-25 [1] CRAN (R 4.0.2)                       
 rjson                  0.2.20     2018-06-08 [1] CRAN (R 4.0.3)                       
 rlang                  0.4.10     2020-12-30 [1] CRAN (R 4.0.3)                       
 rle                    0.9.2      2020-09-25 [1] CRAN (R 4.0.3)                       
 rmarkdown              2.7        2021-02-19 [1] CRAN (R 4.0.3)                       
 RMTstat                0.3        2014-11-01 [1] CRAN (R 4.0.3)                       
 rngtools               1.5        2020-01-23 [1] CRAN (R 4.0.2)                       
 ROCR                   1.0-11     2020-05-02 [1] CRAN (R 4.0.2)                       
 Rook                   1.1-1      2014-10-20 [1] CRAN (R 4.0.3)                       
 rpart                  4.1-15     2019-04-12 [4] CRAN (R 4.0.0)                       
 rprojroot              2.0.2      2020-11-15 [1] CRAN (R 4.0.2)                       
 RSpectra               0.16-0     2019-12-01 [1] CRAN (R 4.0.2)                       
 RSQLite                2.2.3      2021-01-24 [1] CRAN (R 4.0.3)                       
 rsvd                   1.0.3      2020-02-17 [1] CRAN (R 4.0.2)                       
 Rtsne                  0.15       2018-11-10 [1] CRAN (R 4.0.2)                       
 S4Vectors              0.28.1     2020-12-09 [1] Bioconductor                         
 sandwich               3.0-0      2020-10-02 [1] CRAN (R 4.0.2)                       
 scales                 1.1.1      2020-05-11 [1] CRAN (R 4.0.2)                       
 SCAomni                0.0.1.8    2021-02-21 [1] Bioconductor                         
 scattermore            0.7        2020-11-24 [1] CRAN (R 4.0.3)                       
 scde                   2.18.0     2020-10-27 [1] Bioconductor                         
 sctransform            0.3.2      2020-12-16 [1] CRAN (R 4.0.3)                       
 sessioninfo            1.1.1      2018-11-05 [1] CRAN (R 4.0.2)                       
 Seurat                 3.2.3      2020-12-15 [1] CRAN (R 4.0.3)                       
 shape                  1.4.5      2020-09-13 [1] CRAN (R 4.0.3)                       
 shiny                  1.6.0      2021-01-25 [1] CRAN (R 4.0.3)                       
 SIMLR                  1.16.0     2020-10-27 [1] Bioconductor                         
 SingleCellExperiment   1.12.0     2020-10-27 [1] Bioconductor                         
 slam                   0.1-48     2020-12-03 [1] CRAN (R 4.0.2)                       
 sna                    2.6        2020-10-06 [1] CRAN (R 4.0.3)                       
 SparseM                1.81       2021-02-18 [1] CRAN (R 4.0.3)                       
 sparsesvd              0.2        2019-07-15 [1] CRAN (R 4.0.2)                       
 spatstat               1.64-1     2020-05-12 [1] CRAN (R 4.0.2)                       
 spatstat.data          2.1-0      2021-03-21 [1] CRAN (R 4.0.3)                       
 spatstat.utils         2.1-0      2021-03-15 [1] CRAN (R 4.0.3)                       
 statnet.common         4.4.1      2020-10-03 [1] CRAN (R 4.0.3)                       
 stringi                1.5.3      2020-09-09 [1] CRAN (R 4.0.2)                       
 stringr                1.4.0      2019-02-10 [1] CRAN (R 4.0.2)                       
 SummarizedExperiment   1.20.0     2020-10-27 [1] Bioconductor                         
 survival               3.2-7      2020-09-28 [1] CRAN (R 4.0.2)                       
 svglite                2.0.0      2021-02-20 [1] CRAN (R 4.0.3)                       
 systemfonts            1.0.1      2021-02-09 [1] CRAN (R 4.0.3)                       
 tensor                 1.5        2012-05-05 [1] CRAN (R 4.0.2)                       
 testthat             * 3.0.2      2021-02-14 [1] CRAN (R 4.0.3)                       
 tibble               * 3.0.6      2021-01-29 [1] CRAN (R 4.0.3)                       
 tidyr                  1.1.2      2020-08-27 [1] CRAN (R 4.0.2)                       
 tidyselect             1.1.0      2020-05-11 [1] CRAN (R 4.0.2)                       
 uwot                   0.1.10     2020-12-15 [1] CRAN (R 4.0.3)                       
 V8                     3.4.0      2020-11-04 [1] CRAN (R 4.0.3)                       
 vctrs                  0.3.6      2020-12-17 [1] CRAN (R 4.0.3)                       
 VGAM                   1.1-5      2021-01-14 [1] CRAN (R 4.0.3)                       
 viridis                0.5.1      2018-03-29 [1] CRAN (R 4.0.2)                       
 viridisLite            0.3.0      2018-02-01 [1] CRAN (R 4.0.2)                       
 withr                  2.4.1      2021-01-26 [1] CRAN (R 4.0.3)                       
 xfun                   0.21       2021-02-10 [1] CRAN (R 4.0.3)                       
 XML                    3.99-0.5   2020-07-23 [1] CRAN (R 4.0.2)                       
 xml2                   1.3.2      2020-04-23 [1] CRAN (R 4.0.2)                       
 xtable                 1.8-4      2019-04-21 [1] CRAN (R 4.0.2)                       
 XVector                0.30.0     2020-10-27 [1] Bioconductor                         
 yaml                   2.2.1      2020-02-01 [1] CRAN (R 4.0.2)                       
 zlibbioc               1.36.0     2020-10-27 [1] Bioconductor                         
 zoo                    1.8-8      2020-05-02 [1] CRAN (R 4.0.2)                       

[1] /home/dbdimitrov/R/x86_64-pc-linux-gnu-library/4.0
[2] /usr/local/lib/R/site-library
[3] /usr/lib/R/site-library
[4] /usr/lib/R/library
```

