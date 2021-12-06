# LIANA: a LIgand-receptor ANalysis frAmework <img src="liana_logo.png" align="right" height="100">

The continuous developments of single-cell RNA-Seq (scRNA-Seq) have sparked
an immense interest in understanding intercellular crosstalk. Multiple
tools and resources that aid the investigation of cell-cell communication (CCC)
were published recently.
However, these methods and resources are usually in a fixed combination of a
tool and its corresponding resource, but in principle any resource could be
combined with any method.  

    
## Install LIANA  
```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
devtools::install_github('saezlab/liana')
```

If you wish to make use of the CellChat algorithm:
```{r}
devtools::install_github("sqjin/CellChat")
```

  
## Tutorial
See a [tutorial](https://saezlab.github.io/liana/articles/liana_tutorial.html) how to use LIANA to run up to 8 methods and 16 different resources!
The tutorial with the test data takes minutes to complete!    
  
Additional tutorials:  

* [How to combine LIANA and NicheNet](https://saezlab.github.io/liana/articles/liana_nichenet.html)  

* [How to make use of OmniPathR's intracellular component](https://saezlab.github.io/liana/articles/liana_intracell.html)  

* [How to customize OmniPath's intercellular resource](https://saezlab.github.io/liana/articles/liana_intracell.html)  

* [LIANA for developers and benchmarks](https://saezlab.github.io/liana/articles/liana_devel.html)  
  
  
## LIANA Framework
  
To this end, we built a framework to decouple the methods from their corresponding resources.
   
LIANA also goes a step further as it provides:

* A robust and extendable architecture that aims to accelerate method development and benchmarks

* A rank aggregate from the results of different methods

* A customizable plethora of resources
  
![landingpage](ligrec_pipe.png)
  
  
## Tools

The tools implemented in this repository are:

- [CellPhoneDBv2](https://github.com/Teichlab/cellphonedb) (*)
- [CellChat](https://github.com/sqjin/CellChat)
- [NATMI](https://github.com/forrest-lab/NATMI) (*)
- [Connectome](https://github.com/msraredon/Connectome) (`edge_weights`) (*)
- [SingleCellSignalR](https://github.com/SCA-IRCM/SingleCellSignalR) (`LRscores`) (SCA) (*)
- [iTALK](https://github.com/Coolgenome/iTALK)-inspired mean logFC score (`logfc`) (*)
- [CytoTalk](https://advances.sciencemag.org/content/7/16/eabf1356)-inspired `cross-talk` scores (*)
- [RobustRank](https://pubmed.ncbi.nlm.nih.gov/22247279/)-aggregate (`aggregate_rank`) scores (*)
  
*The scoring systems from these methods were re-implemented in LIANA in order to account for multimeric complexes, to simplify the calls to the individual pipelines, or reduce any possible inconsistencies and redundancies in their downstream integration. If you wish to run LIANA with the original tools please see [LIANA++](https://saezlab.github.io/liana/articles/liana_devel.html).
  
  
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
- CellTalkDB
- OmniPath
  
  
### OmniPath
  
All the resources above are retrieved from [OmniPath](https://omnipathdb.org/),
and more specifically [OmnipathR](https://github.com/saezlab/OmnipathR).
However, individual resources retrieved from the OmniPath web service are not to be affected by this, as each resource expected to be identical to its original form, apart from minor processing imperfections.
  
`OmniPath` itself serves as a composite CCC resource combining all the ones listed
above + [more](https://doi.org/10.15252/msb.20209923). `OmniPath` also collects
further information about the roles and localisation of proteins in intercellular communication.
We made use of this information regarding the and by default the `OmniPath`CCC
resource in LIANA is filtered according to the consensus localisation and curation of
ligand-receptor interactions. To obtain more information how we filtered the default CCC `OmniPath`,
as well as to explore custom filter options see [customizing OmniPath resources](https://saezlab.github.io/liana/articles/liana_custom_op.html).  
  
  
## LIANA++
If you are interested in making use of the LIANA architecture for your own method, [this vignette](https://saezlab.github.io/liana/articles/liana_devel.html) provides instructions how to obtain a comprehensive table of LR statistics, which can then be used by custom scoring functions.
In the [same vignette](https://saezlab.github.io/liana/articles/liana_devel.html) are also instructions how to install and run the original methods via a convenient R wrapper, e.g. for their unbiased benchmark.
  
  
  
## Contact
We appreciate any feedback, so please do not hesitate to open an issue on the [liana github page](https://github.com/saezlab/liana)!  
  
  
  
### Citing `LIANA`:
Dimitrov, D., Türei, D., Boys, C., Nagai, J.S., Flores, R.O.R., Kim, H., Szalai, B., Costa, I.G., Dugourd, A., Valdeolivas, A. and Saez-Rodriguez, J., 2021.  Comparison of Resources and Methods to infer Cell-Cell Communication from Single-cell RNA Data.
  bioRxiv. [10.1101/2021.05.21.445160v1](https://www.biorxiv.org/content/10.1101/2021.05.21.445160v1)

#### Also, if you use the `OmniPath` CCC Resource for your analysis, please cite:
Türei, D., Valdeolivas, A., Gul, L., Palacio‐Escat, N., Klein, M., Ivanova, O., Ölbei, M., Gábor, A., Theis, F., Módos, D. and Korcsmáros, T., 2021. Integrated intra‐and intercellular signaling knowledge for multicellular omics analysis. Molecular systems biology, 17(3), p.e9923.
https://doi.org/10.15252/msb.20209923
  
#### Similarly, please consider appropritaly citing any of the methods and/or resources implemented in liana, that were particularly relevant for your research!
