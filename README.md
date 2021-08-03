# LIANA: a LIgand-receptor ANalysis frAmework <img src="liana_logo.png" align="right" height="100">

The continuous developments of single-cell RNA-Seq (scRNA-Seq) have sparked
an immense interest in understanding intercellular crosstalk. Multiple
tools and resources that aid the investigation of cell-cell communication (CCC)
were published recently.
However, these methods and resources are usually in a fixed combination of a
tool and its corresponding resource, but in principle any resource could be
combined with any method.


## LIANA Framework
  
To this end we built a framework to decouple the methods from their corresponding resources.
  
![landingpage](ligrec_pipe.png)
  
  
## Tools

The tools implemented in this repository are:

- CellPhoneDB algorithm (via [Squidpy](https://squidpy.readthedocs.io/en/latest/))
- CellChat
- NATMI*
- Connectome*
- SingleCellSignalR* (SCA)
- iTALK*
  
  
*These methods were re-implemented in LIANA in order to enable them to work with multimeric complexes and to simplify their calls. If you wish to run the original tools please see [LIANA++](https://saezlab.github.io/liana/articles/liana_devel.html).
  
  
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
  
All the resources above are retrieved from [OmniPath](https://omnipathdb.org/),
and more specifically [OmnipathR](https://github.com/saezlab/OmnipathR).
However, individual resources retrieved from the OmniPath web service are not to be affected by this, as each resource expected to be identical to its original form, apart from minor processing imperfections.

`OmniPath` itself serves as a composite CCC resource combining all the ones listed
above + [more](https://doi.org/10.15252/msb.20209923). `OmniPath` also collects
further information about the roles of proteins in intercellular communication
and provides flexible filtering options (see [liana](https://saezlab.github.io/liana/articles/liana_custom_res.html)).  
  
  
  
## Install LIANA  
Installation of all dependencies takes ~5-10 minutes.
  
  
  
## Tutorial
See a [tutorial](https://saezlab.github.io/liana/articles/liana_tutorial.html) how to use LIANA to run all methods and resource from above!
The tutorial with the test data takes minutes to complete!
  
## LIANA++

  
  
## Contact
We appreciate any feedback, so please do not hesitate to open an issue on the [liana github page](https://github.com/saezlab/liana)!  
  
  
### Citing `LIANA`:
Dimitrov, D., Türei, D., Boys, C., Nagai, J.S., Flores, R.O.R., Kim, H., Szalai, B., Costa, I.G., Dugourd, A., Valdeolivas, A. and Saez-Rodriguez, J., 2021.  Comparison of Resources and Methods to infer Cell-Cell Communication from Single-cell RNA Data.
  bioRxiv. [10.1101/2021.05.21.445160v1](https://www.biorxiv.org/content/10.1101/2021.05.21.445160v1)

#### Also, if you use the `OmniPath` CCC Resource for your analysis, please cite:
Türei, D., Valdeolivas, A., Gul, L., Palacio‐Escat, N., Klein, M., Ivanova, O., Ölbei, M., Gábor, A., Theis, F., Módos, D. and Korcsmáros, T., 2021. Integrated intra‐and intercellular signaling knowledge for multicellular omics analysis. Molecular systems biology, 17(3), p.e9923.
https://doi.org/10.15252/msb.20209923
  
#### Similarly, please consider appropritaly citing any of the methods and/or resources implemented in liana, that were particularly relevant for your research!
