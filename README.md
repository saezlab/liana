# Cell-Cell Communication Investigation


## Objectives

The continuous developments of single-cell RNA-Seq (scRNA-Seq) have sparked
an immense interest in understanding multi-cellular crosstalk. Multiple
methods and resources that aid the investigation of cell-cell communication
have been recently published.
However, these methods and resources are usually in a fixed combination of a
method and its corresponding resource, but in principle any resource could be
combined with any statistical method. Yet, it is largely unclear the
difference that the choice of resource and method can have on the predicted
CCC events. Thus, we attempt to shed some light on this topic via a
systematic overview of how different combinations might influence CCC
inference, by decoupling the methods from their corresponding resources.

As such, we compare all combinations between 15 resources and 7 methods and
explore the effect on downstream results.

We provide the resources and methods as a pipeline for further use in this
repository.


## Methods

The methods implemented in this repository are:

- CellPhoneDB algorithm (via Squidpy)
- CellChat
- NATMI
- Connectome
- SingleCellSignalR (SCAomni)
- iTALK
- scTalk - to be implemented
- cell2cell - to be implemented


## Resources

### Ligand-receptor resources

The following intercellular signalling (ligand-receptor interaction)
resources are accessible by this package:

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
processing imperfections. OmniPath as a composite resource we use in four
varieties: the full OmniPath intercellular network, only ligand-receptor
interactions, quality filtered (50 percentile consensus score cut off), and
ligand-receptor only quality filtered.

### Random and default

Moreover, a Randomized resource can be generated via reshuffling any of the
abovementioned using the `BiRewire` package, and each tool can be run with
its 'Default' resource, the dataset used in its original publication.


## Dependencies

### Modified version of `SingleCellSignalR` (SCA)

```{r}
# To avoid an iTALK warning:
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE)

devtools::install_github(
    repo = "https://github.com/CostaLab/SingleCellSignalR_v1",
    subdir = "SingleCellSignalR"
)
```
