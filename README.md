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
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE) # iTALK throws a warning...
### Modified version of SingleCellSignalR (SCA)
devtools::install_github(repo = "https://github.com/CostaLab/SingleCellSignalR_v1", subdir = "SingleCellSignalR")    
