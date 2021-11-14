# LIANA 0.0.4
## New Features
* [`scConnect`](https://academic.oup.com/bioinformatics/article/37/20/3501/6273571#307516453)-inspired
interactions scores were implemented.
* Filtered OmniPath resource was updated


# LIANA 0.0.3
## Improvements
* The R re-implementation of CellPhoneDBv2's permutation algorithm was optimized
to work with sparse matrices (and is now uqicker), and set as the default option
in LIANA (replacing the re-implementation of the same algorithm from squidpy)  

* Custom proportion filtering - Connectome and CytoTalk are now not filtered by
expr_prop as this affects the way that their scores are calculated, since they 
require all clusters/cluster pairs to be present to appropriately scale or
normalize their scores.


## Bug Fixes
* Fixed an issue where logFC was assigned only the value of the ligand   



# LIANA 0.0.2

## New Features
* [`CytoTalk`](https://advances.sciencemag.org/content/7/16/eabf1356)-inspired 
Cross-talk Scores (CTSs) were added.
In contrast to the CytoTalk, in our calculation CTS with ligand or receptor with
PEM of 0 are assigned 0 CTS. Furthermore, we use the inverse of the non-self-talk
scores calculated in CytoTalk to also allow for autocrine signalling interactions,
and thus make cytotalk comparable to the rest of the methods in LIANA.
Finally, as part of LIANA, CytoTalk's re-implemented scores would not take 
account of complexes and we also apply liana-specifc filtering such as according
to `expr_prop`.

## Changes

* Changed `expr_thresh` to 0.1, based on lack of improvement in performance when using 0.2, hence opted out for the less conservative threshold as default   
* Changed the way that default parameters are passed to each method  
* Enabled housekeeping score aggregation for external methods (needed for revisions) via `.score_housekeep`
* Fixed Bug where external methods could not be called with their default DB. The resource is now always decomplexified
* Seurat Testdata is now properly normalized
* liana_aggragate will now by defaul dissociate complexe for CellChat Complexes
* Added tests for changes


# LIANA 0.0.1

## New Features

`liana_wrap` and `liana_aggragate` as the two highest level functions to run all the methods in liana and aggragate them, respectively.

### Re-implemented the following scores in LIANA:   

* logFC  
* NATMI specificity edges  
* Connectome scaled_weights    
* CellPhoneDB algorithm   
* SingleCellSignalR LRScore

each called via `liana_call`, which leverages the statistics provided by `liana_pipe`,

### Others

* Not re-implemented method score names now start with `call_*`

* `decomplexify` and `recomplexify` as functions used to dissociate complexes in resources and 
account for complexes of the re-implemented methods above   

* `liana_aggragate` - a handy wrapper to aggregate results   

* `LIANA` and `LIANA++` are now the user-friendly and benchmark version of LIANA, respectively   

* A webpage with vignettes showing the validity of the re-implemented methods, a developer/benchmark-focused vignette, and a vignette to customize OmniPath  

## Bug fixes

A number of fixes were implemented thanks to the early stage users. Thanks you.
