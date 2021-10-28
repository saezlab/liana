# LIANA 0.0.2

## Changes

* Changed `expr_thresh` to 0.1, based on lack of improvement in performance when using 0.2, hence opted out for the less conservative threshold as default   
* Changed the way that default parameters are passed to each method  
* Enabled housekeeping score aggregation for external methods (needed for revisions) via `.score_housekeep`
* Fixed Bug where external methods could not be called with their default DB. The resource is now always decomplexified
* Seurat Testdata is now properly normalized
* liana_aggragate will now by defaul dissociate complexes - meaning that CellChat Complexes are now appropriately dissociated 
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
