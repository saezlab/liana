# LIANA 0.0.6

## New implementations

- added `liana_dotplot` as a basic, but flexible, dotplot function for LIANA output.

## Changes  

- `global_mean` is now calculated in a more efficient manner.

## Bug Fixes  

- By mistake, `assay.type` in `liana_pipe` was passed to `get_logFC` would
result in using the logcounts, rather than the raw counts for logFC calculation.
Has now been appropriately changed.


# LIANA 0.0.5

## Changes
* I now filter the Crosstalk scores to include only those > 0. Otherwise,
LIANA would return all possible combinations of clusters and interactions, which
would be simply NAs and 0s for Crosstalk scores. Should do the same for Connectome (>0).

* `CellChat` and Crosstalk scores/`cytotalk` will no longer by called by default by liana_wrap.
However, it both are available as an option to be passed to the `method` parameter.

* I now filter all methods by `expr_prop`. This is done in a slightly different manner for Connectome's 
scaled weights and crosstalk scores, since they require all pairs/clusters to be present to appropriately
calculate their scores. Thus, for them we filter after we calculate the scores, while for the others methods
we pre-filter.

* We now provide a tutorial how to make use of [intracellular OmniPath](https://saezlab.github.io/liana/articles/liana_intracell.html) as well as
how to [combine LIANA with NicheNet](https://saezlab.github.io/liana/articles/liana_nichenet.html)


# LIANA 0.0.4

* The OmniPath resource had a major update.

* `CellCall` and `Cellinker` resources were added, while talklr was removed. The OmniPath resources itself was filtered further and 1,000
lower quality interactions were excluded. Further improvements were made to all resources,
most of which were minor. Changes worth mentioning were made to ICELLNET (updated to latest resource version), 
CellPhoneDB (was filtered for ambigous interactions), and CellChatDB was filtered for mislaballed interactions.


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
* `Crosstalk scores` inspired by [Cytotalk]((https://advances.sciencemag.org/content/7/16/eabf1356)) were added.
In contrast to the CytoTalk, in our calculation CTS with ligand or receptor with
PEM of 0 are assigned 0 CTS. Furthermore, we use the inverse of the non-self-talk
scores calculated in CytoTalk to also allow for autocrine signalling interactions,
and thus make Crosstalk scores comparable to the rest of the methods in LIANA.
Finally, as part of LIANA, CytoTalk's re-implemented scores would not take 
account of complexes and we also apply liana-specifc filtering such as according
to `expr_prop`. Worth noting, we only re-implement the cross-talk scores, but we
don't include the intracellular part of Cytotalk.

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
