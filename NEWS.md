# LIANA 0.1.13 (03.11.23)

- Changed the way that `max` rank is imputed when NAs are presented, or when `return_all` is true. 
Essentially, `RobustRankAggregate` will use the max rank in the matrix, rather than the size of the dataframe.

- Fixed a bug related to newer versions of CellChat with unused argument #75

- Bumped `cell2cell` to its latest `0.6.8` version.

- LIANA's `source` and `target` will not inherit the levels of `idents_col` #99


# LIANA 0.1.12 (24.02.23)

- Added `rank_aggregate` to generate both specificity and magnitude rank
aggregates. Essentially runs `liana_aggregate` twice with different `aggregate_how`
parameters and joins.
- Added `invert_specificity`, `invert_magnitude`, `invert_function` parameters
to `liana_dotplot`.
- Added `decompose_tensor` as a function to run only the decomposition on a 
pre-built Tensor.
- Aggregation can now be performed also via `liana_bysample`, takes `aggregate_how` 
parameter, which allows `magnitude`, `specificity`, or `both`.
- Added `preprocess_scores` function that handles the conversion of liana's
scores to Tensor format.
- Added additional tests related to Tensor-cell2cell


# LIANA 0.1.11 (06.02.23)

- Fixed issue with sample_col in `liana_tensor` functions
- Removed redundant scconnect code
- Merged #89 to fix typo in `liana_aggregate` documentation
- Fix bug with column duplicates in `cellchat_formatDB`.
- passing "all" to `select_resource` will now only return human resources.
- add minor condition improvements to `liana_wrap` PR#92

# LIANA 0.1.10 (23.01.23)

- Fixed issues in the `generate_lr_geneset` function, and other typos and
by mistake hardcoded values.
- Bumped up the version of Tensor-cell2cell to the latest one, and it will now
return error values that could be used to estimate the elbow curve.
- Added an example of the elbow plot used to guide the number of ranks to be
considered.
- Extended options by which the cell2cell module can be loaded.
- Fixed an issue where `return_all` was being passed to external methods
- Added an example elbow curve in the tensor tutorial.
- LIANA's doplot will now keep the order of the original dataframe passed to it.


# LIANA 0.1.9 (13.12.22)

## Changes
- `return_all` parameter was included to `liana_wrap`. `return_all` enables
the return of all interactions by liana, not only the ones that pass the
`expr_prop` threshold. Those that don't pass the threshold are assigned the
worst possible score, and a `lrs.to.keep` flag that indicates whether the
interaction passed the threshold.

- `supp_columns` was included to `liana_wrap`, which allows additional
columns to be added to the output of any method.


## Minor changes
- Min prop for complexes is now explicitly assigned to the minimum expr. prop. across all subunits
This deals with edge cases with non-expression ligand/receptor scores (e.g. z-score) where
the lower score subunit is the one with higher expression proportion. Hence,
this is intended to make all methods consistent according to which interactions
are returned, regardless of which subunit has the lower score.
- `aggregate_how` parameter added to `liana_aggregate` to allow the aggregation
by specificity and magnitude scores.


# LIANA 0.1.8 (08.11.22)
## New Implementations 
- Untargeted between-condition (context/sample) decomposition of cell-cell 
communication latent patterns /w `tensor_cell2cell`. Makes use of `basilisk` to
automatically set-up a conda env for liana.
- added `min_cells` parameter to `liana_wrap`, to exclude any cell identity
which does not pass a minimum cells threshold.

## Changes
- Mouse Consensus resource is now provided by default.
- The intracellular OmniPath vignette was removed. An updated and more user-friendly one 
will be provided in next updates. In the meantime, the old one can still be downloaded
from [drive](https://drive.google.com/file/d/1lqxHhmz0Jq7eIuQAe0SxvInGgo2U-RlC/view?usp=share_link)
- Source and Target titles are now plotted by the `liana_dotplot`
- added explicit error if `idents_col` was not found in metadata/colData

# LIANA 0.1.7 (13.10.22)
## Changes

- Changed the way ties are handles in liana_aggregate. Namely, I previously
assigned the minimum rank, but this resulted in ties getting lower p-values
than they should, particularly for scores with a lot of ties (e.g. CPDB p-value).

- Fixed an issue where some subunits of 0 `expr_prop` would not get filtered.
This was observed due to previous changes to `.filt_liana_pipe` in 0.1.6, where
some subunits were filtered before `recomplexifying`.

- Fixed an issue in which some NATMI complexes would be missing due to `recomplexify`
being done on both .expr and .sum columns. These are now seperated into `columns`, which
are the ones for which I account for complexes, and `add_cols`, the ones that
are additional - no need to account for complexes (e.g. also `global_mean`).

## Minor Changes
- I now refer to [SCPubr](https://enblacar.github.io/SCpubr-book/index.html) in liana's tutorial.
- Throw exception for NAs in cell idents
` Remove duplicated rows from orthologous resource

## Minor changes

- I now refer to SCpubr in the tutorial for more advanced plots.

# LIANA 0.1.6 (11.08.22)

## Changes
- Fixed an issue where interactions with complexes will not filtered be according to
`expr_prop` for some methods. I now filter twice - once via `.filt_liana_pipe`
for computational speed, and once after `recomplexify` to also remove the 
complexes with `expr_prop` <= X. Will now also filter `crosstalk_scores` to `expr_prop`.

- In the edgecase of complexes with subunits with equal expression, LIANA's internal
methods will not arbitrarily discard duplicate complex interactions.

- Will now return `expr_prop` for each method. Note that this information is
discarded by `liana_aggregate`.

- `liana_doplot` function is now more explicit in the way interactions are
selected. Will now take `topn` and return the highest ranked interactions.
Size of dots is also more distinguishable by default and the user can now 
pass a customizable value for the size range.

- Added a `rank_method` helper function to rank single methods according to
`specificity` and magnitude.

- Removed ~20 bad quality interactions from the `Consensus` resource.

- Minor changes on filtering SCE object in `liana_pipe` to ensure all
complex subunits are present in the sce


# LIANA 0.1.5 (04.07.22)

## Changes
- Re-implemented the `RRA` method from  Kolde et al., 2012, as a consequence of
the removal of the `RobustRankAggregate` package from CRAN.

- Integrate `generate_homologs` with OmniPath's `homologene` database.
This allows homology conversion by simply passing an organism ID. Also, handles
complicated cases, such as complex subunits with one-to-many mapping homologs.




# LIANA 0.1.4 (18.06.22)

## Changes
- Add `prod_weight` to NATMI's score. This is the weight that both Connectome and
NATMI suggest for between-condition comparisons. Add NATMI to the housekeeping
aggregate ranking.

- Enable weighing of interactions by cell pairs (using a DF in which each 
cell pair has  an assigned weight). This would typically be done by spatial
constraints, etc. These weights can also be used to mask any cell-pair interactions
which are not relevant (by assigning weights of 0). This currently assumes that
the weights would be between 0 to 1 - to be extended. Tutorial on this /w appropriate
spatial weight generation to be written. 


## Minor Changes
- By default, the base for logFC will now be automatically assigned depending 
on the object passed to LIANA, i.e. `.antilog1m` for SCE will use 2 as base,
and Euler's number for Seurat. One could also pass the base they wish to use
via `liana_wrap`.
- Automate website deployment to gh-pages and run R checks on push.


# LIANA 0.1.3 (15.05.22)

## Changes
- Changed the aggregation columns of `liana_aggregate`, as in some cases
methods would assign different subunits as the minimum, which results in
redundancies for the same complex. As such, `liana_aggregate` will now
return only the complex columns, nevertheless, the methods will still return 
both the minimum (lowest expressed subunit) and it's corresponding complex.

- `base` used to calculate logFC (via `get_log2FC`) can now be passed as 
a parameter to `liana_wrap` via `liana_pipe.params`.
Passing `NaN` to base would result in log2FC calculation using the raw counts
without any pre-processing (e.g. no batch correction, etc).

The base is by default set to 2, assuming that log2 transformation is performed
following library size normalization, and thus preserving the normalization, 
while reverting back to ~counts.

## Minor Changes
- Extended the heatmap to allow filtering down to certain cell types.
- Removed redundant/leftover code

## New implementations
- A chord plot for interaction frequncies was included.



# LIANA 0.1.2 (03.05.22)

## New Implementations
- Frequency Heatmap available via the `heat_freq` functions, added due to being common requests.
This heatmap was inspired by CellPhoneDB and CellChat.

## Changes
- Extended basic tutorial to accommodate new heatmap plots.

## Minor changes
- Allow labels to be passed to `liana_dotplot`
- Cleaned up docs, dependencies, examples, and warnings


# LIANA 0.1.1 (26.04.22)

## Changes
- Change the order of  non-expressed genes and empty droplet filtering. 
I now appropriately filter cells in the `sce` object *after* limiting the gene
universe to ligands and receptors in the resource.

## Minor changes
- Appropriately pass `verbose` to `.filter_sce`

- Silence expected warning in `cellchat_formatDB`


# LIANA 0.1.0 (20.04.22)

## Changes

- For `logFC_mean`, rather than normalizing the counts by library size, I instead inverse log the counts and use those to calculate log2FC. This is to preserve any prior correction of the counts,
i.e. mainly for consistency with the rest of the methods.

## Minor Changes

- Flipped `x` and `y` axes dotplot according to feedback.

- Added `.default_fun` parameter to `generate_orthologues`

- Minor code clean up

- Added `examples` to main exported functions docs

- Removed several low quality interactions from the `Consensus` resource

- CellChat will now work with the simplified format of intercell (i.e. Consensus resource).

- CellChat and Squidpy should now be called using `call_cellchat` and `call_squidpy`, respectively.


# LIANA 0.0.9 (23.03.2022)

## New Implementations

- We now provide a tutorial dedicated to `orthology conversion` of the resources in LIANA

## Minor Changes
LIANA will now check if:
- There is enough of a gene intesect between the resource and the data (i.e. if they are both for human)

- There are negative counts

- An assay with normalized counts was been provided (i.e. if the data/logcounts slot is scaled).

- LIANA will now convert non-sparse matrices to sparse.


# LIANA 0.0.8 (05.03.2022)

## Changes
- Reduced dependencies (specifically `Seurat` and `OmniPathR`)

## Minor Changes
- testthat tests external methods only if requested explicitly
- Readme updated - clarified and accordingly describes the `Consensus` resource as default


## New Features
- `idents_col` is can now be explicitly passed to `liana_wrap`, if not provided defaults to the active
idents/colLabels for SCE and Seurat, respectively.
- `verbose` param allows to omit any messages and warnings from LIANA
- `assay` can now be passed explicitly when working with a Seurat object, defaults to the active one otherwise


# LIANA 0.0.7

## Changes

- LIANA will now use the `Consensus` resource by default. This is a highly-literature supported resource, generated using similar
filtering steps as the 'OmniPath' (old default) resource. This resource is similar in size (~4,700 interactions), but contains a 
higher complex and curation content.

- All resources might show some very minor changes related to an update of UniProt IDs and homology-conversion improvements.

- LIANA now uses `mean0` to account for heteromeric complexes, i.e. the mean is computed, unless there is a value of 0, then 0 is returned.
This means that any complex, the subunit of which is not expressed is filtered. LIANA now also appropriately accepts any custom function to
account for complexes.

- `liana_aggregate` now groups by ligand.complex and receptor.complex as well as the subunits,
and hence returns a all of those columns


## Minor Changes

- Added option to show complexes on dotplot and is now the default option

- Documentation improvements

- `decomplexify` function is now exported

- `liana_aggregate` will no longer return a median_rank, it's largely redundant.

- re-arranged the column order of `liana_aggregate` due to the addition of .complex columns

- Replaced min0 (used to obtain closest to 0 value) to min -> relevant for z-scores used in Connectome.

## Bugs

- Complexes with missing subunits are not correctly assigned as 'missing' and hence filtered/treated as non-expressed.

- Fixed a bug where LIANA will return the minimum subunit expression, instead of the mean for some methods. 
  This stemmed from not properly passing the incorrect `complex_policy` to certain methods, i.e. they were getting a hard-coded value instead.
  
- Remove `decomplexify` logical from `liana_call` and `liana_pipe` -> redundant.

- edge case fix: liana_aggregate should now rank interactions with the same subunits, but coming from different complexes seperately


# LIANA 0.0.6

## New implementations

- LIANA has now been optimized in terms of RAM, by swapping all internal function to rely solely on
the BioConductor single-cell framework (for all internal methods).

- LIANA now accepts both `SingleCellExperiment` and `Seurat` objects as input.

- added `liana_dotplot` as a basic, but flexible, dotplot function for LIANA output. (+ tests)

## Changes  

- LIANA will now perform a basic filtering step, where all genes and cells with 0 summed counts are removed.

- `global_mean` is now calculated in a more efficient manner.

- `assay.type` in `liana_pipe` was passed to `get_logFC` would
result in using the logcounts, rather than the library-normalized counts.

## Bug Fixes  

- Fixed a bug where incorrectly passing method names in different cases results in an error.


## Deprecated

- External LIANA methods (i.e. `call_`) are now deprecated. The pipelines will be maintained solely for power users,
who intend to benchmark the original implementations, but will not be the focus of any downstream analyses.
These will be solely developed for the internal (or re-implemented methods). These still rely on a `SeuratObject` as 
interface, but will now accept both sce and seurat as input.


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
