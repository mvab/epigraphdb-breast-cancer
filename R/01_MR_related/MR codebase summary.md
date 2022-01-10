## MR codebase summary

_very rough draft_

1. `mr_epigraphdb_query.R` - MR-Eve data collection for all MR outcomes:

	This uses approach: get everthing, get CIs, keep only those that don't overlap the null, regardless of p-val. get about 35 that have pval >0.05, but most of thsre are random trait and won't make it to the final analysis anyway
	
Output: `bc_all_mr_fromCIs.tsv`
	

2. explore_mr_results.Rmd

takes input from previous step, adds all exposure/outcome things, some minor filetering,  and allows exploring each trait category interactive plots. Equivalent to RShiny app, but less refined. Also drops a bunch or not needed traits in several categories. the output can be used for more directed filtering (ignoring acruak trait names)

output: `tidy_traits_by_cat.tsv`

3. `process_mr_results.Rmd`

* redo manual MR in all reaming traits + sensitivity