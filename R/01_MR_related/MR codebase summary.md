## MR codebase summary

_very rough draft_

1. `mr_epigraphdb_query.R` - MR-EvE data collection for all MR outcomes:

	The approach: get all MR relationship for all BC outcomes, calualate CIs, keep only those exposure for which the effect CIs don't overlap the null, regardless of p-val. _(We get about 35 relationships where pval >0.05, but most of there are random traits and won't make it to the final analysis anyway)._
	
	Re-extract data for only those exposures from all BC outcomes.

	
	Output: `bc_all_mr_fromCIs.tsv`
	

2. `explore_mr_results.Rmd` - tidy up MR-EvE output and split into categories

	This takes input from the previous step, adds exposure/outcome labels, performs minor filtering, and allows exploring each trait category in interactive plots. The interactive plots are equivalent to RShiny app, but less refined. Also drops a bunch of not needed traits in several categories. The output df can be used for more directed filtering (ignoring actual trait names) (done in the next script).
	
	Output: `tidy_traits_by_cat.tsv`
	

3. `process_mr_results.Rmd` - Traits validation

	* Extract traits with consistent effect (2/3 main datasets)
	
		Output: `trait_for_followup.tsv` (n = 314)
	
	* Perform manual MR in all remaining traits + sensitivity analysis
	
		Output: 
		
		- `redone_MR_fulloutput.tsv` 
		- `redone_MR_fulloutput_sens.tsv`
		- `redone_MR_subsetoutput_ivw.tsv` # only ivw version
	 
	* Perform manual MR on subtype outcome [in other script `mr_results_for_bcac_subtypes.Rmd` on all traits in `trait_for_followup.tsv`] 
	
	* read in subtype MR results produced in the next script (`all_traits_MR_vs_BCAC2020.tsv`)
	
	* join the results for old and new outcomes in a wide table of effect direction
	
	* There are some heatmap / hierarchical clustering plotting  code(currently not reviewed)
	
	Output: `trait_manual_ivw_subtypes_merged.tsv`


4. `mr_results_for_bcac_subtypes.Rmd` - Manual MR on subtypes

	* Perform manual MR on subtype outcomes (new BCAC 2020) on all traits `trait_for_followup.tsv`
	
	* saves data to `mr_evidence_outputs/mr_subtypes/` for each trait individually 
	
	* joins data into full table, which is later read by the previous script. 
	
	Output:
	
	`all_traits_MR_vs_BCAC2020.tsv`
	 


5. `query_confounders.R`