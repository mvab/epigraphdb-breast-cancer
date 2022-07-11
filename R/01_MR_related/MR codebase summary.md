## MR codebase summary

## File summary

MR workflow scripts (details below):

```
├── 01_mr_epigraphdb_query.R
├── 02_explore_mr_results.Rmd
├── 02_explore_mr_results.nb.html
├── 03_process_mr_results.Rmd
├── 03_process_mr_results.nb.html
├── 03sub_validate_mr.Rmd
├── 04_query_mreve_mediators.R
├── 04sub_mreve_mediators_validation.R
├── 05_case_study_report.Rmd
├── mr_related_functions.R
```

Apps:

```
├── app1_MR-EvE_app
│   ├── app.R
│   ├── functions.R
├── app2_heatmaps_app
│   ├── app.R
│   ├── functions_copy_from_mreveapp.R (copy! - must store in the app)
│   ├── heatmap_functions.R
│   ├── heatmap_static.R (can be run to generate the static figure and the RDS input to the app)
```

Legacy / supplementary:

```
├── adhoc_mr_testing.R
├── protein_names_testing.R
└── review_opengwas_riskfactors.Rmd

├── app1_MR-EvE_app_2files (not maintained)

```

## Workflow scripts


1. **MR-EvE data collection for all breast cancer outcomes** `01_mr_epigraphdb_query.R` 

	The approach: get all MR relationships for all BC outcomes, calculate CIs, keep only those exposures for which the effect CIs don't overlap the null, regardless of the p-value. _(We get about 35 relationships where pval >0.05, but most of there are random traits and won't make it to the final analysis anyway)._
	
	Re-extract data for only those exposures from all BC outcomes.

	
	Output: `bc_all_mr_fromCIs.tsv` (N= 2332 -> 1970)
	

2. **Tidy up MR-EvE output and split it into categories** `02_explore_mr_results.Rmd` 

	Using the output from the previous step, we add exposure/outcome labels, perform minor filtering/exclusions, and explore each trait category in interactive plots. The interactive plots are equivalent to the RShiny app, but less refined. The output df can be used for more directed filtering (ignoring actual trait names) (done in the next script).
	
	Exclusions:
	* male-only sample traits
	* 'raw' UKB traits
	* anthropometric traits that are limb measurements, older versions of the same data with smaller sample sizes, Neale lab UK Biobank GWAS if MRC-IEU version was available
	* anything that does not fall into 12 categories defined in functions
	
	Output: `tidy_traits_by_cat.tsv`(n = 1643 -> 905)
	

3. **Traits processing and validation summary**
 `03_process_mr_results.Rmd` 
	* Extract traits with consistent effect (2/3 main datasets)
	
		Output: `trait_for_followup.tsv` (n = 309)
		
		**~~ MR validation is done in a separate script:** `03sub_validate_mr.R`
	
	* read in MR validation results produced in the `03sub_validate_mr.R` script for BCAC 2017 and 2020 and join them in a wide format.
		- `redone_MR_fulloutput.tsv` 
		- `all_traits_MR_vs_BCAC2020.tsv`
	
	The script produces multiple outputs:
	
	
	* `table1_counts_by_exposure_cat.tsv` - count by trait category and outcomes, at various stages of processing - basis for the dynamically made Table 1 in the paper (made with _flextable_ in the Rmd)
	
	* `trait_manual_ivw_subtypes_merged.tsv` - table of effect direction in wide format for all outcomes

	* Venn diagrams 309 / 213 / 171 / 168
	* `passed_multiple_testing.tsv` marking traits that passed MTC
	* `all_data_with_sens_filters.xlsx` - XLS file (a sheet per outcome) IVW/MR main results including sens tests + MTC
	* `all_data_validation.tsv` - all results (other MR methods as a single table)

	
	 
4. **MR results validation** `03sub_validate_mr.R` 

	Perform MR as validation on 309 traits on all BCAC 2017 and 2020 outcomes (+ sensitivity analyses)
	
	Input: `trait_for_followup.tsv`  from `03_process_mr_results.Rmd` 
	
	* Perform MR in BCAC 2017
			
		- `redone_MR_fulloutput.tsv` 
		- `redone_MR_fulloutput_sens.tsv`
		- `redone_MR_subsetoutput_ivw.tsv` (this is used in mediaotr validation)
	 
	* Perform MR in BCAC 2020 (saved separately and joined in a single table as separate step)
		
		- `all_traits_MR_vs_BCAC2020.tsv`	
		- `all_traits_sensMR_vs_BCAC2020.tsv`	

5. **Query and process potential mediators from MR-EvE** `04_query_mreve_mediators.R`

	For the final set of traits in `trait_manual_ivw_subtypes_merged.tsv` (across BCAC 2017 outcomes only), we run MR-EvE queries to extract confounders, mediators, colliders, reverse intermediates. 
	
	Initially, we extract all relationships with a high p-value threshold (all results will be manually validated later, so it does not matter). Also not restricting search by pval of med->out, as will be using validation from the previous script to filter those.
	
	When we get a list of all potential meds per trait, we validate their effect in 
	`04sub_mreve_mediators_validation.R` script. 
	
	Inputs:
	`redone_MR_subsetoutput_ivw.tsv` from 03sub - trait-BC validated results
	`redone_MRmeds_subsetoutput_ivw.tsv` from 04sub - trait-mediator validated results
	
	For all identified exp-med-out relationships in MR-EvE, we pull numbers from validation tables, merge, and export as validated. 
	
	
	Outputs:
	`med_extracted_all_r3.csv` - all potential mediators, not validated
	`mediators_counts_per_traits.csv` - counts; validated
	`med-table-validated.xlsx` - validated mediator results 



6. **Run validation for identified mediators** `04sub_mreve_mediators_validation.R`

	We manually re-run MR for all identified mediators for each risk factor trait.
	
	Input: 
	`med_extracted_all_r3.csv`
	Output:
	`redone_MRmeds_fulloutput.tsv`
	`redone_MRmeds_fulloutput_sens.tsv`
	`redone_MRmeds_subsetoutput_ivw.tsv`



7. **Case study report** `05_case_study_report.Rmd` 


	The report is split into 3(+1) parts:
	
	1. MR results for case study trait for all breast cancer outcomes (BCAC 2017 and 2020)
	2. Overview of potential mediators identified from MR-EvE data; their validation with two-step and multivariable MR (MVMR)
	3. Overview of potential mediators identified from literature-mined data; their validation with two-step, bidirectional and MVMR
	4. _Optional_ validation of results with female-only exposure data (if available)

	The script depends on these input files:
	
	- `01_MR_related/results/mr_evidence_outputs/all_data_with_sens_filters.xlsx`
	- `02_literature_related/results/literature_outputs/lit_space_stats.tsv`
	- `02_literature_related/results/literature_outputs/sankey_terms_storage/....` - specific case study file in that location
	- local GWAS data in a different project if the optional step is run