## MR codebase summary

## File summary

MR workflow scripts (details below) - NB all V3:

```
├── 01_mr_epigraphdb_query_V3.R
├── 02_explore_mr_resultsV3.Rmd
├── 03sub_validate_mrV3.Rmd
├── 04_process_mr_resultsV3.Rmd
├── 05_query_mreve_mediatorV3s.R
├── 06sub_mreve_mediators_validationV3.R
├── mr_related_functions.R
```

Apps:

```


├── app2_heatmaps_app
│   ├── app.R
│   ├── functions_copy_from_mreveapp.R (copy! - must store in the app)
│   ├── heatmap_functions.R
│   ├── heatmap_static.R (can be run to generate the static figure and the RDS input to the app)
│   ├── supl_forestplots.R - adhoc create forest plot of FDR passed results
```

Legacy / supplementary:

```
├── make_supl_data1_V3.R
├── make_abbs_supl.R
├── metadat_collect.R

├── adhoc_mr_testing.R
├── add_protein_gene_info.R
├── protein_names_testing.R
└── review_opengwas_riskfactors.Rmd

├── app1_MR-EvE_app (legacy)
├── app1_MR-EvE_app_2files (not maintained)

```

## Workflow scripts


1. **MR-EvE data collection for all breast cancer outcomes** `01_mr_epigraphdb_queryV3.R` 

	Extract data for BCAC 2017 3 BC outcomes + do FDR

	
	Output: 
		
	* `bc_all_mr_fromCIs.tsv` 
	* `all_mreve_bc_results.csv` - full res, to be used as Supl data 1, with extra columns added.
	

2. **Tidy up MR-EvE output and split it into categories** `02_explore_mr_resultsV3.Rmd` 

	Using the output from the previous step, we add exposure/outcome labels, perform minor filtering/exclusions, and explore each trait category in interactive plots. The interactive plots are equivalent to the RShiny app, but less refined. The output df can be used for more directed filtering (ignoring actual trait names) (done in the next script).
	
	Exclusions:
	* male-only sample traits
	* 'raw' UKB traits
	* anthropometric traits that are limb measurements, older versions of the same data with smaller sample sizes, Neale lab UK Biobank GWAS if MRC-IEU version was available
	* anything that does not fall into 12 categories defined in functions
	
	Output: `tidy_traits_by_cat.tsv`
	


3. **Traits processing and validation summary**
 `04_process_mr_results.Rmd` 
	* Extract traits with consistent effect (2/3 main datasets)
	
		Output: `trait_for_followup.tsv` 
		
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

	Perform MR as validation on  BCAC 2017 and 2020 outcomes (+ sensitivity analyses)
	
	Input: `trait_for_followup.tsv`  from `04_process_mr_results.Rmd` 
	
	* Perform MR in BCAC 2017
			
		- `redone_MR_fulloutput.tsv` 
		- `redone_MR_fulloutput_sens.tsv`
		- `redone_MR_subsetoutput_ivw.tsv` (this is used in mediaotr validation)
	 
	* Perform MR in BCAC 2020 (saved separately and joined in a single table as separate step)
		
		- `all_traits_MR_vs_BCAC2020.tsv`	
		- `all_traits_sensMR_vs_BCAC2020.tsv`	
		
		**UPD** this step is done prior to the one above

5. **Query and process potential mediators from MR-EvE** `05_query_mreve_mediators.R`

	For the final set of traits in `trait_manual_ivw_subtypes_merged.tsv` (across BCAC 2017 outcomes only), we run MR-EvE queries to extract confounders, mediators, colliders, reverse intermediates. 
	
	Initially, we extract all relationships with a high p-value threshold (all results will be manually validated later, so it does not matter). Also not restricting search by pval of med->out, as will be using validation from the previous script to filter those.
	
	When we get a list of all potential meds per trait, we validate their effect in 
	`06sub_mreve_mediators_validation.R` script. 
	
	Inputs:
	`redone_MR_subsetoutput_ivw.tsv` from 03sub - trait-BC validated results
	`redone_MRmeds_subsetoutput_ivw.tsv` from 06sub - trait-mediator validated results
	
	For all identified exp-med-out relationships in MR-EvE, we pull numbers from validation tables, merge, and export as validated. 
	
	
	Outputs:
	`med_extracted_all_r3.csv` - all potential mediators, not validated
	`mediators_counts_per_traits.csv` - counts; validated
	`med-table-validated.xlsx` - validated mediator results 



6. **Run validation for identified mediators** `06sub_mreve_mediators_validation.R`

	We manually re-run MR for all identified mediators for each risk factor trait.
	
	Input: 
	`med_extracted_all_r3.csv`
	Output:
	`redone_MRmeds_fulloutput.tsv`
	`redone_MRmeds_fulloutput_sens.tsv`
	`redone_MRmeds_subsetoutput_ivw.tsv`