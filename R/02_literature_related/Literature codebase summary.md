# Literature codebase summary

1. `01_extract_lit_space.Rmd` 

	**Part 1:** breast cancer lit space
	
	- we collect breast cancer literature space from two traits `ieu-a-1126` (Breast cancer) and `finn-a-C3_BREAST` (Malignant neoplasm of the breast)
	
	- this queries the PRODUCTION graph using epigraphdb package. This lit space is more limited that dev version (v1.2)
	
	
	Output: `breast_cancer_litspace_prod.csv`
	
	
	**Part 2:** risk factor traits lit spaces
	
	Input: `trait_manual_ivw_subtypes_merged.tsv` list of 213 traits from MR analysis
	
	- we exclude traits that we know beforehand that are not good for lit space extraction, or with duplicate, similar traits. In total, we query 154 traits (some will become combined lit spaces)
	- we also identify zero-spaces, and tidy spaces and save triple/pair counts
	
	Output: 
	
	* `traits_marked_for_lit_analysis.tsv` - 154/213 marked for analysis
	* `lit_spaces_finalset.RData` - raw lit space for those 154 traits
	* `lit_spaces_finalset_tidy.RData` - tidy space for 90 traits with non-empty spaces
	* `traits_marked_for_lit_analysis_with_size.tsv` - same as first, but included lit spaces sizes, highlighting those qith zero-size lit spaces. 
	
	
	Then we get lit spaces summary sizes by categories:
	
	* `available_litspace_counts_by_exposure_cat.tsv`-- Table 2 in paper
	
		
	We also combine similar traits' lit spaces to present them as a single trait:
	
	```
	1. Height, sitting height, standing height
	2. Weight, Body mass index (BMI)
	3. Overweight, Obesity class 2, Obesity class 1, Extreme body mass index
	4. Hip circumference, hip circumference, Waist-to-hip ratio
	5. Fresh fruit intake, Dried fruit intake, Cherry intake
	6. Had menopause, Age at menopause (last menstrual period)
	7. Age started hormone-replacement therapy (HRT), Ever used hormone-replacement therapy (HRT)
	8. Number of live births, Age at last live birth
	9. VLDL particles
	10. HDL particles
	```
	
	Output:
	
	* `lit_spaces_combined_traits.RData`
	* `lit_spaces_combined_traits_tidy.RData`
	* `traits_marked_for_lit_analysis_combined.tsv` 
	* `lit_space_stats.tsv` --- lit space counts (individual and combined) 
	* `all_tidy_spaces.xls` xls that combines `lit_spaces_finalset_tidy.RData` and `lit_spaces_combined_traits_tidy.RData` to be used as supl data 9
		
	* Also making a figure	
		
	
2. `02_literature_overlap.Rmd`

	* Loads breast cancer literature space `breast_cancer_litspace_prod.csv` created in 01.. script, does space tidy up (like for all traits in 01), then extract 4 sequential triples.

	* For each molecular trait (i.e. trait with anchor):
		 - load tidy lit space form `lit_spaces_finalset_tidy.RData`
		 - extract and link trait triples
		 - make mini sankey of trait triples (for exploratory purposes)
		 - overlap literature spaces of trait and breast cancer
		 - the output has the full_sankey or subset_sankey (excludes n=1 (or more)) 
		 - save overlap intermediate terms


	* For lifestyle traits, the process is similar.

		 - we load tidy trait space, by specifying whether it is a single-ID trait, or a combined trait: data is loaded from the selected RData
		 - for lifestyle trait there is no anchor, so we first only extract a single triple from the trait (it is allowes to loop back on itself)
		 - we also make a trait-only sankey 
		 - then we connect BC space to trait space via that one triple A-B  and allow trait triple to backwards connect to another triple X-A
		 - lit overlap sankey
		 - save intermediates


	* Overlap intermediates are saved in `results/literature_outputs/sankey_terms_storage/lit_terms_TRAITNAME.csv` These will later be reviewed and manually mapped to OpenGWAS traits. Some reusing of the mapping is done in `merge_term_to_gwas_tables.R`
	
	
	* All functions used by this script are stored within the Shinyapp `app3_sankey_app/functions_literature.R`



	
3. Legacy scripts in `legacy/`


	* `review_literature_mapping.Rmd` - legacy; exploratory -- will reuse parts or GWAS linking analysis
	* `explore_literature_by_years.Rmd` - breast cancer space exploration by years
	* `explore_literature_networks_viz.Rmd` - drawing lit space networkd with D3 (legacy now)



<br><br><br>
	Breast cancer lit space comparison based on prod vs v1.2:
		
	|                | prod  | v1.2  |   
	|----------------|-------|-------|
	| unique triples | 65738 | 74333 | 
	| unique PMIDs    | 23809 | 26392 |   
	| missing PMIDs from the other version | 3032  | 449   |  
	
	 ~1000 out of 3032  from 2020/2019 