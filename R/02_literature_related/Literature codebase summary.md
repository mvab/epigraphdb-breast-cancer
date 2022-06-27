# Literature codebase summary

1. `01_extract_lit_space.Rmd` **add summary**

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
* `lit_space_stats.tsv` --- individual 

	
Also making a figure	
	








	
2. `review_literature_mapping.Rmd` - legacy; exploratory -- will reuse parts or GWAS linking analysis

	Breast cancer lit space comparison based on prod vs v1.2:
	
|                | prod  | v1.2  |   
|----------------|-------|-------|
| unique triples | 65738 | 74333 | 
| unique PMIDs    | 23809 | 26392 |   
| missing PMIDs from the other version | 3032  | 449   |  

 ~1000 out of 3032  from 2020/2019 