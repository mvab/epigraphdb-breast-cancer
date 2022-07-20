### <span style="color:white"> Integrating Mendelian randomization</span>

#### R/Shiny apps

* [MR-EvE app](https://mvab.shinyapps.io/brca-miner/)

  Breast cancer MR-EvE (Mendelian Randomization "Everything-vs-Everything") app was created to facilitate the exploration of exposures that could be breast cancer risk or protective factors. The app includes 905 exposure traits split into 12 categories, with MR results available for multiple breast cancer outcome GWAS.
  
  <img src="content/figs/app1.png" width="275"/>

* [Heatmaps app](https://mvab.shinyapps.io/MR_heatmaps/)

   The MR effect direction heatmap app includes the results for 309 traits that were selected with filtering steps from the initial MR-EvE exploration dataset. 213 traits (or 105 after FDR correction) have eveidence of effect on at least one validation breast cancer outcome. 

  ![Image](content/figs/app2.png)


* [Literature overlap Sankey plot app](https://mvab.shinyapps.io/literature_overlap_sankey/)

  The Sankey plot app provides a visualisation of literature spaces overlap of risk factor traits (case studies only) and breast cancer. 

  <img src="content/figs/app3.png" width="275"/>


#### Case study reports

  <img src="content/figs/case_studies.png"/>
  
  In the article, we explore four  breast cancer risk/protective factors as case studies in more detail. The investigation of potential mediators identified for these traits from MR-EvE and literature-mined data is presented in invidual report, which were generated with Rmd in an automated way. 

* [Cardiotrophin-1](content/case_study_report_Cardiotrophin-1.html)

* [IGF-1](content/case_study_report_IGF-1.html)

* [Childhood body size](content/case_study_report_Childhood_body_size.html)

* [Age at menopause](content/case_study_report_Age_at_menopause.html)


#### EpiGraphDB queries examples 

[https://github.com/mvab/epigraphdb_mr_literature_queries](https://github.com/mvab/epigraphdb_mr_literature_queries)





