# TRIX Shiny App – Bodrum & Milas (2013-2022)

Interactive R/Shiny app for TRIX analysis in Bodrum & Milas:
- Yearly summaries and 2013→2022 delta
- Station-level trends (Mann-Kendall + Sen)
- Dataset-level trends
- ITS (break at 2020) and DiD

## Citation 
https://doi.org/10.5281/zenodo.17102904

## How to run
```r
install.packages(c("shiny","shinythemes","shinycssloaders","DT","readxl","dplyr",
                   "tidyr","stringr","purrr","ggplot2","scales","broom","janitor",
                   "zoo","zyp","Kendall","trend","fixest","zip","readr"))
shiny::runApp()
