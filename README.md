![Kegg pathway analysis and AbHAC network analysis](http://www.cng.fr/cagekid/img/cagekid_logo.jpg)



## Note


The data and results folder are not made public before publication.

To access those folders, please contact : mehran dot karimzadeh at uhnresearch dot ca

## Dependencies

The analysis relies on the following R packages:

[AbHAC](https://github.com/mehrankr/AbHAC)

[KeggPA](https://github.com/mehrankr/KeggPA)

[ggplot2](http://cran.r-project.org/web/packages/ggplot2/index.html)

[pheatmap](http://cran.r-project.org/web/packages/pheatmap/index.html)


## Installing dependencies

```R
install.packages("devtool")
require(devtools)
install_github("AbHAC", username="mehrankr")
install_guthub("KeggPA", username="mehrankr")
install.packages("ggplot2")
install.packages("pheatmap")
```

## Running the codes

In linux/mac environment, change your directory to the src folder.

You would need data and results folder (not provided publicly) in order to run the scripts.


```R
Rscript pathway_analysis.R 0.05 1 
#The first number of the FDR cutoff threshold
#The second number is the minimum number of patients in gene must be mutated in so it is included in the analysis
Rscript abhac_analysis.R 0.05 1
#Parameters are similar to before
#You can also try:
Rscript abhac_analysis.R -h
Rscript pathway_analysis.R -h
```
