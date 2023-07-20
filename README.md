countprop: Calculate Model-Based Metrics of Proportionality on Count-Based Compositional Data

This package allows estimation of several types of proportionality metrics for count-based compositional data such as 16S, metagenomic, and single-cell sequencing data. The package includes functions that allow standard empirical estimates of these proportionality metrics, as well as estimates based on the multinomial logit-normal model.

To install this package from github, run the following code:
```r
if (!require(devtools)) {
  install.packages("devtools")
  library(devtools)
}
install_github("kevinmcgregor/countprop", dependencies=TRUE, build_vignettes=TRUE)
```

Once the package is loaded into R, you can view the vignette:
```r
library(countprop)
vignette("countprop")
```