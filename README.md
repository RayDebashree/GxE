### Description
GxE is intended to be a suite of R functions for implementing gene-environment interaction (GxE) tests on a genome-wide scale. 

Currently, it has one function `JMA2df` that implements a joint 2-df SNP and SNPxE association test using GWAS summary statistics. It is based on the statistical method proposed in Manning et al (2011), and as used in Zhang et al (2021+). Please refer/cite *both* articles if this function is used:

1. Manning, A.K., LaValley, M., Liu, C.T., ... , Dupuis, J. (2011) "Meta-analysis of gene-environment interaction: joint estimation of SNP and SNP x environment regression coefficients". *Genet Epidemiol* 35(1), 11-18, https://doi.org/10.1002/gepi.20546

2. Zhang, W., Venkataraghavan, S., Hetmanski, J.B., ... , Ray, D., Beaty, T.H. (2021+, *under review*) "Detecting gene-environment interaction for maternal exposures using case-parent trios ascertained through an orofacial cleft case". 

**Key Words:** Gene-environment test; GWAS summary statistics; Interaction test; Meta-analysis

### Requirements
R (>= 3.0.1)


### How to Install within R
```{r}
require(devtools)
source_url("https://github.com/RayDebashree/PLACO/blob/master/GxE_v0.1.R?raw=TRUE")
```
It is recommended to download/copy the stand-alone R program in this repository, save it in your local directory of choice and `source()` it from your local directory. When a new version of the software is available, older versions may be removed from this repository, and the above `devtools::source_url()` technique may not work.


### Changes
Version 0.1 - February 28, 2021
> First public release of the software.
