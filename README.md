### Description
GxE is intended to be a suite of R functions for implementing gene-environment interaction (GxE) tests on a genome-wide scale. 

Currently, it has one function `JMA2df` that implements a joint 2-df SNP and SNPxE association test using GWAS summary statistics. It is based on the statistical method proposed in Manning et al (2011), and as used in Zhang et al (2021+). Please refer/cite *both* articles if this function is used:

1. Manning, A.K., LaValley, M., Liu, C.T., ... , Dupuis, J. (2011) "Meta-analysis of gene-environment interaction: joint estimation of SNP and SNP x environment regression coefficients". *Genet Epidemiol* 35(1), 11-18, https://doi.org/10.1002/gepi.20546

2. Zhang, W., Venkataraghavan, S., Hetmanski, J.B., ... , Ray, D., Beaty, T.H. (2021+, *under review*) "Detecting gene-environment interaction for maternal exposures using case-parent trios ascertained through a case with nonsyndromic orofacial cleft". 

**Key Words:** Gene-environment test; GWAS summary statistics; Interaction test; Meta-analysis

### Requirements
R (>= 3.0.1)


### How to Install within R
```{r}
require(devtools)
source_url("https://github.com/RayDebashree/GxE/blob/main/joint_metaanalysis_2df_GxE.R?raw=TRUE")
```
It is recommended to download/copy the stand-alone R program in this repository, save it in your local directory of choice and `source()` it from your local directory. When a new version of the software is available, older versions may be removed from this repository, and the above `devtools::source_url()` technique may not work.


### Usage

#### Simple example
```{r}
JMA2df(SNP.coef, GXE.coef, SNP.se, GXE.se, SNP.GXE.cov=NULL, trio.package.input=FALSE)
```
#### Arguments
| Input | Description |
| ---: | --- |
| `SNP.coef` | The vector of estimated SNP effects of size `k` from a 2 degree of freedom (df) gene-environment model, where `k`(>1) is the number of studies from which the SNP effects are obtained. In other words, a model with main effect G and interaction effect GxE are fit for each of the `k` independent studies to be meta-analyzed, and `SNP.coef` contains all the coefficient estimates for G. |
| `GXE.coef` | The vector of estimated GxE effects, from the afore-mentioned model, of size `k` . |
| `SNP.se` | The vector of estimated standard errors of SNP effects, from the afore-mentioned model, of size `k`. |
| `GXE.se` | The vector of estimated standard errors of GxE effects, from the afore-mentioned model, of size `k`. |
| `SNP.GXE.cov` | The vector of estimated covariance between SNP and GxE effects, from the afore-mentioned model, of size `k`. |
| `trio.package.input` | If the SNP and GxE effect estimates are obtained from the [R trio package](https://www.bioconductor.org/packages/release/bioc/html/trio.html) (as used in Zhang et al, 2021+), then no input should be given for `SNP.GXE.cov`. |

#### Value
| Output | Description |
| ---: | --- |
| `SNP.coef.JMA` | Jointly meta-analyzed SNP effect for a 2-df model with main effect G and gene-environment interaction effect GxE. |
| `GXE.coef.JMA` | Jointly meta-analyzed GxE effect for the afore-mentioned 2-df model. |
| `SNP.se.JMA` | Jointly meta-analyzed standard error of SNP effect estimate for the afore-mentioned 2-df model. |
| `GXE.se.JMA` | Jointly meta-analyzed standard error of GxE effect estimate for the afore-mentioned 2-df model. |
| `SNP.GXE.cov.JMA` | Jointly meta-analyzed covariance between SNP and GxE effect estimates for the afore-mentioned 2-df model. |
| `wald2df.stat.JMA` | 2-df Wald test statistic for the joint test of SNP and GxE effects. |
| `wald2df.pval.JMA` | P-value of the 2-df Wald test. P-value below chosen threshold (usually the traditional genome-wide threshold of 5e-7) means either or both the SNP and the GxE effects are statistically significantly different from 0. |

