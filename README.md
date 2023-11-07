# MR-SimSS: Mendelian Randomisation method that combats Winner's Curse using a simulated sample splitting approach

**Note: `mr.simss` is still under active development.** 

<br>

This R package, `mr.simss`, has been designed to allow users to easily implement the Mendelian Randomisation (MR) method, **MR-SimSS**. This method aims to eliminate bias induced by ***Winner's Curse*** in MR causal effect estimates, using **summary statistics** from genome-wide association studies (GWASs). The approach is based on a **repeated simulated sample splitting** procedure and works in combination with existing MR methods, such as IVW and MR-RAPS. It also takes into account any potential sample overlap between the exposure and outcome GWASs and thus, can be applied in both one-sample and two-sample cases. 

The ***main function*** in this package for executing **MR-SimSS** is `mr_simss`. `mr_simss` has several parameters which users can adjust based on their desired form of method implementation. Further discussion regarding these parameters can be viewed in the first vignette, ['Performing MR-SimSS'](https://amandaforde.github.io/winnerscurse/articles/winners_curse_methods.html). As well as estimates of the SNP-exposure and SNP-outcome associations and corresponding standard errors, the method also requires *(a)* knowledge of the number of samples in both outcome and exposure GWASs, the number of overlapping samples between the two GWASs, and the correlation between the exposure and the outcome, or *(b)* an estimate of the correlation between the SNP-outcome and SNP-exposure effect sizes. Due to this, this R package also contains a ***supplementary function***, `est_lambda`, which allows users to obtain an *unbiased* estimate for *lambda*, the correlation between the SNP-outcome and SNP-exposure effect sizes, using a conditional log-likelihood approach. It is recommended to use this function when the fraction of overlap and the correlation between exposure and outcome are unknown. Merely the GWAS **summary statistics** of a large set of *unpruned* SNPs are required in order to obtain an estimate for *lambda*.


### Installation

You can install the current version of `mr.simss` from **GitHub** with:

``` r
install.packages("remotes")
remotes::install_github("amandaforde/mr.simss")
library(mr.simss)
```


### Winner's Curse in MR 


*small paragraph here detailing the issue of winner's curse in MR*


### Weak Instrument Bias in MR

*small paragraph here detailing the issue of weak instrument bias in MR*


### MR-SimSS Overview 

*small paragraph here providing brief description of method.. building on what was already stated in the introduction*
