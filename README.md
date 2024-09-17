# Simulated sample splitting to overcome Winner’s Curse in summary-level Mendelian Randomization

<span style="color:darkred;">**Note: `mr.simss` is still under active development.**</span>
 

<br>

This R package, `mr.simss`, has been designed to allow users to easily implement the Mendelian Randomisation (MR) method, **MR-SimSS**. This method aims to eliminate bias induced by ***Winner's Curse*** in MR causal effect estimates, using **summary statistics** from genome-wide association studies (GWASs). The approach is based on a **repeated simulated sample splitting** procedure and works in combination with existing MR methods, such as IVW and MR-RAPS. It also takes into account any potential ***sample overlap*** between the exposure and outcome GWASs and thus, can be applied in both one-sample and two-sample cases. 

The ***main function*** in this package for executing **MR-SimSS** is `mr_simss`. `mr_simss` has several parameters which users can adjust based on their desired form of method implementation. Further discussion regarding these parameters can be viewed in the first vignette, ['Performing MR-SimSS'](https://amandaforde.github.io/winnerscurse/articles/winners_curse_methods.html). As well as estimates of the variant-exposure and variant-outcome associations and corresponding standard errors, the method also requires *(a)* knowledge of the number of samples in both outcome and exposure GWASs, the number of overlapping samples between the two GWASs, and the correlation between the exposure and the outcome, or *(b)* an estimate of the correlation between the SNP-outcome and SNP-exposure effect sizes. Due to this, this R package also contains a ***supplementary function***, `est_lambda`, which allows users to obtain an *unbiased* estimate for *lambda*, the correlation between the variant-outcome and variant-exposure effect sizes, using a conditional log-likelihood approach. It is recommended to use this function when the fraction of overlap and the correlation between exposure and outcome are unknown. Merely the GWAS **summary statistics** of a large set of *unpruned* variants or single nucleotide polymorphisms (SNPs) are required in order to obtain an estimate for *lambda*.


### Installation

You can install the current version of `mr.simss` from **GitHub** with:

``` r
install.packages("remotes")
remotes::install_github("amandaforde/mr.simss")
library(mr.simss)
```


### Winner's Curse in MR 

As stated in [Forde *et al.* (2023)](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1010546), it has been observed that in general, the effect size of a SNP tends to be lower in a replication study than in the GWAS that discovered the SNP-trait association. This observation is due to the phenomenon known as *Winner’s Curse*. In the context of a single discovery GWAS, the term *Winner’s Curse* describes how the estimates of association strength for SNPs that have been deemed most significant are very likely to be exaggerated compared with their true underlying values.

*Winner’s Curse* bias can have many practical consequences, especially with respect to techniques which are reliant on SNP-trait association estimates obtained from GWASs. Thus, as MR is a statistical framework which uses genetic variants as instrumental variables to estimate the magnitude of the causal effect of an exposure on an outcome, *Winner’s Curse* bias can also prove problematic in the MR setting. 

In the case of two-sample MR, if the same GWAS is used to identify instrument SNPs and estimate their effects relative to the exposure, *Winner’s Curse* will result in the overestimation of these SNP-exposure associations. This bias will then propagate into the causal estimate, resulting in a deflation of this estimate. On the other hand, if instrument SNPs are discovered in the same GWAS as that used to estimate the SNP-outcome associations, the causal estimate will be inflated, as discussed in [Jiang *et al.* (2022)](https://academic.oup.com/ije/article/52/4/1209/6961569?login=false). In addition, *Winner’s Curse* has been shown to greatly increase the magnitude of weak instrument bias in these MR analyses. This is explored in detail in [Sadreev et al. (2021)](https://www.medrxiv.org/content/10.1101/2021.06.28.21259622v1.full).


[Jiang *et al.* (2022)](https://academic.oup.com/ije/article/52/4/1209/6961569?login=false)
This bias can have a direct impact on Mendelian randomization estimates calculated using association estimates derived from the discovery GWAS. With a single
genetic variant, the Mendelian randomization estimate can
be expressed as the ratio of the genetic association with the
outcome divided by the genetic association with the exposure.12 With multiple genetic variants, the standard combined estimate (the inverse-variance weighted estimate) is a
weighted mean of these ratio estimates calculated for each
variant.13 Hence, winner’s curse in the ‘exposure’ association estimates would be expected to result in a ‘deflation’
in the Mendelian randomization estimate, whereas winner’s curse in the ‘outcome’ association estimates would be
expected to result in an ‘inflation’ in the Mendelian randomization estimate.
Winner’s curse can be alleviated by selecting genetic
variants and estimating genetic associations in non-overlapping data sets.


[Zheng et al. (2017)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5711966/)
In the case of single-sample MR, utilizing the same sample as a discovery analysis for genetic instruments is not a good idea because estimates of the SNP-exposure association will be biased upwards.
In the case of two-sample MR, genetic associations published in discovery GWAS may overestimate the SNP-trait association, particularly if the GWAS is
underpowered to detect the particular loci. In the
case of GWAS of exposures, this will overestimate
the effect of the genetic instrument relative to the
exposure and result in bias of the causal estimate
towards the null. Likewise, in the case of GWAS of
outcomes, winner’s curse will overestimate the
association between the genetic instrument and the
outcome and lead to a bias in causal estimates away
from the null.


image?


### Weak Instrument Bias in MR

[Zheng et al. (2017)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5711966/)
The use of weak instruments
can bias MR estimates towards the confounded observational estimate in one-sample MR settings and towards the null in two-sample MR settings (with nonoverlapping samples).

[Sadreev et al. (2021)](https://www.medrxiv.org/content/10.1101/2021.06.28.21259622v1.full)
The problem of weak instruments for MR studies has been described previously5–8. MR is
predicated on determining the relationship between exposure and outcome through the lens of
only the known genetic factors for the exposure, which typically represent a small fraction of
the total variance of the phenotype, though this is a fraction of the variance which is less
susceptible to confounding and reverse causation. When the core assumptions are met, the
difference in the outcome across values of the genotype (instrument) can be ascribed to
differences in the exposure. However, if the instrument-exposure association is weak, then a
fraction of the observed association is likely to be due to chance associations with
confounding factors. The weaker the instrument-exposure association, the larger the fraction
of the error term from the SNP-exposure estimation that correlates with confounders. The
consequence of this phenomenon is that MR associations obtained using weak instruments
will be biased. In the case that the exposure and outcome effects are estimated in a single (i.e.
fully overlapping) sample, the bias will be in the direction of the (confounded) observational
association. When SNP-exposure and SNP-outcome estimates are obtained from two nonoverlapping samples, the bias is in the direction of the null. Partial sample overlap will lead
to an estimate that lies between these two extremes. Weak instrument bias is conceptually
similar to regression dilution bias: in the single sample setting due to non-differential
measurement error and in the two sample setting due to classical measurement error9,10. Since
the latter bias is always conservative, studies often go to substantial lengths to avoid any
sample overlap between the exposure and outcome.

[Jiang *et al.* (2022)](https://academic.oup.com/ije/article/52/4/1209/6961569?login=false)
For a ‘one-sample’ Mendelian randomization analysis, in
which the same data set is used for estimating the genetic
associations with the exposure and the outcome, these
chance correlations affect associations with both the exposure and the outcome in a related way, and bias (known as
weak instrument bias) is ‘towards the observational association’ between the exposure and outcome. For a ‘two-sample’ Mendelian randomization analysis, in which genetic
associations with the exposure and outcome are obtained
from independent samples, these chance correlations will
differ between the data sets and so affect associations with
the exposure and the outcome independently, and weak instrument bias is ‘towards the null’.

image?


### MR-SimSS Overview 

*small paragraph here providing brief description of method.. building on what was already stated in the introduction*
