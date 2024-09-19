# Simulated sample splitting to overcome *Winner’s Curse* in summary-level Mendelian Randomization

<br>

[$\star$]{style="color: darkred;"} <span style="color:darkred;">**Note:** </span> **`mr.simss` is still under <u> active
development</u>.**

<br>


This R package, `mr.simss`, has been designed to allow users to easily implement
the Mendelian Randomization (MR) framework known as **MR Simulated Sample Splitting  (MR-SimSS)**. The MR-SimSS paradigm aims to <u>eliminate bias induced by *Winner's Curse*</u> in MR causal effect estimates, using only **summary statistics** from genome-wide association studies (GWASs). The approach is based on a <span style="color:darkblue;">**repeated simulated sample splitting**</span> procedure and works in combination with existing summary-level MR methods, such as IVW and MR-RAPS. In addition, the framework is designed to <u>assist in circumventing biases (e.g. weak instrument bias) arising due to sample overlap</u> between exposure and outcome GWASs, and thus, can be applied in both one-sample and two-sample cases. The <u>main function</u> in this package for executing MR-SimSS is
`mr_simss`. `mr_simss` has several parameters which users can adjust based on
their desired form of method implementation. Further discussion regarding these
parameters can be viewed in the first vignette, ['Performing
MR-SimSS'](https://amandaforde.github.io/mr.simss/articles/perform-MR-SimSS.html), while a detailed derivation of the MR-SimSS framework can be found in the second vignette, ['Performing
MR-SimSS'](https://amandaforde.github.io/mr.simss/articles/derive-MR-SimSS.html).


### Installation

You can install the current version of `mr.simss` from **GitHub** with:
``` r
install.packages("remotes")
remotes::install_github("amandaforde/mr.simss")
library(mr.simss)
```

<br>




### MR-SimSS Overview
*small paragraph here providing brief description of method.. building on what was already stated in the introduction*


The method in question, namely MR Simulated Sample Splitting (MR-SimSS), aims to
extend the benefits of the ‘threesample’ design approach to the single sample
setting, ensuring minimal weak instrument bias and avoiding Winner’s Curse in
the process. In our work, we have explicitly focused on establishing an MR
method that can be used when only GWAS summary-level data, in the form of
variant-exposure and variant-outcome association estimates and standard errors,
is available. It is assumed that this data is derived from two large exposure
and outcome GWASs, which have been conducted with two sets of non-overlapping,
partially overlapping or fully overlapping samples. Being a summary-level
approach, the proposed method benefits from a key advantage, in the sense that
it is very computationally efficient. MR-SimSS works by imitating the splitting
of an individual-level data set into three portions, with the first fraction
reserved for instrument selection and the other two used to independently
estimate the variant-exposure and variant-outcome associations. Taking advantage
of asymptotic conditional distributions, MR-SimSS repetitively simulates
association estimates in each fraction conditional on the known estimates in the
full data set. On each iteration, instrument variants are selected according to
the simulated estimates in the first fraction, while a two-sample MR method,
such as IVW or MR-RAPS, is then applied to the estimates from the other two
fractions to obtain a causal effect estimate. MR-SimSS then averages over the MR
estimates generated on each iteration to provide a final exposure-outcome causal
effect estimate. This act of averaging over a great number of iterations ensures
reduced variance in the final estimate.
