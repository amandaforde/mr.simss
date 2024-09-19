# Simulated sample splitting to overcome *Winner’s Curse* in summary-level Mendelian Randomization

<br>

[$\star$]{style="color: darkred;"} <span style="color:darkred;">**Note:** </span> **`mr.simss` is still under <u> active
development</u>.**

<br>


This R package, `mr.simss`, has been designed to allow users to easily implement
the Mendelian Randomization (MR) framework known as **MR Simulated Sample Splitting  (MR-SimSS)**. The MR-SimSS paradigm aims to <u>eliminate bias induced by *Winner's Curse*</u> in MR causal effect estimates, using only **summary statistics** from genome-wide association studies (GWASs). The approach is based on a <span style="color:darkblue;">**repeated simulated sample splitting**</span> procedure and works in combination with existing summary-level MR methods, such as IVW and MR-RAPS. In addition, the framework is designed to <u>assist in circumventing biases (e.g. weak instrument bias) arising due to sample overlap</u> between exposure and outcome GWASs, and thus, can be applied in both one-sample and two-sample cases. The <u>main function</u> in this package for executing MR-SimSS is
`mr_simss`. `mr_simss` has several parameters which users can adjust based on
their desired form of method implementation. Further discussion regarding these
parameters can be viewed in the vignette titled ['Performing
MR-SimSS'](https://amandaforde.github.io/mr.simss/articles/perform-MR-SimSS.html). In addition, ['MR-SimSS: The algorithm'](https://amandaforde.github.io/mr.simss/articles/MR-SimSS-algorithm.html) describes the MR-SimSS algorithm as well as the intuition for the approach, while a detailed derivation of the necessary formulae used in the construction of MR-SimSS can be found in ['Deriving
MR-SimSS'](https://amandaforde.github.io/mr.simss/articles/derive-MR-SimSS.html).

<br> 

#### Installation

You can install the current version of `mr.simss` from **GitHub** with:
``` r
install.packages("remotes")
remotes::install_github("amandaforde/mr.simss")
library(mr.simss)
```

<br>

[$\star$]{style="color: darkblue;"} **Note:** A detailed description of how *Winner's Curse* and weak instrument bias can impact MR causal effect estimates can be found [here](https://amandaforde.github.io/mr.simss/articles/wc_wib.html).

<br>

#### Overview of the MR-SimSS framework


As mentioned in ['Winner's Curse and weak instrument bias in MR'](https://amandaforde.github.io/mr.simss/articles/wc_wib.html), the most popularised solution to the problem of *Winner's Curse* in MR studies is to follow the <span style="color:darkblue;">**'three-sample' design**</span> approach, in which three independent samples are used to select instrument variants, estimate variant-exposure associations and obtain variant-outcome association estimates. In addition to ensuring that neither sets of association estimates suffer from *Winner's Curse* bias, this approach guarantees that any weak instrument bias in the MR causal effect estimate will be in the direction of the null. Given the benefits of the <span style="color:darkblue;">**‘three-sample’ design**</span> approach, **MR Simulated Sample Splitting (MR-SimSS) aims to extend these benefits to the one- and two-sample settings**, <u>ensuring minimal weak instrument bias and avoiding *Winner’s Curse*</u> in the process. 


A key advantage of MR-SimSS is that it has been designed for use in a setting in which **only GWAS summary-level data**, in the form of variant-exposure and variant-outcome association estimates and standard errors, is available. The MR-SimSS framework assumes that this data is derived from two large exposure and outcome GWASs, which have been conducted with <u>two sets of non-overlapping, partially overlapping or fully overlapping samples</u>. Being a summary-level approach, MR-SimSS is considered to be very **computationally efficient**. 

MR-SimSS works by **imitating the splitting of an individual-level data set into three portions**, with the first fraction reserved for instrument selection and the other two used to independently estimate the variant-exposure and variant-outcome associations. Taking advantage of <span style="color:darkred;"><u>**asymptotic conditional distributions**</u></span>, MR-SimSS repetitively simulates association estimates in each fraction conditional on the known estimates in the full data set. On each iteration, instrument variants are selected according to the simulated estimates in the first fraction, while a two-sample MR method, such as IVW or MR-RAPS, is then applied to the estimates from the other two fractions to obtain a causal effect estimate. MR-SimSS then <u>averages over the MR estimates generated on each iteration</u> to provide a final exposure-outcome causal effect estimate. This act of averaging over a great number of iterations ensures reduced variance in the final estimate.
