# Prevalence estimation with imperfect tests

Code and aggregate results for the manuscript *Prevalence Estimation in Infectious Diseases with Imperfect Tests: A Comparison of Frequentist and Bayesian Logistic Regression Methods with Misclassification Correction*.

## Data

The individual-level data (HIV, syphilis and hepatitis B screening, n = 11,452) belong to the Secretaría de Salud Pública y Seguridad Social de Pereira (Colombia). They are not anonymised and cannot be shared here. Access can be requested from the corresponding author, subject to agreement with the data owner.

The scripts expect a file `vih.rds` in the repository root with one row per individual and these variables:

| Variable | Description |
|---|---|
| `edad_gr` | Age in years |
| `sex` | Sex (1 = male) |
| `vih_1` | HIV result, serial algorithm Determine → Bioline (0/1) |
| `reactivo_sif1` | Syphilis rapid test result (0/1) |
| `reactivo_hepb` | Hepatitis B rapid test result (0/1) |
| `tip_pob_HSH`, `tip_pob_LGTBIQ`, `tip_pob_OTRAS POBLACIONES`, `tip_pob_TRABAJADOR SEXUAL` | Key population indicators (reference: general population) |

No model objects are included, because objects fitted with `glm()` or `stan_glm()` store a copy of the data.

## Structure

```
Nuevos_resultados/
  Script/                        R scripts
    reanalisis_revision.R        Revised analysis (use this one)
    vih1.R, sifilis.R, ...       Scripts of the original submission
  Modelo_logistico_Bayesiano/    Stan models (bic_logit.stan for the revision)
  Resultados/revision/           Aggregate results of the revised analysis
  Resultados/*.rds               Aggregate results of the original submission
  Modelos/traceplot_*.png        Trace plots of the original submission
function_predict.R               Helper used by the original scripts
```

## Reproducing the revised analysis

With `vih.rds` in the repository root, run from the root:

```bash
Rscript Nuevos_resultados/Script/reanalisis_revision.R vih
Rscript Nuevos_resultados/Script/reanalisis_revision.R sif
```

Requires R ≥ 4.4 and the packages `tidyverse`, `rstan`, `rstanarm`, `posterior`, `numDeriv`, `MASS` and `binom`. Each run takes about ten minutes; outputs are written to `Nuevos_resultados/Resultados/revision/`.

The revised script:

- applies the Rogan–Gladen correction as (p − (1 − Sp)) / (Se + Sp − 1);
- fits the Liu et al. model by maximum likelihood with and without constraining the misclassification rates to [0, 0.5);
- computes intervals from each model's own uncertainty;
- reports the marginal standardized prevalence (and, as a complement, the prevalence at the mean covariate profile);
- refits the Bayesian model with internal correction under four prior specifications;
- reports convergence diagnostics (R̂, bulk and tail ESS, divergences, tree depth, E-BFMI) for every parameter.
