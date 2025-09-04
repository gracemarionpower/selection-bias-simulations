# Selection Bias Simulations: Childhood Body Size and Breast Cancer Risk

**Authors:** Grace M. Power, Gibran Hemani
**Start date:** 30 June 2025  
**Last updated:** 4 September 2025

This repository contains R code to simulate and visualise how selection bias could produce the observed protective Mendelian Randomization (MR) effect of childhood body size on breast cancer risk (observed OR = 0.59; log(OR) ≈ −0.527), assuming no true causal effect.

## Overview

We simulate correlated genetic risk scores (GRS) for childhood and adult body size and generate latent traits where adult tracks childhood (ϕ_track). Selection into the analysed sample depends on:

- mean-centered childhood body-size category (0/1/2),
- breast cancer status,
- their interaction (child × cancer).

> By default, adult body size is *not* directly selected on (adult selection coefficient = 0), but we still estimate and report an adult effect. Selection on child (and the interaction) plus tracking can induce associations with adult in the selected sample.

Estimation uses a **2SRI (control-function) logistic model**:
1) First-stage linear models of ordered body-size categories on both GRS.  
2) Second-stage logit of cancer on predicted child/adult and their residuals.

---

## Simulation design

- **Selection grids**
  - Child selection: `{0, −0.25, −0.5}`
  - Cancer selection: `{0, −0.25, −0.5}`
  - Interaction selection: `{0, −0.25, −0.5}`
- **Replicates:** 500 (configurable)
- **Per-replicate N:** 246,511
- **Key parameters:** `rg = 0.67`, `cancer_prev = 1/7`, `ϕ_track = 0.35`, `h2_child = 0.10`, `h2_adult_target = 0.10`
- **Seeds:** sims `set.seed(1407)`, R² block `set.seed(20250901)`

---

## Outputs

1. **Main plot:** mean log(OR) across replicates for both effects  
   - **Child effect** (solid) and **Adult effect** (short-dashed)  
   - Error bars: mean ± 1.96×SD across replicates  
   - Faceted by interaction selection; colour encodes cancer selection  
   - Red dashed reference at the observed MR log(OR) (label shown once)

2. **Selection R² decomposition table** (Shapley-averaged McFadden R²)  
   - Shares attributed to: child term, cancer term, child×cancer interaction  
   - Adult selection fixed at 0, so R² shares do not depend on adult.

3. **Final table** merging **R² shares** with **child & adult estimates**  
   - Per parameter cell: `Mean logOR`, `SD logOR`, `Mean SE` for **Child** and **Adult**, plus three R² percentages.

---

## How to run

```r
# Run the script top-to-bottom.
# Avoid dplyr/MASS verb conflicts (select, filter) by either:
library(MASS); library(dplyr)  # load dplyr after MASS
# or prefix dplyr verbs: dplyr::select(), dplyr::distinct(), etc.

# To change the number of replicates:
sim_results <- purrr::map_dfr(1:500, function(i) { ... })
