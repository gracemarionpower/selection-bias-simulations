# Selection Bias Simulations: Childhood Body Size and Breast Cancer Risk

**Authors:** Grace M. Power, Gibran Hemani
**Date:** 30 June 2025  

This repository contains R code to simulate and visualise how selection bias could produce the observed protective Mendelian Randomization (MR) effect of childhood body size on breast cancer risk (observed OR = 0.59), assuming no true causal effect.

## Overview

The script generates simulated data with genetic instruments and performs IV regression analysis after selection based on body size and cancer status. Two models are considered:

- Additive selection: where selection depends on childhood body size and cancer.
- Interaction-only selection: where selection depends on the interaction between childhood body size and cancer.

The simulations vary selection parameters and examine whether bias alone could account for the observed effect size.

## Output

The code produces two main plots:

1. Estimated bias under additive selection on body size and cancer.
2. Estimated bias under interaction-dependent selection (childhood body size × cancer).

Each plot includes a reference line corresponding to the observed MR effect (log(OR) ≈ -0.527).

