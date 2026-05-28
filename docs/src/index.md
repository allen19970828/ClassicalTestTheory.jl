# CTT.jl Documentation

Welcome to the documentation for **CTT.jl**, a modern and high-performance Julia package for Classical Test Theory and Psychometrics.

## Overview

`CTT.jl` provides an extensive suite of tools for:
- **Item Analysis**: Difficulty, point-biserial, biserial, and item-deleted reliability.
- **Reliability Estimation**: Cronbach's Alpha, Guttman's Lambdas, and McDonald's Omega (PAF based).
- **Validity & Scaling**: Criterion validity, disattenuation, Z-scores, and T-scores.
- **Advanced Correlations**: Polychoric, tetrachoric, polyserial, and biserial correlations via 2-step MLE.
- **Exploratory Factor Analysis**: Principal Axis Factoring, parallel analysis, and KMO tests.

## Why CTT.jl?

Built in pure Julia, it is significantly faster than standard R equivalents (`CTT`, `psych`), while providing the exact same (or better) level of psychometric rigor. Crucially, it comes with built-in defense mechanisms like **pairwise deletion for missing values** and **correlation matrix smoothing** to ensure your analyses don't crash on messy, real-world data.

See the [API Reference](@ref) for detailed documentation of all functions.
