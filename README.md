# ClassicalTestTheory.jl

A high-performance, robust, and modern Julia package for **Classical Test Theory (CTT)** and **Psychometrics**. 

Inspired by the widely used R packages `CTT` and `psych`, `CTT.jl` provides an extensive suite of tools for item analysis, reliability and validity estimation, advanced correlations (e.g., polychoric), and exploratory factor analysis (EFA). Designed with real-world data in mind, it safely handles missing values through pairwise deletion and correlation matrix smoothing.

## Installation

You can install `ClassicalTestTheory.jl` from the Julia REPL. Press `]` to enter the Pkg REPL mode and run:

```julia
pkg> add ClassicalTestTheory
```

*(Note: If the package is not yet registered, you may need to add it via its GitHub URL.)*

## Key Features

- **Item Analysis**: Difficulty, point-biserial, biserial correlation, and Cronbach's alpha if deleted.
- **Reliability**: 
  - Cronbach's Alpha
  - Guttman's Lambda (1 through 6)
  - McDonald's Omega (powered by Principal Axis Factoring)
- **Validity & Scaling**
- **Criterion Validity (`criterion_validity`)**
- **Scale Correlations (`scale_cor`)**
- **Disattenuated Correlation (`disattenuated_cor`)**
- **Standard Scores**: `zscore`, `tscore`
- **Advanced Correlations**: Fast and pure-Julia implementation of Biserial, Polyserial, Tetrachoric, and Polychoric correlations using 2-step Maximum Likelihood Estimation (MLE).
- **Exploratory Factor Analysis (EFA)**: 
  - Principal Axis Factoring (PAF) to avoid Heywood cases.
  - Seamless integration with `FactorRotations.jl` for orthogonal (e.g., Varimax) and oblique (e.g., Promax) rotations.
  - Parallel Analysis via Monte Carlo simulations to suggest the optimal number of factors.
  - Kaiser-Meyer-Olkin (KMO) test for sampling adequacy.
- **Robust Missing Data Handling**: Built-in pairwise deletion and automatic correlation matrix smoothing (`cor_smooth`) to ensure non-positive definite matrices don't crash your analysis.

## Quick Start

### 1. Basic Item Analysis and Reliability

```julia
using ClassicalTestTheory

# Prepare a response matrix (Rows = Subjects, Columns = Items)
responses = [
    1 1 0 1;
    1 0 1 1;
    0 0 1 0
]

# 1. Item Analysis
item_res = item_analysis(responses)
println("Difficulty: ", item_res.difficulty)
println("Point-biserial: ", item_res.pbis)

# 2. Reliability (Cronbach's Alpha & McDonald's Omega)
alpha = cronbach_alpha(responses)
omega = mcdonald_omega(responses)
println("Cronbach's Alpha: ", alpha)
println("McDonald's Omega: ", omega)

# 3. Exploratory Factor Analysis (EFA) - Principal Axis Factoring (PAF)
efa_res = efa(responses, 1, rotation="varimax")
println("Factor Loadings:\n", efa_res.loadings)
```

### 2. Advanced Correlations (Polychoric)

If you are dealing with Likert scale data, Pearson correlation assumes continuous variables and often underestimates the true relationship. You can easily compute the Polychoric correlation matrix:

```julia
# Compute the pairwise polychoric correlation matrix
R_poly = polychoric_matrix(X)
```

### 3. Exploratory Factor Analysis (EFA)

Stop guessing the number of factors. Use Parallel Analysis and run a robust EFA:

```julia
using FactorRotations

# 1. Determine the number of factors using Parallel Analysis
n_factors = parallel_analysis(X; cor_type=:polychoric)

# 2. Check Sampling Adequacy
overall_kmo, item_kmo = kmo(pairwise_cor(X))
println("Overall KMO: ", overall_kmo)

# 3. Run EFA with Varimax rotation (automatically smooths bad matrices)
efa_res = efa(X, n_factors; cor_type=:polychoric, rotation=Varimax())

# Extract factor loadings and communality
println(efa_res.loadings)
println(efa_res.communality)
```

## Supported Types

Functions in `CTT.jl` accept `AbstractMatrix` inputs, meaning they seamlessly support matrices of type `Matrix{Float64}`, `Matrix{Int}`, or arrays containing missing values like `Matrix{Union{Missing, Float64}}`.

## License

This project is licensed under the MIT License.
