# COMP_PLS-SEM

A PLS-SEM (Partial Least Squares Structural Equation Modeling) framework for **compositional data**, using the ILR (Isometric Log-Ratio) transformation to handle compositional variables within the PLS-SEM pipeline.

## Getting Started

**Use `tutorial_code.r` as the main entry point to understand the full analysis pipeline.** It walks through every step from data preprocessing to model evaluation in a clear, sequential manner.

## Data Description

The dataset contains four constructs based on the Theory of Planned Behavior (TPB):

- **Attitude**, **PBC** (Perceived Behavioral Control), **BI** (Behavioral Intention) — compositional latent variables. Each is measured by three-part compositions (every 3 columns form one composition). Raw data must undergo closure and ILR transformation before analysis.
- **Subject Norms (SNnew)** — not a latent variable. The two columns are already ILR-transformed (see Kogovšek et al., 2013).

## File Descriptions

| File | Description |
|------|-------------|
| **`tutorial_code.r`** | **Main tutorial script (recommended entry point).** Demonstrates the complete workflow: data preprocessing → ILR transformation → PLS-SEM estimation → PLSpredict cross-validation → compositional block-level metrics → reliability evaluation (CR, Cronbach's α, R², Bootstrap). Uses `data_sub.xlsx`. |
| `Compositional Code.R` | Full analysis script, similar to the tutorial but uses `TPB.xlsx`. Additionally includes Q² evaluation, compositional loadings via ILR inverse transformation, and 10-fold cross-validation comparing PLS / LM1 / LM2 predictions at both LV and MV levels. |
| `comp_function.R` | Compositional data utility library: closure, Helmert matrix, CLR/ILR transforms and their inverses, compositional regression (real↔compositional, compositional↔compositional), Aitchison inner product, etc. |
| `seminr.R` | Custom PLS-SEM modeling functions: measurement model definition (`constructs`, `composite`, `reflective`), structural model definition (`relationships`, `paths`), inner weighting schemes (path weighting / factorial), path coefficient estimation, HTMT, total effects, etc. |
| `PLSpredict.R` | PLS-PM estimation and prediction: `evaluate_simplePLS` (simple PLS-PM algorithm), `PLSpredict` (predict new data from a trained model), `validatePredict` (K-fold cross-validation computing RMSE, MAPE, MAD for both PLS and LM benchmarks). |
| `evaluation.r` | Model evaluation functions: Cronbach's α, Composite Reliability (CR, with support for compositional grouping via COMP-CR), AVE, HTMT, R², Q² (Blindfolding), and Bootstrap confidence intervals. |
| `generation.R` | Simulation data generation for Monte Carlo studies. Constructs covariance matrices under formative-formative and reflective-reflective settings. |
| `real_outer_loadings.csv` | Outer loadings matrix exported from the PLS model. |
| `real_outer_loadings-manual.csv` | Manually reorganized loadings used for compositional loadings ILR inverse transformation in `Compositional Code.R`. |

## Dependencies

**R packages:** `readxl`, `dplyr`, `psych`, `caret`, `Metrics`, `lmtest`, `MASS`, `mvtnorm`, `seminr`

**Data files (not included, must be prepared separately):**
- `data_sub.xlsx` — used by `tutorial_code.r`
- `TPB.xlsx` — used by `Compositional Code.R`

## References

- Kogovšek, T., et al. (2013). ILR transformation for compositional data in social science.
