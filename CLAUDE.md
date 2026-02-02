# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

Fourier Mixture Regression (FMR) estimates the distribution of common-variant effect sizes from GWAS summary statistics and Fourier LD scores. The method is described in O'Connor 2021, *Nature Genetics*.

There are two implementations:
- **MATLAB** (`MATLAB/`): the original, fully functional implementation used in the paper
- **Python** (`FMR/`): a work-in-progress port, not yet ready for use

## Architecture

### MATLAB Implementation

The core algorithm lives in `MATLAB/FMR/`:

- `run_FMR.m` — Main entry point. Takes chi-squared stats, Fourier LD scores (`lF`), LD scores (`l2`), LD 4th moments (`l4`), mixture variance grid (`ss`), and sampling times (`tt`). Returns estimated mixture weights (`ww`), scaled variance parameters (`ss`), and LD4M output.
- `LD4M.m` — Stratified LD 4th moments regression. Estimates heritability, kurtosis, and effective number of independent SNPs (`Me`). Used internally by `run_FMR` to determine rescaling.
- `nonnegative_regression_jk.m` / `linear_regression_jk.m` — Block-jackknife weighted regression (nonneg uses MATLAB's `lsqlin` from Optimization Toolbox).

Postprocessing functions in `MATLAB/otherfunctions/`:
- `estimate_h2_quantiles.m` — Estimates quantiles of the heritability distribution from FMR output
- `compute_polygenicity.m` — Computes polygenicity measures (entropy-based, effective, softmax) from mixture weights
- `predict_future_GWAS.m` — Predicts results of larger future GWAS from FMR output

Scripts in `MATLAB/scripts/` reproduce paper figures and run FMR on specific datasets.

### Python Implementation

Located in `FMR/` (a separate Python project with its own `pyproject.toml`). Uses `uv` for dependency management.

- `FMR/src/fmr/cli.py` — CLI entry point (`fmr` command): reads data, runs FMR, outputs `.tall.csv` and `.jk.csv`
- `FMR/src/fmr/io.py` — Data I/O. Reads Fourier scores from Parquet (with scalar metadata) and sumstats from TSV. Merges on SNP column. Defines `DataFMR` dataclass.
- `FMR/src/fmr/regression.py` — Core regression logic, ported from MATLAB

### Data Flow

1. Input: GWAS summary statistics (chi-squared or Z-scores) + precomputed Fourier LD scores
2. LD4M regression estimates heritability and kurtosis to determine rescaling
3. Construct moment equations: Fourier moments (`yyF`), second moments (`yy2`), fourth moments (`yy4`)
4. Stack into a single regression: `X = [lF; l2; l4]`, `yy = [yyF; yy2; yy4]`
5. Nonnegative jackknife regression produces mixture weights
6. Weights are normalized to sum to 1

### Key Conventions

- Mixture components are parameterized by variance `ss` (typically `2^(-7:5)`) and sampling times `tt = sqrt(ss)`
- Jackknife standard errors use 100 blocks by default; SE = SD * sqrt(n_blocks + 1)
- The `lF` matrix columns correspond to `ss`-`tt` pairs with a diagonal structure indexed by `index_matrix`
- `RescaleParam` controls how the mixture grid aligns with the summary statistics scale

## Commands

### Python

```bash
cd FMR && uv sync              # install dependencies
cd FMR && uv run pytest        # run tests
cd FMR && uv run fmr <scores.parquet> <sumstats.tsv> <output_prefix>  # run FMR
```

### MATLAB

Scripts assume `matfiles/` exists in the repo root with precomputed data (Fourier LD scores, LD matrices, sumstats). Download from the Dropbox link in README.md. Run scripts from `MATLAB/scripts/` after adding `../FMR` and `../otherfunctions` to the path.

## Data

- **Fourier LD scores** (Python): Parquet files with columns `SNP`, `l2`, `l4`, `lF_1`...`lF_K`, plus `scalars` in file metadata
- **Summary statistics** (Python): TSV with `SNP` and `Z` columns; SNPs must start with `rs`
- **MATLAB data**: `.mat` files in `matfiles/` containing `lF`, `l2`, `l4`, `SNPs`, `sumstats` structs

## Notebook

This project uses a separate notebook repository for analysis logs. See `notebook/INDEX.md` for a summary of past work.
