# GaussianEarth

Code for "An EOF-Based Emulator of Means and Covariances of Monthly Climate Fields".

This repo supports two workflows:
1) Regenerate the emulator HDF5 artifacts from slightly pre-processed CMIP data.
2) Load precomputed HDF5s and use the emulator directly.

It also includes a one-shot script to generate all figures in the manuscript.

## Quick Start

- Regenerate emulator artifacts: run `generate_emulator.jl` and `optimization.py`.
- Use precomputed artifacts: point `save_directory` at the HDF5 bundle and run the relevant scripts.
- Generate paper figures: run `generate_figures.jl`.

Both workflows expect `save_directory` to be defined once at the top of the workflow script.

## Paths Configuration

Set `save_directory` in exactly one place per workflow:

- Emulator generation workflow:
  - `generate_emulator.jl` (define `save_directory` once and pass it into the included scripts)
- Figure generation workflow:
  - `generate_figures.jl` (define `save_directory` once before includes)

In the generation workflow, update the `data_directory` in the same workflow entry point. This keeps all included scripts consistent.
The generation workflow will read from the data directory and save HDF5 files to the save directory. The emulator and figure creation workflow will read these files from the save directory. 

## A Note on The Data

## Workflow A: Regenerate the Emulator HDF5s

1) Set paths in `generate_emulator.jl` by editing each included script's `save_directory`
   (and `data_directory` where applicable).
2) Run:

```bash
julia --project=. generate_emulator.jl
```

This will:
- Compute EOF bases.
- Project CMIP fields onto those bases.
- Build the regression variable.
- Fit mean regressions.
- Export covariances for JAX optimization.
- Create pattern-scaling baselines.

### JAX covariance fitting

The low-rank covariance factors used by the emulator are fit in `jax_code`.
After running `generate/export_covariances.jl`, use the JAX workflow to run `optimization.py` and produce:
- `<field>_covariances_model.hdf5`

Then run `emulator.jl` / `emulator_hurs.jl` to assemble:
- `<field>_gaussian_model.hdf5`

## Workflow B: Use Precomputed HDF5s

Point `save_directory` at your HDF5 bundle and use any of the emulator scripts
(`emulator.jl`, `emulator_hurs.jl`, `generate_figures.jl`, etc.).

### HDF5 bundle overview (what each file is)

Per field (`tas`, `hurs`):
- `<field>_basis.hdf5`
  - EOF basis, singular values, lat/lon grids, and metric weights.
- `<field>_<scenario>_projection.hdf5`
  - EOF projections, global mean, ensemble mean, ensemble std for each scenario.
- `<field>_mean_regression.hdf5`
  - Linear regression coefficients of EOF modes vs global mean temperature.
- `<field>_mean_regression_quadratic.hdf5`
  - Quadratic regression coefficients of EOF modes vs global mean temperature.
- `<field>_covariances.hdf5`
  - Empirical covariances used as input to covariance fitting.
- `<field>_covariances_model.hdf5`
  - Fitted low-rank covariance factors (from JAX).
- `<field>_gaussian_model.hdf5`
  - Assembled emulator (means, covariance factors, basis, scale factor).
- `<field>_pattern_scaling.hdf5`
  - Pattern-scaling baseline fits (used in appendix/diagnostics).
- `<field>_full_errors.hdf5`
  - Optional precomputed error arrays (used by `figure_C2.jl`).

Shared:
- `tas_ensemble_yearly_average.hdf5`
  - Global mean temperature time series per scenario.
- `land_sea_mask.hdf5`
  - Land/sea mask used in `figure_7.jl`.
- `new_scenarios.hdf5`
  - Custom scenario used in the case study figures.

Scenarios typically used: `historical`, `ssp119`, `ssp245`, `ssp585`.

## Generate Paper Figures

Set `save_directory` (and optionally `data_directory`) in:
- `generate_figures.jl`

Then run:

```bash
julia --project=. generate_figures.jl
```

This script loads the emulator and executes all figure scripts in `paper_plots/`.

## Using The Emulator