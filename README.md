# EarthCovar

Code for "An EOF-Based Emulator of Means and Covariances of Monthly Climate Fields".

This repo supports two workflows:
1) Regenerate the emulator HDF5 artifacts from slightly pre-processed CMIP data. 
2) Load precomputed HDF5s and use the emulator directly.

It also includes a one-shot script to generate all figures in the manuscript.

## Quick Start

- Regenerate emulator artifacts from CMIP data: run `generate_emulator.jl` and `optimization.py`.
- Use precomputed artifacts to e.g. run the emulator: point `save_directory` at the HDF5 bundle and run the relevant scripts.
- Regenerate (from CMIP data) the ground truth metrics needed to reproduce paper figures: run `precompute_ground_truth_bundle.jl`.
- Generate paper figures: run `generate_figures.jl`.

## Paths Configuration

Set `save_directory` at the top of each of the following files. This will house all HDF5 save-out files.

- Emulator generation workflow:
  - `generate_emulator.jl` 
  - `precompute_ground_truth_bundle.jl` 
- Figure generation workflow:
  - `generate_figures.jl` 

In the generation workflow, update the `data_directory` in the same two scripts. 
The generation workflow will read from the data directory and save HDF5 files to the save directory. The emulator and figure creation workflow will read these files from the save directory. 

## A Note on The Data

In order to work with CMIP data directly, the scripts expect the data directory to be structured as follows: 
`${data_directory}/${scenario_name}/${variable_name}/r**i1p1f1_${scenario_name}_${variable_name}.nc`
where each file in the leaf directories is the full history of a given variable in a given scenario for a given ensemble member.

The MPI-ESM1.2-LR  data used in this study may be downloaded directly from the CMIP6 archive. Files may need to be concatenated to span the scenario timeperiods.


Some figures use ground-truth statistics from the underlying CMIP ensemble (means, standard
deviations, and histogram samples). To avoid redistributing the full CMIP archive, we have generated a compact bundle of the necessary ground truth metrics. You can
generate the same bundle using:

```bash
julia --project=. precompute_ground_truth_bundle.jl
```

Set `save_directory` and `data_directory` at the top of that script. It writes
`ground_truth_bundle.hdf5` into `save_directory`. If the file exists, the figure scripts
will automatically load from the bundle instead of reading raw CMIP netCDFs.

The bundle includes:
- `grid/metric` and `grid/mask` (if present in `save_directory`)
- `tas` ensemble mean/std monthly + annual mean for `historical`, `ssp119`, `ssp245`, `ssp585`
- Figure-specific samples:
  - `samples/figure_2/...` (ssp119 low/high windows)
  - `samples/figure_6/...` (latitude samples) + `stats/figure_6/...` (zonal mean/std)
  - `samples/figure_7/...` and `samples/figure_8/...` for `historical`, `ssp245`, and `ssp585` (for Figures D2 and D3)
  - `stats/figure_A1/selected_indices` + `samples/figure_A1/point_*`

It uses HDF5 compression to keep size reasonable.

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
- `<field>_model.hdf5`

## Workflow B: Use Precomputed HDF5s

Point `save_directory` at your HDF5 bundle and use any of the emulator scripts
(`emulator.jl`, `emulator_hurs.jl`, `generate_figures.jl`, etc.).

### HDF5 artifacts overview (what each file is)

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
- `<field>_model.hdf5`
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

Set `save_directory` (and optionally `data_directory`, in the absence of the ground truth bundle file) in:
- `generate_figures.jl`

Then run:

```bash
julia --project=. generate_figures.jl
```

This script loads the emulator and executes all figure scripts in `paper_plots/`.

## Using The Emulator


The emulator produces mean/variance (and samples) of monthly climate fields conditioned
on a global mean temperature. Typical usage is:

```julia
include("utils.jl")
include("emulator.jl")  # loads `emulator` from save_directory

month = 1
T = regression_variable("ssp245")[10]  # or any scalar temperature you want
emulator.month[1] = month
emulator.global_mean_temperature[1] = T

μ = mean(emulator)          # length 192*96 (flattened grid)
σ2 = variance(emulator)     # per-grid variance, same length
x = rand(emulator)          # one sample realization

μ_grid = reshape(μ, 192, 96)
σ_grid = reshape(σ2, 192, 96)
```

Notes and common operations:
- The emulator is stored in `<field>_model.hdf5`. You can construct it directly with
  `CovarEmulator(μmodel, Lmodel, basis)` or `CovarEmulator(save_directory * field)`. Note that both the linear and quadratic mean fits are stored, but the emulator defaults to the quadratic fit when constructured.
- Outputs are flattened in the basis grid order (192×96 for tas/hurs); `reshape` to map form.
- `mode_mean(emulator)` / `mode_variance(emulator)` give mean/variance in EOF space.
- `emulator_variance(emulator; modes, month, global_mean_temperature)` returns the EOF-space
  covariance for a specific month/temperature.
- `emulator_mean_variance_linear_functionals(observables, emulator)` lets you compute
  mean/variance of linear summaries (e.g., area averages) without sampling.

For scenario-driven runs, use `regression_variable(scenario)` to get the temperature
time series and loop over time/months, as done in `paper_plots/`. In general, see figure generation scripts for detailed examples of emulator usage.


Reach out to G. Geogdzhayev at geogdzh (at) gmail (dot) com with bugs/questions/etc.
