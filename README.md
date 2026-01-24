# GaussianEarth
Code for "An EOF-Based Emulator of Means and Covariances of Monthly Climate Fields"


## Data and Setup



## Reproducibility pipeline (generated artifacts and order)

All derived products are written to `${GAUSSIAN_EARTH_WORK}` (default: `data/work/`).

### Emulator Generation

**Stage 1 — EOF bases (per variable)**
- Script: `compute_basis.jl`
- Reads: CMIP NetCDF: `CMIP_ROOT/{scenario}/{var}/*.nc`
- Writes: `WORK/{var}_basis.hdf5`

**Stage 2 — Projections and ensemble stats (per variable × scenario)**
- Script: `project_onto_basis.jl`
- Reads: `WORK/{var}_basis.hdf5` + CMIP NetCDF
- Writes: `WORK/{var}_{scenario}_projection.hdf5` (datasets: `projection`, `global mean`, `ensemble mean`, `ensemble standard deviation`)

**Stage 3 — Regression driver (GMT time series)**
- Script: `regression_variable.jl`
- Reads: `WORK/tas_{scenario}_projection.hdf5` (`global mean`)
- Writes: `WORK/tas_ensemble_yearly_average.hdf5` (datasets: `{scenario}`, `{scenario}_normalized`, `normalization`)

**Stage 4 — Mean regression in EOF amplitude space (per field)**
- Script: `mean_regression.jl`
- Reads: `WORK/{field}_{scenario}_projection.hdf5` + `WORK/tas_ensemble_yearly_average.hdf5`
- Writes: `WORK/{field}_mean_regression.hdf5` (and optionally `WORK/{field}_mean_regression_quadratic.hdf5`)

**Stage 5 — Export covariance targets for optimizer (per field)**
- Script: `export_covariances.jl`
- Reads: `WORK/{field}_{scenario}_projection.hdf5` + `WORK/tas_ensemble_yearly_average.hdf5`
- Writes: `WORK/{field}_covariances.hdf5` (datasets: `temperature`, `scale`, `covariances 1`…`12`)

**Stage 6 — Covariance model fit (JAX/Optax, per field)**
- Script: `jax_code/optimization.py`
- Reads: `WORK/{field}_covariances.hdf5`
- Writes: `WORK/{field}_covariances_model.hdf5` (datasets: `L1`…`L12`)

**Stage 7 — Assemble emulator artifact (per field)**
- Script: `emulator.jl` (tas) / `emulator_hurs.jl` (hurs)
- Reads: `WORK/{field}_basis.hdf5` + `WORK/{field}_mean_regression.hdf5` + `WORK/{field}_covariances_model.hdf5`
- Writes: `WORK/{field}_gaussian_model.hdf5` (datasets: `mean`, `L model`, `basis`, `scale factor`)

**Additional — Pattern scaling baseline**
- Script: `pattern_scaling.jl`
- Reads: CMIP NetCDF + `WORK/tas_ensemble_yearly_average.hdf5`
- Writes: `WORK/{field}_pattern_scaling.hdf5`

To run the full pipeline (skipping existing outputs):
```bash
julia --project=. reproduce.jl
```

### Reproduce paper plots
