using NCDatasets, LinearAlgebra, Statistics, HDF5, ProgressBars

include("utils.jl")

save_directory = "PLEASE/SET/YOUR/SAVE/PATH/HERE/"
data_directory = "PLEASE/SET/YOUR/DATA/PATH/HERE/"

# Determine EOFs in the Historical Period for 1 Realization
include("generate/compute_basis.jl")
# Project onto the EOFs
include("generate/project_onto_basis.jl")
# Compute Regression Variable (global, ensemble, and yearly averaged surface temperature)
include("generate/regression_variable.jl")
# Compute Mean Regression Coefficients
include("generate/mean_regression.jl")
# Export covariances for optimization
include("generate/export_covariances.jl")
### at this point, switch to python/JAX to do the covariance optimization and save regression coefficients ###
# Create pattern scaling baseline
include("generate/pattern_scaling.jl")