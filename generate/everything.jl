# Determine EOFs in the Historical Period for 1 Realization
include("compute_basis.jl")
# Project onto the EOFs
include("project_onto_basis.jl")
# Compute Regression Variable (global, ensemble, and yearly averaged surface temperature)
include("regression_variable.jl")
# Compute Mean Regression Coefficients
include("mean_regression.jl")
# Export covariances for optimization
include("export_covariances.jl")

# Create pattern scaling baseline
include("pattern_scaling.jl")