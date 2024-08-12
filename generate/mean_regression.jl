using NCDatasets, LinearAlgebra, Statistics, HDF5, ProgressBars

include("utils.jl")
##
save_directory = "/net/fs06/d3/sandre/GaussianEarthData/"
data_directory = "/net/fs06/d3/mgeo/CMIP6/interim/"
scenario_directories = readdir(data_directory)
current_path = joinpath(data_directory, scenario_directories[1])
variable_directories = readdir(current_path)

##
field = "tas"
modes, temperature = concatenate_regression(field, ["historical", "ssp585"])

month = 1
order = 1
regression_coefficients = Float32.(regression(modes, temperature, month; order))

i = 1 # pick a mode
j = 1 # pick a year
rc = regression_coefficients[i, :]
model = sum([rc[k+1]* (temperature[j])^k for k in 0:order])
truth = mean(modes[i, month:12:end, :], dims = 2)[j]
relative_error_percent = abs.((model - truth) / truth) * 100
