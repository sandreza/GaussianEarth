using NCDatasets, LinearAlgebra, Statistics, HDF5, ProgressBars

include("utils.jl")
##
save_directory = "/net/fs06/d3/sandre/GaussianEarthData/"
data_directory = "/net/fs06/d3/mgeo/CMIP6/interim/"
scenario_directories = readdir(data_directory)
current_path = joinpath(data_directory, scenario_directories[1])
variable_directories = readdir(current_path)

##
for field in ProgressBar(["tas"])
    modes, temperature = concatenate_regression(field, ["historical", "ssp585"])
    order = 1
    hfile = h5open(save_directory * field * "_mean_regression.hdf5", "w") 
    for month in ProgressBar(1:12)
        regression_coefficients = Float32.(regression(modes, temperature, month; order))
        hfile["regression_coefficients $month"] = regression_coefficients
    end
    close(hfile)
end

for field in ProgressBar(["hurs"])
    modes, temperature = concatenate_regression(field, ["historical", "ssp585"])
    ensemble_members = 1:30
    order = 1
    hfile = h5open(save_directory * field * "_mean_regression.hdf5", "w") 
    for month in ProgressBar(1:12)
        regression_coefficients = Float32.(regression(modes, temperature, month; order, ensemble_members))
        hfile["regression_coefficients $month"] = regression_coefficients
    end
    close(hfile)
end
