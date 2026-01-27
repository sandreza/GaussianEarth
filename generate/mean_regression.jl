scenario_directories = readdir(data_directory)
current_path = joinpath(data_directory, scenario_directories[1])
variable_directories = readdir(current_path)

## linear
for field in ProgressBar(["tas", "hurs"])
    modes, temperature = concatenate_regression(field, ["historical", "ssp585"])
    order = 1
    hfile = h5open(save_directory * field * "_mean_regression.hdf5", "w") 
    for month in ProgressBar(1:12)
        regression_coefficients = Float32.(regression(modes, temperature, month; order))
        hfile["regression_coefficients $month"] = regression_coefficients
    end
    close(hfile)
end

## quadratic
for field in ProgressBar(["tas", "hurs"])
    modes, temperature = concatenate_regression(field, ["historical", "ssp585"])
    order = 2
    hfile = h5open(save_directory * field * "_mean_regression_quadratic.hdf5", "w") 
    for month in ProgressBar(1:12)
        regression_coefficients = Float32.(regression(modes, temperature, month; order))
        hfile["regression_coefficients $month"] = regression_coefficients
    end
    close(hfile)
end