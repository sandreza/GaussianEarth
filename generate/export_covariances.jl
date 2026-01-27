scenario_directories = readdir(data_directory)
current_path = joinpath(data_directory, scenario_directories[1])
variable_directories = readdir(current_path)

##
field = "tas"
modes, temperature = concatenate_regression(field, ["historical", "ssp585"])
##
regression_indices = [1:139..., 141:251...]
covsave = h5open(save_directory * field * "_covariances.hdf5", "w")
covsave["temperature"] = temperature[regression_indices]
covsave["scale"] = 273

inds = 1:30
Σs = zeros(Float32, 1000, 1000, length(1:12:size(modes,2)))
for month in ProgressBar(1:12)
    for (i,j) in ProgressBar(enumerate(month:12:size(modes,2)))
        Σs[:,:,i] .= cov(modes[1:1000,j, inds], dims = 2)
    end
    covsave["covariances $month"] = Σs[:,:, regression_indices]
end

close(covsave)

##
field = "hurs"
modes, temperature = concatenate_regression(field, ["historical", "ssp585"])
##
regression_indices = [1:139..., 141:251...]
covsave = h5open(save_directory * field * "_covariances.hdf5", "w")
covsave["temperature"] = temperature[regression_indices]
covsave["scale"] = 273

inds = 1:30
Σs = zeros(Float32, 1000, 1000, length(1:12:size(modes,2)))
for month in ProgressBar(1:12)
    for (i,j) in ProgressBar(enumerate(month:12:size(modes,2)))
        Σs[:,:,i] .= cov(modes[1:1000,j, inds], dims = 2)
    end
    covsave["covariances $month"] = Σs[:,:, regression_indices]
end

close(covsave)
