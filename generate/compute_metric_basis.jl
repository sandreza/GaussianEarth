using NCDatasets, LinearAlgebra, Statistics, HDF5, ProgressBars

save_directory = "/net/fs06/d3/sandre/GaussianEarthData/"
data_directory = "/net/fs06/d3/mgeo/CMIP6/interim/"

scenario_directories = readdir(data_directory)
scenario_directory = scenario_directories[1]

current_path = joinpath(data_directory, scenario_directory)
variable_directories = readdir(current_path)

@info "Computing basis for all variables"
for variable_directory in ProgressBar(variable_directories)
    local_current_path = joinpath(current_path, variable_directory)
    file_names = readdir(local_current_path)

    file_name = file_names[1] # pick the first file for computing a basis
    file_path = joinpath(local_current_path, file_name)

    ds = Dataset(file_path)
    field_name = keys(ds)[end] # always the last key for variable of interest
    field = Float32.(ds[field_name][:,:,:])
    M, N, L = size(field)

    latitude = Float32.(ds["lat"][:])
    longitude = Float32.(ds["lon"][:])
    Δϕ = reshape(2π / M * ones(M), (M, 1, 1))
    Δθ = reshape(π / N * ones(N) .* cos.(deg2rad.(latitude)), (1, N, 1))
    metric = (Δθ .* Δϕ) / (4π)

    metric_weighted_field = sqrt.(metric) .* field
    reshaped_field = reshape(field, M * N, L) .+ eps(eltype(field)(1.0))
    U, S, V = svd(reshaped_field)

    basis = U
    singular_values = S



    hfile = h5open(save_directory * field_name * "_metric_basis.hdf5", "w")
    hfile["basis"] = basis
    hfile["singular_values"] = singular_values
    hfile["latitude"] = latitude
    hfile["longitude"] = longitude
    hfile["metric"] = metric
    close(hfile)
end