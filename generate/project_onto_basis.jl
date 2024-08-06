using NCDatasets, LinearAlgebra, Statistics, HDF5, ProgressBars

save_directory = "/net/fs06/d3/sandre/GaussianEarthData/"
data_directory = "/net/fs06/d3/mgeo/CMIP6/interim/"

scenario_directories = readdir(data_directory)
@info "Computing Projection for all variables and all scenarios"
for scenario in ProgressBar(scenario_directories[1:end])

    current_path = joinpath(data_directory, scenario)
    variable_directories = readdir(current_path)

    for variable_directory in ProgressBar(variable_directories)
        local_current_path = joinpath(current_path, variable_directory)
        file_names = readdir(local_current_path)

        file_name = file_names[1] # pick the first file for obtaining varaibles
        file_path = joinpath(local_current_path, file_name)

        ds = Dataset(file_path)
        field_name = keys(ds)[end]
        field = Float32.(ds[field_name][:,:,:])
        M, N, L = size(field)
        close(ds)

        hfile = h5open(save_directory * field_name * "_basis.hdf5", "r")
        basis = read(hfile["basis"])
        metric = read(hfile["metric"])
        close(hfile)

        projection = zeros(Float32, size(basis)[2], L, length(file_names))
        global_mean = zeros(Float32, L, length(file_names))

        for (i, file_name) in ProgressBar(enumerate(file_names))
            file_path = joinpath(local_current_path, file_name)
            ds = Dataset(file_path)
            field = Float32.(ds[field_name][:,:,:])
            M, N, L = size(field)
            reshaped_field = reshape(field, M * N, L)
            projected_field = basis' * reshaped_field
            projection[:, :, i] .= projected_field
            global_mean[:, i] .= sum(metric .* field, dims = (1, 2))[:]
        end

        hfile = h5open(save_directory * field_name * "_" * scenario * "_projection.hdf5", "w")
        hfile["projection"] = projection
        hfile["global mean"] = global_mean
        close(hfile)
    end
end