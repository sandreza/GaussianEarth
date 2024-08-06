using NCDatasets, LinearAlgebra, Statistics

data_directory = "/net/fs06/d3/mgeo/CMIP6/interim/"
scenario_directories = readdir(data_directory)
scenario_directory = scenario_directories[1]

current_path = joinpath(data_directory, scenario_directory)
variable_directories = readdir(current_path)
variable_directory = variable_directories[end]

current_path = joinpath(current_path, variable_directory)
file_names = readdir(current_path)

file_name = file_names[1]
file_path = joinpath(current_path, file_name)

ds = Dataset(file_path)
field_name = keys(ds)[end]
field = ds[field_name][:,:,:]
M, N, L = size(field)
reshaped_field = reshape(field, M * N, L)
U, S, V = svd(reshaped_field)