using CairoMakie
using NCDatasets, LinearAlgebra, Statistics, HDF5, ProgressBars
using LinearAlgebra

include("utils.jl")
##
save_directory = "/net/fs06/d3/sandre/GaussianEarthData/"
data_directory = "/net/fs06/d3/mgeo/CMIP6/interim/"
scenario_directories = readdir(data_directory)
current_path = joinpath(data_directory, scenario_directories[1])
variable_directories = readdir(current_path)

##
field = "tas"
covsave = h5open(save_directory * field * "_covariances_model.hdf5", "r")
Ls = []
for month in 1:12
    L = read(covsave["L" * "$month"])
    push!(Ls, L)
end
close(covsave)

rC = []
hfile = h5open(save_directory * field * "_mean_regression.hdf5", "r") 
for month in 1:12
    regression_coefficients = read(hfile["regression_coefficients $month"])
    push!(rC, regression_coefficients)
end
close(hfile)

hfile = h5open(save_directory * field * "_basis.hdf5", "r")
basis = read(hfile["basis"] )
close(hfile)

hfile = h5open(save_directory * field * "_gaussian_model.hdf5", "w")
Lmodel = zeros(Float32, size(Ls[1])..., 12)
[Lmodel[:, :, :, i] .= Ls[i] for i in 1:12]
μmodel = zeros(Float32, size(rC[1])..., 12)
[μmodel[:, :, i] .= rC[i] for i in 1:12]
hfile["mean"] = μmodel
hfile["L model"] = Lmodel
hfile["basis"] = basis
hfile["scale factor"] = 273 
close(hfile)

##
emulator = GaussianEmulator(μmodel, Lmodel, basis)
mean_field = mean(emulator)
variance_field = variance(emulator)
mean_modes = mode_mean(emulator)
variance_modes = mode_variance(emulator)