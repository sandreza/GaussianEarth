figure_directory = "figures/"
process_data = true

using CairoMakie, Printf, GeoMakie
using NCDatasets, LinearAlgebra, Statistics, HDF5, ProgressBars
using LinearAlgebra, Distributions

# save_directory = "/net/fs06/d3/sandre/GaussianEarthData/"
save_directory = "/net/fs06/d3/mgeo/GaussianEarthData/"
data_directory = "/net/fs06/d3/mgeo/CMIP6/interim/"
scenario_directories = readdir(data_directory)
current_path = joinpath(data_directory, scenario_directories[1])
variable_directories = readdir(current_path)

include("utils.jl")

#load in saved out tas emulator
field = "tas"
hfile = h5open(save_directory * field * "_gaussian_model.hdf5", "r")
μmodel = read(hfile["mean"])
μmodel_quadratic = read(hfile["mean_quadratic"])
Lmodel = read(hfile["L model"])
basis = read(hfile["basis"])
close(hfile)
emulator = GaussianEmulator(μmodel_quadratic, Lmodel, basis)
mean_field = mean(emulator)
variance_field = variance(emulator)
mean_modes = mode_mean(emulator)
variance_modes = mode_variance(emulator)

include("figure_1.jl")
include("figure_2.jl")
include("figure_3.jl")
include("figure_4.jl")
include("figure_5.jl")
include("figure_6.jl")
include("figure_7.jl")
include("figure_8.jl")
include("figure_9_land_sea_redux.jl")
include("figure_10_redux.jl")
include("parametric_assumption_redux.jl")
include("generate_new_scenario_redux.jl")
include("case_study_redux3.jl")