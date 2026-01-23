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

include("figure_1.jl") # gmt in mpi ensemble ✅
include("figure_2.jl") # parametric assumption ✅
include("figure_3.jl") # regression check ✅
include("figure_4.jl") # emulator error maps ✅
include("figure_5.jl") # std maps ✅
include("figure_6.jl") # zonal average + location histograms ✅
include("figure_7.jl") # histograms ssp245 land/ocean ✅
include("figure_8.jl") # histograms ssp245 locations ✅
include("generate_new_scenario.jl") ## MOVE TO GENERATE ✅
include("figures_case_study.jl")
include("figure_A1.jl") # check model fits
include("figure_B1.jl")
include("figure_B2.jl") # linear emulator error maps ✅
include("figure_C1.jl") # linear mean error  ✅
include("figure_C2.jl") # regression std error
include("figure_D1.jl") # global + hemisphere seasonal ✅ #relies on pre-loaded "emulator" --> should probably standardize to not this, but performance?
include("figures_D45.jl") # emulated vs mpi comparison ✅