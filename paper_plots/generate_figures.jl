figure_directory = "figures/"
process_data = false

using CairoMakie, Printf, GeoMakie
using NCDatasets, LinearAlgebra, Statistics, HDF5, ProgressBars
using LinearAlgebra, Distributions

save_directory = "/net/fs06/d3/sandre/GaussianEarthData/"
data_directory = "/net/fs06/d3/mgeo/CMIP6/interim/"
scenario_directories = readdir(data_directory)
current_path = joinpath(data_directory, scenario_directories[1])
variable_directories = readdir(current_path)



# include("figure_2.jl")
# include("figure_3.jl")
# include("figure_4.jl")
include("figure_5.jl")
# include("figure_6.jl")
# include("figure_7.jl")
# include("figure_8.jl")
# include("figure_9_land_sea_redux.jl")
# include("figure_10_redux.jl")
# include("parametric_assumption_redux.jl")