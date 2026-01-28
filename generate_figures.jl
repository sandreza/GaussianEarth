figure_directory = "figures/"
if !isdir(figure_directory)
    mkpath(figure_directory)
end

using CairoMakie, Printf, GeoMakie
using NCDatasets, LinearAlgebra, Statistics, HDF5, ProgressBars
using LinearAlgebra, Distributions, Random

# save_directory = "PLEASE/SET/YOUR/SAVE/PATH/HERE/"
# data_directory = "PLEASE/SET/YOUR/DATA/PATH/HERE/"
save_directory = "/net/fs06/d3/mgeo/GaussianEarthData/"
# data_directory = "/net/fs06/d3/mgeo/CMIP6/interim/"

# scenario_directories = readdir(data_directory) ## this stuff should only be called if no ground truth bundle
# current_path = joinpath(data_directory, scenario_directories[1])
# variable_directories = readdir(current_path)

include("utils.jl")

use_ground_truth_bundle = has_ground_truth_bundle()
ground_truth_bundle_path = get_ground_truth_bundle_path()

include("paper_plots/figure_1.jl") # gmt in mpi ensemble
include("paper_plots/figure_2.jl") # parametric assumption check
include("paper_plots/figure_3.jl") # regression check 
include("paper_plots/figure_4.jl") # emulator error maps
include("paper_plots/figure_5.jl") # std maps 
include("paper_plots/figure_6.jl") # zonal average + location histograms
include("paper_plots/figure_7.jl") # histograms ssp245 land/ocean and figure D2 (ssp585)
include("paper_plots/figure_8.jl") # histograms ssp245 locations and figure D3 (ssp585)
include("paper_plots/generate_new_scenario.jl") # generate scenario for case study
include("paper_plots/figures_case_study.jl") # case study figures
include("paper_plots/figure_A1.jl") # check model fits
include("paper_plots/figure_B1.jl") # linear vs quadratic errors ### STILL VERY MESSY
include("paper_plots/figure_B2.jl") # linear emulator error maps 
include("paper_plots/figure_C1.jl") # linear mean error
include("paper_plots/figure_C2.jl") # regression std error
include("paper_plots/figure_D1.jl") # global + hemisphere seasonal
include("paper_plots/figures_D45.jl") # emulated vs mpi comparison ### WHICH FIGS, WHICH SCENARIOS?
