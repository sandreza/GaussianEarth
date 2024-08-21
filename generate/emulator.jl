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
using Statistics, Random
import Statistics.mean, Statistics.std, Random.rand

struct GaussianEmulator{M, L, B, A1, A2}
    mean::M
    decomposition::L
    basis::B
    month::A1
    global_mean_temperature::A2
end

function GaussianEmulator(data_directory; modes = 1000)
    hfile = h5open(data_directory * "_gaussian_model.hdf5", "r")
    mean = read(hfile["mean"])[1:modes, :, :]
    decomposition = read(hfile["L model"])[1:modes, 1:modes, :, :]
    basis = read(hfile["basis"])[:, 1:modes]
    close(hfile)
    return GaussianEmulator(mean, decomposition, basis, [1], [1.0555532f0])
end

function GaussianEmulator(μmodel, Lmodel, basis)
    return GaussianEmulator(μmodel, Lmodel, basis, [1], [1.0555532f0])
end

emulator = GaussianEmulator(μmodel, Lmodel, basis)

function mean(emulator::GaussianEmulator; modes = 1000)
    M, N = size(emulator.basis)
    month = emulator.month[1]
    global_mean_temperature = emulator.global_mean_temperature[1]
    ensemble_mean = zeros(Float32, M)
    for i in 1:modes
        ensemble_mean .+= emulator.basis[:, i] *(emulator.mean[i, 1, month] + emulator.mean[i, 2, month]* global_mean_temperature)
    end
    return ensemble_mean
end

function mode_mean(emulator; modes = 1000)
    global_mean_temperature = emulator.global_mean_temperature[1]
    month = emulator.month[1]
    return emulator.mean[1:modes, 1, month] + emulator.mean[1:modes, 2, month]* global_mean_temperature
end

function mode_variance(emulator; modes = 1000)
    global_mean_temperature = emulator.global_mean_temperature[1]
    month = emulator.month[1]
    L = emulator.decomposition[:, :, 1, month] + emulator.decomposition[:, :, 2, month]*global_mean_temperature
    Σ = (L' * L)[1:modes, 1:modes]
    return [Σ[i,i] for i in 1:modes]
end

function variance(emulator::GaussianEmulator; modes = 1000)
    M, N = size(emulator.basis)
    N = minimum([N, modes])
    global_mean_temperature = emulator.global_mean_temperature[1]
    month = emulator.month[1]
    L = emulator.decomposition[:, :, 1, month] + emulator.decomposition[:, :, 2, month]*global_mean_temperature
    Σ = (L' * L)[1:modes, 1:modes]
    variance_field = zeros(Float32, M)
    for i in ProgressBar(1:M)
        v = emulator.basis[i, 1:N]
        variance_field[i] = v' * (Σ * v)
    end
    return variance_field
end

function emulator_variance(emulator::GaussianEmulator; modes = 1000)
    global_mean_temperature = emulator.global_mean_temperature[1]
    month = emulator.month[1]
    L = emulator.decomposition[:, :, 1, month] + emulator.decomposition[:, :, 2, month]*global_mean_temperature
    Σ = (L' * L)[1:modes, 1:modes]
    return Σ
end

mean_field = mean(emulator)
variance_field = variance(emulator)
mean_modes = mode_mean(emulator)
variance_modes = mode_variance(emulator)