using NCDatasets, LinearAlgebra, Statistics, HDF5, ProgressBars
using Random
import Statistics.mean, Statistics.std, Random.rand


const DEFAULT_SAVE_DIRECTORY = "PLEASE/SET/YOUR/SAVE/PATH/HERE/"
const DEFAULT_DATA_DIRECTORY = "PLEASE/SET/YOUR/DATA/PATH/HERE/"
const DEFAULT_GROUND_TRUTH_BUNDLE = "ground_truth_bundle.hdf5"

# Compatibility helpers for HDF5.jl (precompute_ground_truth_bundle.jl expects these names)
g_create(parent, name; kwargs...) = create_group(parent, name; kwargs...)
d_create(parent, name, dtype, dspace; kwargs...) = create_dataset(parent, name, dtype, dspace; kwargs...)

function get_save_directory()
    return @isdefined(save_directory) ? save_directory : DEFAULT_SAVE_DIRECTORY
end
function get_data_directory()
    return @isdefined(data_directory) ? data_directory : DEFAULT_DATA_DIRECTORY
end

function get_ground_truth_bundle_path()
    return joinpath(get_save_directory(), DEFAULT_GROUND_TRUTH_BUNDLE)
end

function has_ground_truth_bundle(bundle_path = get_ground_truth_bundle_path())
    return ispath(bundle_path)
end

function read_ground_truth(dataset; bundle_path = get_ground_truth_bundle_path())
    h5open(bundle_path, "r") do h5
        return read(h5[dataset])
    end
end

# grabbing data
function regression_variable(scenario; directory = get_save_directory(), normalized = true)
    file_name = "/tas_ensemble_yearly_average.hdf5"
    hfile = h5open(directory * file_name, "r")
    if normalized
        global_mean_temperature = read(hfile[scenario * "_normalized"])
    else
        global_mean_temperature = read(hfile[scenario])
    end
    close(hfile)
    return global_mean_temperature
end

function eof_modes(variable, scenario; directory = get_save_directory())
    file_name = "/$(variable)_$(scenario)_projection.hdf5"
    hfile = h5open(directory * file_name, "r")
    projection = read(hfile["projection"]) # has all 12 months
    close(hfile)
    return projection
end

function eof_basis(variable; directory = get_save_directory())
    file_name = "/$(variable)_basis.hdf5"
    hfile = h5open(directory * file_name, "r")
    basis = read(hfile["basis"])
    close(hfile)
    return basis
end

function concatenate_regression(variable, scenarios; directory = get_save_directory(), normalized = true)
    modes = eof_modes(variable, scenarios[1]; directory)
    temperature = regression_variable(scenarios[1]; directory, normalized)
    for scenario in scenarios[2:end]
        modes = cat(modes, eof_modes(variable, scenario), dims = 2)
        temperature = cat(temperature, regression_variable(scenario), dims = 1)
    end
    return modes, temperature
end

function regression(modes, temperature, month; order = 1, ensemble_members = 1:40)
    # horizontal axis
    order > 5 ? error("Order must be less than 5") : nothing
    rX = ones(Float64, size(modes[1, 1:12:end, ensemble_members])..., 1 + order)
    for i in 1:size(modes[1, 1:12:end, :])[1]
        for j in 1:order
            rX[i, :, 1+j] .= Float64(temperature[i])^j
        end
    end
    number_of_modes = size(modes, 1)
    regression_coefficients = zeros(Float64, number_of_modes, 1 + order)
    for i in 1:number_of_modes
        rY = Float64.(modes[i, month:12:end, ensemble_members])
        regression_coefficients[i, :] .= reshape(rX, length(rY), order + 1) \ reshape(rY, length(rY))
    end
    return regression_coefficients
end

function common_array(scenario, variable; data_directory = get_data_directory(), ensemble_members = 45)
    current_path = joinpath(data_directory, scenario)
    local_current_path = joinpath(current_path, variable)
    file_names = readdir(local_current_path)
    fields = Array{Float32, 3}[]
    if length(file_names) > 0 # sometimes the directory is empty
        file_name = file_names[1] # pick the first file for obtaining varaibles
        file_path = joinpath(local_current_path, file_name)
        for (i, file_name) in ProgressBar(enumerate(file_names[1:ensemble_members]))
            file_path = joinpath(local_current_path, file_name)
            ds = Dataset(file_path)
            field = Float32.(ds[variable][:,:,:])
            push!(fields, field)
        end
    end
    time1 = size(fields[1])[end]
    monthtime1 = time1 ÷ 12
    all_together = zeros(Float32, size(fields[1])[1], size(fields[1])[2], monthtime1, 12, ensemble_members);
    for ω in ProgressBar(1:ensemble_members)
        for month in 1:12
            @inbounds all_together[:,:, 1:monthtime1, month, ω] .= fields[ω][:, :, month:12:end]
        end
    end
    return all_together
end

function ensemble_averaging(
    scenario,
    variable;
    data_directory = get_data_directory(),
    ensemble_members = 45,
    return_std = false,
    use_bundle = true,
    bundle_path = get_ground_truth_bundle_path()
) 
    if use_bundle && ispath(bundle_path) && variable == "tas"
        try
            if return_std
                std_monthly = read_ground_truth("tas/std_monthly/$scenario"; bundle_path)
                std_mean = mean(std_monthly, dims = 4)[:,:,:,1]
                return std_monthly, std_mean
            else
                mean_monthly = read_ground_truth("tas/mean_monthly/$scenario"; bundle_path)
                annual_mean = read_ground_truth("tas/annual_mean/$scenario"; bundle_path)
                return mean_monthly, annual_mean
            end
        catch err
            @warn "Ground-truth bundle missing data; falling back to raw data" err
        end
    end

    current_path = joinpath(data_directory, scenario)
    local_current_path = joinpath(current_path, variable)
    file_names = readdir(local_current_path)
    fields = Array{Float32, 3}[]
    if length(file_names) > 0 # sometimes the directory is empty
        file_name = file_names[1] # pick the first file for obtaining varaibles
        file_path = joinpath(local_current_path, file_name)
        for (i, file_name) in ProgressBar(enumerate(file_names[1:ensemble_members]))
            file_path = joinpath(local_current_path, file_name)
            ds = Dataset(file_path)
            field = Float32.(ds[variable][:,:,:])
            push!(fields, field)
        end
    end
    time1 = size(fields[1])[end]
    monthtime1 = time1 ÷ 12
    all_together = zeros(Float32, size(fields[1])[1], size(fields[1])[2], monthtime1, 12, ensemble_members);
    for ω in ProgressBar(1:ensemble_members)
        for month in 1:12
            @inbounds all_together[:,:, 1:monthtime1, month, ω] .= fields[ω][:, :, month:12:end]
        end
    end
    if return_std
        # returns both the average std and the monthly stds
        return std(all_together, dims = 5)[:,:,:,:,1], mean(std(all_together, dims =4)[:,:,:,:, 1], dims=4)[:,:,:,1]
    else
        return mean(all_together, dims = 5)[:,:,:,:,1], mean(all_together, dims = (4, 5))[:,:,:,1, 1]
    end
end

# computing statistics 

function kurtosis(x)
    n = length(x)
    μ = mean(x)
    σ = std(x)
    return sum((x .- μ).^4) / (n * σ^4) - 3
end
function skewness(x)
    n = length(x)
    μ = mean(x)
    σ = std(x)
    return sum((x .- μ).^3) / (n * σ^3)
end

gaussian(x, μ, σ) = exp(-0.5 * ((x .- μ) ./ σ).^2) ./ sqrt(2π * σ^2)

function gaussian_grid(μ, σ; factor = 4, n = 100)
    xs = range(μ - factor * σ, stop = μ + factor * σ, length = n)
    ys = gaussian.(xs, μ, σ)
    return collect(xs), ys
end

##
struct CovarEmulator{M, L, B, A1, A2}
    mean::M
    decomposition::L
    basis::B
    month::A1
    global_mean_temperature::A2
end

function CovarEmulator(data_directory; modes = 1000)
    hfile = h5open(data_directory * "_model.hdf5", "r")
    mean = read(hfile["mean"])[1:modes, :, :]
    decomposition = read(hfile["L model"])[1:modes, 1:modes, :, :]
    basis = read(hfile["basis"])[:, 1:modes]
    close(hfile)
    return CovarEmulator(mean, decomposition, basis, [1], [1.0555532f0]) 
end

function CovarEmulator(μmodel, Lmodel, basis)
    return CovarEmulator(μmodel, Lmodel, basis, [1], [1.0555532f0])
end

function mean(emulator::CovarEmulator; modes = 1000)
    M, N = size(emulator.basis)
    k = size(emulator.mean)[2] - 1
    month = emulator.month[1]
    global_mean_temperature = emulator.global_mean_temperature[1]
    ensemble_mean = zeros(Float32, M)
    for i in 1:modes
        if k == 1
            ensemble_mean .+= emulator.basis[:, i] *(emulator.mean[i, 1, month] + emulator.mean[i, 2, month]* global_mean_temperature)
        elseif k == 2
            ensemble_mean .+= emulator.basis[:, i] *(emulator.mean[i, 1, month] + emulator.mean[i, 2, month]* global_mean_temperature + emulator.mean[i, 3, month]* global_mean_temperature^2)
        end
    end
    return ensemble_mean
end

function rand(emulator::CovarEmulator; modes = 1000)
    M, N = size(emulator.basis)
    month = emulator.month[1]
    k = size(emulator.mean)[2] - 1
    global_mean_temperature = emulator.global_mean_temperature[1]
    if k == 1
        μ = emulator.mean[1:modes, 1, month] .+ emulator.mean[1:modes, 2, month]* global_mean_temperature
    elseif k == 2
        μ = emulator.mean[1:modes, 1, month] .+ emulator.mean[1:modes, 2, month]* global_mean_temperature + emulator.mean[1:modes, 3, month]* global_mean_temperature^2
    end
    L = emulator.decomposition[:, :, 1, month] + emulator.decomposition[:, :, 2, month]*global_mean_temperature
    Σ = (L' * L)[1:modes, 1:modes]
    a = rand(MvNormal(μ, Σ))
    realization = zeros(Float32, M)
    for i in 1:modes
        realization .+= emulator.basis[:, i] * a[i]
    end
    return realization
end

function mode_mean(emulator; modes = 1000)
    global_mean_temperature = emulator.global_mean_temperature[1]
    month = emulator.month[1]
    k = size(emulator.mean)[2] - 1
    if k == 1
        return emulator.mean[1:modes, 1, month] + emulator.mean[1:modes, 2, month]* global_mean_temperature
    elseif k == 2
        return emulator.mean[1:modes, 1, month] + emulator.mean[1:modes, 2, month]* global_mean_temperature + emulator.mean[1:modes, 3, month]* global_mean_temperature^2
    end
end

function mode_variance(emulator; modes = 1000)
    global_mean_temperature = emulator.global_mean_temperature[1]
    month = emulator.month[1]
    L = emulator.decomposition[:, :, 1, month] + emulator.decomposition[:, :, 2, month]*global_mean_temperature
    Σ = (L' * L)[1:modes, 1:modes]
    return [Σ[i,i] for i in 1:modes]
end

function variance(emulator::CovarEmulator; modes = 1000)
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

function emulator_variance(emulator::CovarEmulator; modes = 1000, month = nothing, global_mean_temperature = nothing)
    if isnothing(global_mean_temperature)
        global_mean_temperature = emulator.global_mean_temperature[1]
    end
    if isnothing(month)
        month = emulator.month[1]
    end
    L = emulator.decomposition[:, :, 1, month] + emulator.decomposition[:, :, 2, month]*global_mean_temperature
    Σ = (L' * L)[1:modes, 1:modes]
    return Σ
end

function emulator_mean_variance_linear_functionals(observables, emulator::CovarEmulator; modes = 1000, show_progress = true)
    μs = zeros(Float32, length(observables))
    σs = zeros(Float32, length(observables))
    mean_modes = mode_mean(emulator; modes)
    Σ = emulator_variance(emulator; modes)
    for (i, observable) in (show_progress ? ProgressBar(enumerate(observables)) : enumerate(observables))
        observable_basis = [observable(emulator.basis[:, j]) for j in 1:modes]
        σ = sqrt(observable_basis' * (Σ * observable_basis)) 
        μ = mean_modes' * observable_basis
        μs[i] = Float32(μ)
        σs[i] = Float32(σ)
    end
    return μs, σs
end

function kl_div(mu1, sigma1, mu2, sigma2)
    # KL Divergence between N(mu1, sigma1^2) and N(mu2, sigma2^2)
    # measures the divergence of D1 from D2
    return log(sigma2 / sigma1) + (sigma1^2 + (mu1 - mu2)^2) / (2 * sigma2^2) - 0.5
end
