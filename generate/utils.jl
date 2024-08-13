using NCDatasets, LinearAlgebra, Statistics, HDF5, ProgressBars

# grabbing data
function regression_variable(scenario; directory = "/net/fs06/d3/sandre/GaussianEarthData", normalized = true)
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

function eof_modes(variable, scenario; directory = "/net/fs06/d3/sandre/GaussianEarthData")
    file_name = "/$(variable)_$(scenario)_projection.hdf5"
    hfile = h5open(directory * file_name, "r")
    projection = read(hfile["projection"]) # has all 12 months
    close(hfile)
    return projection
end

function log_eof_modes(variable, scenario; directory = "/net/fs06/d3/sandre/GaussianEarthData")
    file_name = "/$(variable)_$(scenario)_log_projection.hdf5"
    hfile = h5open(directory * file_name, "r")
    projection = read(hfile["projection"]) # has all 12 months
    close(hfile)
    return projection
end

function metric_eof_modes(variable, scenario; directory = "/net/fs06/d3/sandre/GaussianEarthData")
    file_name = "/$(variable)_$(scenario)_metric_projection.hdf5"
    hfile = h5open(directory * file_name, "r")
    projection = read(hfile["projection"]) # has all 12 months
    close(hfile)
    return projection
end

function concatenate_regression(variable, scenarios; directory = "/net/fs06/d3/sandre/GaussianEarthData", normalized = true)
    modes = eof_modes(variable, scenarios[1]; directory)
    temperature = regression_variable(scenarios[1]; directory, normalized)
    for scenario in scenarios[2:end]
        modes = cat(modes, eof_modes(variable, scenario), dims = 2)
        temperature = cat(temperature, regression_variable(scenario), dims = 1)
    end
    return modes, temperature
end

function concatenate_log_regression(variable, scenarios; directory = "/net/fs06/d3/sandre/GaussianEarthData", normalized = true)
    modes = log_eof_modes(variable, scenarios[1]; directory)
    temperature = regression_variable(scenarios[1]; directory, normalized)
    for scenario in scenarios[2:end]
        modes = cat(modes, log_eof_modes(variable, scenario), dims = 2)
        temperature = cat(temperature, regression_variable(scenario), dims = 1)
    end
    return modes, temperature
end

function concatenate_metric_regression(variable, scenarios; directory = "/net/fs06/d3/sandre/GaussianEarthData", normalized = true)
    modes = metric_eof_modes(variable, scenarios[1]; directory)
    temperature = regression_variable(scenarios[1]; directory, normalized)
    for scenario in scenarios[2:end]
        modes = cat(modes, metric_eof_modes(variable, scenario), dims = 2)
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
