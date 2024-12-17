using CairoMakie
using NCDatasets, LinearAlgebra, Statistics, HDF5, ProgressBars
using LinearAlgebra
using Printf

ts = 60
xls = 70 
yls = 70
tls = 70
legend_ls = 60
resolution = (400, 150) .* 12
lw = 7
global_common_options = (; titlesize = ts, xlabelsize = xls, ylabelsize = yls, xticklabelsize = tls, yticklabelsize = tls)

if process_data
    include("emulator.jl")
    include("emulator_hurs.jl")
    _, temperatures = concatenate_regression("tas", ["historical", "ssp585"])
end
##
month = 1
field = "tas"

hfile = h5open(save_directory * field * "_basis.hdf5", "r")
latitude = read(hfile["latitude"])
longitude = read(hfile["longitude"])
metric = read(hfile["metric"])
close(hfile)
fmetric = reshape(metric, (192*96, 1))
##
function global_mean_value(field)
    global_mean_value = sum(field .* fmetric, dims = 1)[1]
    return global_mean_value
end
function global_mean_upper(field)
    global_mean_value_upper = sum(reshape(field .* fmetric, (192, 96))[:, 49:end, :], dims = (1, 2)) * 2
    return global_mean_value_upper[1]
end
function global_mean_lower(field)
    global_mean_value_lower = sum(reshape(field .* fmetric, (192, 96))[:, 1:48, :], dims = (1, 2)) * 2
    return global_mean_value_lower[1]
end

index_1 = [157, 46]
index_2 = [54, 48]
index_3 = [19, 87]

function location_1(field)
    rfield = reshape(field, (192, 96))
    return rfield[157, 46]
end

function location_2(field)
    rfield = reshape(field, (192, 96))
    return rfield[54, 48]
end

function location_3(field)
    rfield = reshape(field, (192, 96))
    return rfield[19, 87]
end

##
min_temp = minimum(temperatures) # 1.0544952f0
max_temp = maximum(temperatures) # 1.0720718f0
month = 1
##
observables = [location_1, location_3]
##
if process_data
    historical_tas = common_array("historical", "tas"; ensemble_members = 45)
    historical_hurs = common_array("historical", "hurs"; ensemble_members = 29)
    ssp585_tas = common_array("ssp585", "tas"; ensemble_members = 45)
    ssp585_hurs = common_array("ssp585", "hurs"; ensemble_members = 29)
    historical_temperatures = regression_variable("historical")
    ssp585_temperatures = regression_variable("ssp585")
end
##
indmin = 100
min_temp = historical_temperatures[indmin]
indmax = 43
max_temp = ssp585_temperatures[indmax]
##
window = 2 
increment = 2 * window + 1
month = 1
fig = Figure(; resolution)
observable_names_base = [@sprintf("%.2fᵒ, %.2fᵒ ", longitude[tmp[1]], latitude[tmp[2]]) for tmp in [index_1, index_2, index_3]]
observable_names = observable_names_base .* ["(Jan)"]
for i in eachindex(observables)
    emulator.month[1] = month
    emulator.global_mean_temperature[1] = min_temp
    emulator_hurs.month[1] = month
    emulator_hurs.global_mean_temperature[1] = min_temp
    tasmean, tasstd = emulator_mean_variance_linear_functionals(observables, emulator)
    hursmean, hursstd = emulator_mean_variance_linear_functionals(observables, emulator_hurs)

    if i == 1
        common_options_1 = (; xlabel = "Temperature (K)", ylabel = "PDF")
        common_options_2 = (; xlabel = "Relative Humidity (%)", ylabel = "PDF")
    else
        common_options_1 = (; xlabel = "Temperature (K)")
        common_options_2 = (; xlabel = "Relative Humidity (%)")
    end
    ax1 = Axis(fig[1, i]; title = observable_names[i], common_options_1..., global_common_options...)
    rfield = reshape(historical_tas[:, :, indmin-window:indmin+window, month, :], (192 * 96, 45 * increment))
    tas_hist = [observables[i](rfield[:, j]) for j in 1:45*increment]
    xs, ys = gaussian_grid(tasmean[i], tasstd[i])
    hist!(ax1, tas_hist, bins = 20, color = (:royalblue, 0.2), normalization = :pdf)
    lines!(ax1, xs, ys, color = :blue, linewidth = lw)
    ax2 = Axis(fig[2, i]; title = observable_names[i], common_options_2..., global_common_options...)
    rfield = reshape(historical_hurs[:, :, indmin-window:indmin+window, month, :], (192 * 96, 29*increment))
    hurs_hist = [observables[i](rfield[:, j]) for j in 1:29*increment]
    hist!(ax2, hurs_hist, bins = 20, color = (:royalblue, 0.2), normalization = :pdf, label = "1950 (Data)")
    xs, ys = gaussian_grid(hursmean[i], hursstd[i])
    lines!(ax2, xs, ys, color = :blue, label = "1950 (Emulator)", linewidth = lw)

    emulator.month[1] = month
    emulator.global_mean_temperature[1] = max_temp
    emulator_hurs.month[1] = month
    emulator_hurs.global_mean_temperature[1] = max_temp
    tasmean, tasstd = emulator_mean_variance_linear_functionals(observables, emulator)
    hursmean, hursstd = emulator_mean_variance_linear_functionals(observables, emulator_hurs)

    xs, ys = gaussian_grid(tasmean[i], tasstd[i])
    lines!(ax1, xs, ys, color = :orange, linewidth = lw)
    xs, ys = gaussian_grid(hursmean[i], hursstd[i])
    lines!(ax2, xs, ys, color = :orange, label = "2050 (Emulator)", linewidth = lw)

    rfield = reshape(ssp585_tas[:, :, indmax-window:indmax+window, month, :], (192 * 96, 45*increment))
    tas_hist = [observables[i](rfield[:, j]) for j in 1:45*increment]
    rfield = reshape(ssp585_hurs[:, :, indmax-window:indmax+window, month, :], (192 * 96, 29*increment))
    hurs_hist = [observables[i](rfield[:, j]) for j in 1:29*increment]
    hist!(ax1, tas_hist, bins = 20, color = (:orangered2, 0.2), normalization = :pdf, label = "2050 (Data)")
    hist!(ax2, hurs_hist, bins = 20, color = (:orangered2, 0.2), normalization = :pdf, label = "2050 (Data)")
    if i == 1
        axislegend(ax2, position = :lt, labelsize = legend_ls)
    end
end

month = 7
observable_names_base = [@sprintf("%.2fᵒ, %.2fᵒ ", longitude[tmp[1]], latitude[tmp[2]]) for tmp in [index_1, index_2, index_3]]
observable_names = observable_names_base .* ["(Jul)"]
for i in eachindex(observables)
    emulator.month[1] = month
    emulator.global_mean_temperature[1] = min_temp
    emulator_hurs.month[1] = month
    emulator_hurs.global_mean_temperature[1] = min_temp
    tasmean, tasstd = emulator_mean_variance_linear_functionals(observables, emulator)
    hursmean, hursstd = emulator_mean_variance_linear_functionals(observables, emulator_hurs)

    common_options_1 = (; xlabel = "Temperature (K)")
    common_options_2 = (; xlabel = "Relative Humidity (%)")

    ax1 = Axis(fig[1, i+2]; title = observable_names[i], common_options_1..., global_common_options...)
    rfield = reshape(historical_tas[:, :, indmin-window:indmin+window, month, :], (192 * 96, 45 * increment))
    tas_hist = [observables[i](rfield[:, j]) for j in 1:45*increment]
    xs, ys = gaussian_grid(tasmean[i], tasstd[i])
    hist!(ax1, tas_hist, bins = 20, color = (:royalblue, 0.2), normalization = :pdf)
    lines!(ax1, xs, ys, color = :blue, linewidth = lw)
    ax2 = Axis(fig[2, i+2]; title = observable_names[i], common_options_2..., global_common_options...)
    rfield = reshape(historical_hurs[:, :, indmin-window:indmin+window, month, :], (192 * 96, 29*increment))
    hurs_hist = [observables[i](rfield[:, j]) for j in 1:29*increment]
    hist!(ax2, hurs_hist, bins = 20, color = (:royalblue, 0.2), normalization = :pdf)
    xs, ys = gaussian_grid(hursmean[i], hursstd[i])
    lines!(ax2, xs, ys, color = :blue, linewidth = lw)

    emulator.month[1] = month
    emulator.global_mean_temperature[1] = max_temp
    emulator_hurs.month[1] = month
    emulator_hurs.global_mean_temperature[1] = max_temp
    tasmean, tasstd = emulator_mean_variance_linear_functionals(observables, emulator)
    hursmean, hursstd = emulator_mean_variance_linear_functionals(observables, emulator_hurs)

    xs, ys = gaussian_grid(tasmean[i], tasstd[i])
    lines!(ax1, xs, ys, color = :orange, linewidth = lw)
    xs, ys = gaussian_grid(hursmean[i], hursstd[i])
    lines!(ax2, xs, ys, color = :orange, linewidth = lw)

    rfield = reshape(ssp585_tas[:, :, indmax-window:indmax+window, month, :], (192 * 96, 45*increment))
    tas_hist = [observables[i](rfield[:, j]) for j in 1:45*increment]
    rfield = reshape(ssp585_hurs[:, :, indmax-window:indmax+window, month, :], (192 * 96, 29*increment))
    hurs_hist = [observables[i](rfield[:, j]) for j in 1:29*increment]
    hist!(ax1, tas_hist, bins = 20, color = (:orangered2, 0.2), normalization = :pdf)
    hist!(ax2, hurs_hist, bins = 20, color = (:orangered2, 0.2), normalization = :pdf)
end

display(fig)

save(figure_directory * "climate_change_shifts_hurs_tas_with_data_jan_july_points.png", fig)