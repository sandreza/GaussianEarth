ts = 60
xls = 70 
yls = 70
tls = 70
legend_ls = 70
resolution = (400, 150) .* 12
lw = 7
global_common_options = (; titlesize = ts, xlabelsize = xls, ylabelsize = yls, xticklabelsize = tls, yticklabelsize = tls)

##

use_bundle = @isdefined(use_ground_truth_bundle) ? use_ground_truth_bundle : false
bundle_path = @isdefined(ground_truth_bundle_path) ? ground_truth_bundle_path : get_ground_truth_bundle_path()

include("../emulator.jl")
include("../emulator_hurs.jl")
_, temperatures = concatenate_regression("tas", ["ssp119"])

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
if use_bundle && has_ground_truth_bundle(bundle_path)
    function f2_samples(month, key, level)
        return read_ground_truth("samples/figure_2/month_$month/$key/$level"; bundle_path)
    end
end
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
observables = [location_1, location_2, location_3]
##
if !(use_bundle && has_ground_truth_bundle(bundle_path))
    historical_tas = common_array("historical", "tas"; ensemble_members = 45)
    historical_hurs = common_array("historical", "hurs"; ensemble_members = 29)
    ssp119_tas = common_array("ssp119", "tas"; ensemble_members = 45)
    ssp119_hurs = common_array("ssp119", "hurs"; ensemble_members = 29)
end
historical_temperatures = regression_variable("historical")
ssp119_temperatures = regression_variable("ssp119")

##
window = 2 
increment = 2 * window + 1
ts = 2015:2100
indmin = 6
min_temp = mean(ssp119_temperatures[indmin-window:indmin+window])
indmax = 81
max_temp = mean(ssp119_temperatures[indmax-window:indmax+window])
ts[indmin], ts[indmax]
##
month = 1
fig = Figure(; resolution)
observables = [location_1, location_3]
observable_names_base = [@sprintf("%.2fᵒ, %.2fᵒ ", longitude[tmp[1]], latitude[tmp[2]]) for tmp in [index_1, index_3]]
observable_names = observable_names_base .* ["(Jan)"]
for i in eachindex(observables)
    emulator.month[1] = month
    emulator.global_mean_temperature[1] = min_temp
    if i == 1
        common_options_1 = (; xlabel = "Temperature (K)", ylabel = "PDF")
    else
        common_options_1 = (; xlabel = "Temperature (K)")
    end
    ax1 = Axis(fig[1, i]; title = observable_names[i], common_options_1..., global_common_options...)
    if use_bundle && has_ground_truth_bundle(bundle_path)
        key = i == 1 ? "loc_157_46" : "loc_19_87"
        tas_hist = f2_samples(month, key, "low")
    else
        rfield = reshape(ssp119_tas[:, :, indmin-window:indmin+window, month, :], (192 * 96, 45 * increment))
        tas_hist = [observables[i](rfield[:, j]) for j in 1:45*increment]
    end
    hist!(ax1, tas_hist, bins = 20, color = (:royalblue, 0.2), normalization = :pdf, label = "2020 (SSP1-1.9)")
    emulator.month[1] = month
    emulator.global_mean_temperature[1] = max_temp

    if use_bundle && has_ground_truth_bundle(bundle_path)
        key = i == 1 ? "loc_157_46" : "loc_19_87"
        tas_hist = f2_samples(month, key, "high")
    else
        rfield = reshape(ssp119_tas[:, :, indmax-window:indmax+window, month, :], (192 * 96, 45*increment))
        tas_hist = [observables[i](rfield[:, j]) for j in 1:45*increment]
    end
    hist!(ax1, tas_hist, bins = 20, color = (:orangered2, 0.2), normalization = :pdf, label = "2095 (SSP1-1.9)")
    if i == 1
        axislegend(ax1, position = :rt, labelsize = legend_ls)
    end
end

observables = [global_mean_upper, global_mean_lower]
observable_names_base = ["Northern Hemisphere Mean ", "Southern Hemisphere Mean "]
observable_names = observable_names_base .* ["(Jan)"]
for i in eachindex(observables)
    emulator.month[1] = month
    emulator.global_mean_temperature[1] = min_temp
    if i == 1
        common_options_1 = (; xlabel = "Temperature (K)", ylabel = "PDF")
    else
        common_options_1 = (; xlabel = "Temperature (K)")
    end
    ax1 = Axis(fig[2, i]; title = observable_names[i], common_options_1..., global_common_options...)
    if use_bundle && has_ground_truth_bundle(bundle_path)
        key = i == 1 ? "upper" : "lower"
        tas_hist = f2_samples(month, key, "low")
    else
        rfield = reshape(ssp119_tas[:, :, indmin-window:indmin+window, month, :], (192 * 96, 45 * increment))
        tas_hist = [observables[i](rfield[:, j]) for j in 1:45*increment]
    end
    hist!(ax1, tas_hist, bins = 20, color = (:royalblue, 0.2), normalization = :pdf, label = "2020 (SSP1-1.9)")
    emulator.month[1] = month
    emulator.global_mean_temperature[1] = max_temp

    if use_bundle && has_ground_truth_bundle(bundle_path)
        key = i == 1 ? "upper" : "lower"
        tas_hist = f2_samples(month, key, "high")
    else
        rfield = reshape(ssp119_tas[:, :, indmax-window:indmax+window, month, :], (192 * 96, 45*increment))
        tas_hist = [observables[i](rfield[:, j]) for j in 1:45*increment]
    end
    hist!(ax1, tas_hist, bins = 20, color = (:orangered2, 0.2), normalization = :pdf, label = "2095 (SSP1-1.9)")
    if i == 1
        # axislegend(ax1, position = :rt, labelsize = legend_ls)
    end
end

month = 7
observables = [location_1, location_3]
observable_names_base = [@sprintf("%.2fᵒ, %.2fᵒ ", longitude[tmp[1]], latitude[tmp[2]]) for tmp in [index_1, index_3]]
observable_names = observable_names_base .* ["(Jul)"]
for i in eachindex(observables)
    emulator.month[1] = month
    emulator.global_mean_temperature[1] = min_temp
    if i == 1
        common_options_1 = (; xlabel = "Temperature (K)", ylabel = "PDF")
    else
        common_options_1 = (; xlabel = "Temperature (K)")
    end
    ax1 = Axis(fig[1, i+2]; title = observable_names[i], common_options_1..., global_common_options...)
    if use_bundle && has_ground_truth_bundle(bundle_path)
        key = i == 1 ? "loc_157_46" : "loc_19_87"
        tas_hist = f2_samples(month, key, "low")
    else
        rfield = reshape(ssp119_tas[:, :, indmin-window:indmin+window, month, :], (192 * 96, 45 * increment))
        tas_hist = [observables[i](rfield[:, j]) for j in 1:45*increment]
    end
    hist!(ax1, tas_hist, bins = 20, color = (:royalblue, 0.2), normalization = :pdf)
    emulator.month[1] = month
    emulator.global_mean_temperature[1] = max_temp

    if use_bundle && has_ground_truth_bundle(bundle_path)
        key = i == 1 ? "loc_157_46" : "loc_19_87"
        tas_hist = f2_samples(month, key, "high")
    else
        rfield = reshape(ssp119_tas[:, :, indmax-window:indmax+window, month, :], (192 * 96, 45*increment))
        tas_hist = [observables[i](rfield[:, j]) for j in 1:45*increment]
    end
    hist!(ax1, tas_hist, bins = 20, color = (:orangered2, 0.2), normalization = :pdf, label = "2050 (Data)")
end

observables = [global_mean_upper, global_mean_lower]
observable_names_base = ["Northern Hemisphere Mean ", "Southern Hemisphere Mean "]
observable_names = observable_names_base .* ["(Jul)"]
for i in eachindex(observables)
    emulator.month[1] = month
    emulator.global_mean_temperature[1] = min_temp
    if i == 1
        common_options_1 = (; xlabel = "Temperature (K)", ylabel = "PDF")
    else
        common_options_1 = (; xlabel = "Temperature (K)")
    end
    ax1 = Axis(fig[2, i+2]; title = observable_names[i], common_options_1..., global_common_options...)
    if use_bundle && has_ground_truth_bundle(bundle_path)
        key = i == 1 ? "upper" : "lower"
        tas_hist = f2_samples(month, key, "low")
    else
        rfield = reshape(ssp119_tas[:, :, indmin-window:indmin+window, month, :], (192 * 96, 45 * increment))
        tas_hist = [observables[i](rfield[:, j]) for j in 1:45*increment]
    end
    hist!(ax1, tas_hist, bins = 20, color = (:royalblue, 0.2), normalization = :pdf)
    emulator.month[1] = month
    emulator.global_mean_temperature[1] = max_temp

    if use_bundle && has_ground_truth_bundle(bundle_path)
        key = i == 1 ? "upper" : "lower"
        tas_hist = f2_samples(month, key, "high")
    else
        rfield = reshape(ssp119_tas[:, :, indmax-window:indmax+window, month, :], (192 * 96, 45*increment))
        tas_hist = [observables[i](rfield[:, j]) for j in 1:45*increment]
    end
    hist!(ax1, tas_hist, bins = 20, color = (:orangered2, 0.2), normalization = :pdf, label = "2050 (Data)")
    if i == 1
        # axislegend(ax2, position = :lt, labelsize = legend_ls)
    end
end
# display(fig)

save(figure_directory * "figure_2_parametric_assumption_point_global_statistics_tas_jan_july.png", fig)
@info "Figure 2 generated"
