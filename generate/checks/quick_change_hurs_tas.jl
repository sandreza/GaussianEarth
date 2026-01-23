indmin = 100
min_temp = historical_temperatures[indmin]
indmax = 43
max_temp = ssp585_temperatures[indmax]
##
window = 2 
increment = 2 * window + 1
month = 1
fig = Figure(resolution = (3000, 600))
observable_names = ["Global Mean (Jan)", "Global Upper Hemisphere Mean (Jan)", "Global Lower Hemisphere Mean (Jan)"]
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
    ax1 = Axis(fig[1, i]; title = observable_names[i], common_options_1...)
    rfield = reshape(historical_tas[:, :, indmin-window:indmin+window, month, :], (192 * 96, 45 * increment))
    tas_hist = [observables[i](rfield[:, j]) for j in 1:45*increment]
    xs, ys = gaussian_grid(tasmean[i], tasstd[i])
    hist!(ax1, tas_hist, bins = 20, color = (:royalblue, 0.2), normalization = :pdf)
    lines!(ax1, xs, ys, color = :blue)
    ax2 = Axis(fig[2, i]; title = observable_names[i], common_options_2...)
    rfield = reshape(historical_hurs[:, :, indmin-window:indmin+window, month, :], (192 * 96, 29*increment))
    hurs_hist = [observables[i](rfield[:, j]) for j in 1:29*increment]
    hist!(ax2, hurs_hist, bins = 20, color = (:royalblue, 0.2), normalization = :pdf, label = "1950 (Data)")
    xs, ys = gaussian_grid(hursmean[i], hursstd[i])
    lines!(ax2, xs, ys, color = :blue, label = "1950 (Emulator)")

    emulator.month[1] = month
    emulator.global_mean_temperature[1] = max_temp
    emulator_hurs.month[1] = month
    emulator_hurs.global_mean_temperature[1] = max_temp
    tasmean, tasstd = emulator_mean_variance_linear_functionals(observables, emulator)
    hursmean, hursstd = emulator_mean_variance_linear_functionals(observables, emulator_hurs)

    xs, ys = gaussian_grid(tasmean[i], tasstd[i])
    lines!(ax1, xs, ys, color = :orange)
    xs, ys = gaussian_grid(hursmean[i], hursstd[i])
    lines!(ax2, xs, ys, color = :orange, label = "2050 (Emulator)")

    rfield = reshape(ssp585_tas[:, :, indmax-window:indmax+window, month, :], (192 * 96, 45*increment))
    tas_hist = [observables[i](rfield[:, j]) for j in 1:45*increment]
    rfield = reshape(ssp585_hurs[:, :, indmax-window:indmax+window, month, :], (192 * 96, 29*increment))
    hurs_hist = [observables[i](rfield[:, j]) for j in 1:29*increment]
    hist!(ax1, tas_hist, bins = 20, color = (:orangered2, 0.2), normalization = :pdf, label = "2050 (Data)")
    hist!(ax2, hurs_hist, bins = 20, color = (:orangered2, 0.2), normalization = :pdf, label = "2050 (Data)")
    if i == 1
        axislegend(ax2, position = :lt)
    end
end

display(fig)

save("climage_change_shifts_hurs_tas_with_data.png", fig)


month = 7


observable_names = ["Global Mean (Jul)", "Global Upper Hemisphere Mean (Jul)", "Global Lower Hemisphere Mean (Jul)"]
for i in eachindex(observables)
    emulator.month[1] = month
    emulator.global_mean_temperature[1] = min_temp
    emulator_hurs.month[1] = month
    emulator_hurs.global_mean_temperature[1] = min_temp
    tasmean, tasstd = emulator_mean_variance_linear_functionals(observables, emulator)
    hursmean, hursstd = emulator_mean_variance_linear_functionals(observables, emulator_hurs)

    common_options_1 = (; xlabel = "Temperature (K)")
    common_options_2 = (; xlabel = "Relative Humidity (%)")

    ax1 = Axis(fig[1, i+3]; title = observable_names[i], common_options_1...)
    rfield = reshape(historical_tas[:, :, indmin-window:indmin+window, month, :], (192 * 96, 45 * increment))
    tas_hist = [observables[i](rfield[:, j]) for j in 1:45*increment]
    xs, ys = gaussian_grid(tasmean[i], tasstd[i])
    hist!(ax1, tas_hist, bins = 20, color = (:royalblue, 0.2), normalization = :pdf)
    lines!(ax1, xs, ys, color = :blue)
    ax2 = Axis(fig[2, i+3]; title = observable_names[i], common_options_2...)
    rfield = reshape(historical_hurs[:, :, indmin-window:indmin+window, month, :], (192 * 96, 29*increment))
    hurs_hist = [observables[i](rfield[:, j]) for j in 1:29*increment]
    hist!(ax2, hurs_hist, bins = 20, color = (:royalblue, 0.2), normalization = :pdf)
    xs, ys = gaussian_grid(hursmean[i], hursstd[i])
    lines!(ax2, xs, ys, color = :blue)

    emulator.month[1] = month
    emulator.global_mean_temperature[1] = max_temp
    emulator_hurs.month[1] = month
    emulator_hurs.global_mean_temperature[1] = max_temp
    tasmean, tasstd = emulator_mean_variance_linear_functionals(observables, emulator)
    hursmean, hursstd = emulator_mean_variance_linear_functionals(observables, emulator_hurs)

    xs, ys = gaussian_grid(tasmean[i], tasstd[i])
    lines!(ax1, xs, ys, color = :orange)
    xs, ys = gaussian_grid(hursmean[i], hursstd[i])
    lines!(ax2, xs, ys, color = :orange)

    rfield = reshape(ssp585_tas[:, :, indmax-window:indmax+window, month, :], (192 * 96, 45*increment))
    tas_hist = [observables[i](rfield[:, j]) for j in 1:45*increment]
    rfield = reshape(ssp585_hurs[:, :, indmax-window:indmax+window, month, :], (192 * 96, 29*increment))
    hurs_hist = [observables[i](rfield[:, j]) for j in 1:29*increment]
    hist!(ax1, tas_hist, bins = 20, color = (:orangered2, 0.2), normalization = :pdf)
    hist!(ax2, hurs_hist, bins = 20, color = (:orangered2, 0.2), normalization = :pdf)
end

display(fig)

save("climage_change_shifts_hurs_tas_with_data_jan_july.png", fig)