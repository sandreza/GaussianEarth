using GeoMakie 

ts = 60
xls = 60 
yls = 60
tls = 60
legend_ls = 60
xlp = 20 
ylp = 20 
resolution = (1800 * 2, 1200 *3) 
common_options = (; titlesize = ts, xlabelsize = xls, ylabelsize = yls, xticklabelsize = tls, yticklabelsize = tls, xlabelpadding = xlp, ylabelpadding = ylp)

if process_data 
    field_name = "tas"
    hfile = h5open(save_directory * field_name * "_basis.hdf5", "r")
    latitude = read(hfile["latitude"])
    longitude = read(hfile["longitude"])
    metric = read(hfile["metric"])
    close(hfile)
    sqrt_f_metric = sqrt.(reshape(metric, 192 * 96))

    tas_fields = []
    temperatures = []
    scenarios = ["historical", "ssp585", "ssp119", "ssp245"]
    for scenario in scenarios
        historical_temperatures = regression_variable(scenario)
        push!(temperatures, historical_temperatures)
    end

    emulated_truth = zeros(Float32, 192, 96, length(scenarios))
    emulated_truth_std = zeros(Float32, 192, 96, length(scenarios))
    for scenario_index in ProgressBar(eachindex(scenarios))
        emulator.global_mean_temperature[1] = temperatures[scenario_index][end]
        for month in 1:12
            emulator.month[1] = month
            emulated_truth[:,:,scenario_index] .+= reshape(mean(emulator), 192, 96)
            if month == 1
                emulated_truth_std[:,:,scenario_index] .+= reshape(sqrt.(variance(emulator)), 192, 96)
            end
        end
        emulated_truth[:,:,scenario_index] ./= 12
    end
end

fig = Figure(; resolution)
nlongitude = longitude .- 180
cmap = :plasma
cmap_std = :viridis
crange = extrema(emulated_truth[:,:,2])
crange_std = extrema(emulated_truth_std[:,:,2])
scenario_index = 2 
ax = GeoAxis(fig[1,1]; title = "SSP5-8.5", common_options...)
field = emulated_truth[:, :, scenario_index]
shifted_field = circshift(field, (96, 0))
surface!(ax, nlongitude, latitude, shifted_field; colormap = cmap, colorrange = crange, shading = NoShading) 
hidedecorations!(ax)
ax = GeoAxis(fig[1,2]; title = "SSP5-8.5", common_options...)
field = emulated_truth_std[:, :, scenario_index]
shifted_field = circshift(field, (96, 0))
surface!(ax, nlongitude, latitude, shifted_field; colormap = cmap_std, colorrange = crange_std, shading = NoShading) 
hidedecorations!(ax)
scenario_index = 3
ax = GeoAxis(fig[3,1]; title = "SSP1-1.9", common_options...)
field = emulated_truth[:, :, scenario_index]
shifted_field = circshift(field, (96, 0))
surface!(ax, nlongitude, latitude, shifted_field; colormap = cmap, colorrange = crange, shading = NoShading)
hidedecorations!(ax)
ax = GeoAxis(fig[3,2]; title = "SSP1-1.9", common_options...)
field = emulated_truth_std[:, :, scenario_index]
shifted_field = circshift(field, (96, 0))
surface!(ax, nlongitude, latitude, shifted_field; colormap = cmap_std, colorrange = crange_std, shading = NoShading) 
hidedecorations!(ax)
scenario_index = 4
ax = GeoAxis(fig[2,1]; title = "SSP2-4.5", common_options...)
field = emulated_truth[:, :, scenario_index]
shifted_field = circshift(field, (96, 0))
surface!(ax, nlongitude, latitude, shifted_field; colormap = cmap, colorrange = crange, shading = NoShading)
hidedecorations!(ax)
ax = GeoAxis(fig[2,2]; title = "SSP2-4.5", common_options...)
field = emulated_truth_std[:, :, scenario_index]
shifted_field = circshift(field, (96, 0))
surface!(ax, nlongitude, latitude, shifted_field; colormap = cmap_std, colorrange = crange_std, shading = NoShading)
hidedecorations!(ax)
Colorbar(fig[4,1], colormap=cmap, colorrange=crange, width = Relative(3/4), height=Relative(1/30), vertical=false, label = "Temperature (K)", labelsize = legend_ls, ticklabelsize = legend_ls)
Colorbar(fig[4,2], colormap=cmap_std, colorrange=crange_std, width = Relative(3/4), height=Relative(1/30), vertical=false, label = "Standard Deviation (K)", labelsize = legend_ls, ticklabelsize = legend_ls)
# display(fig)
save(figure_directory * "emulated_truth_unweighted_scenarios_tas.png", fig)