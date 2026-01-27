ts = 60
xls = 60 
yls = 60
tls = 60
legend_ls = 60
xlp = 20 
ylp = 20 
resolution = (1800, 1200) .* 2
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
        historical_tas = common_array(scenario, "tas"; ensemble_members = 45)
        push!(tas_fields, historical_tas)
        historical_temperatures = regression_variable(scenario)
        push!(temperatures, historical_temperatures)
    end

    include("../emulator.jl")

    emulated_truth = zeros(Float32, 192, 96)
    truth = zeros(Float32, 192, 96)
    emulated_truth_truth_error = zeros(Float32, 192, 96, length(scenarios))

    for scenario_index in ProgressBar(eachindex(scenarios))
        for year in ProgressBar(eachindex(temperatures[scenario_index]))
            N = length(temperatures[scenario_index])
            for month in 1:12
                # tas 
                emulator.global_mean_temperature[1] = temperatures[scenario_index][year]
                emulator.month[1] = month
                emulated_truth .= reshape(mean(emulator), 192, 96)
                truth .= mean(tas_fields[scenario_index][:, :, year, month, :], dims=3)[:,:, 1]
                emulated_truth_truth_error[:, :, scenario_index] .+= abs.(emulated_truth .- truth) / ( N * 12)
            end
        end
    end
end


fig = Figure(; resolution)
nlongitude = longitude .- 180
scenario_index = 1
cmap = Reverse(:linear_tritanopic_krjcw_5_98_c46_n256) 
ax = GeoAxis(fig[1,1]; title = "Historical", common_options...)
field = emulated_truth_truth_error[:, :, scenario_index] 
shifted_field = circshift(field, (96, 0))
surface!(ax, nlongitude, latitude, shifted_field; colormap = cmap, colorrange = (0, 1), shading = NoShading)
hidedecorations!(ax)
scenario_index = 2 
ax = GeoAxis(fig[1,2]; title = "SSP5-8.5", common_options...)
field = emulated_truth_truth_error[:, :, scenario_index] 
shifted_field = circshift(field, (96, 0))
surface!(ax, nlongitude, latitude, shifted_field; colormap = cmap, colorrange = (0, 1), shading = NoShading)
hidedecorations!(ax)
scenario_index = 3
ax = GeoAxis(fig[2,1]; title = "SSP1-1.9", common_options...)
field = emulated_truth_truth_error[:, :, scenario_index] 
shifted_field = circshift(field, (96, 0))
surface!(ax, nlongitude, latitude, shifted_field; colormap = cmap, colorrange = (0, 1), shading = NoShading)
hidedecorations!(ax)
scenario_index = 4
ax = GeoAxis(fig[2,2]; title = "SSP2-4.5", common_options...)
field = emulated_truth_truth_error[:, :, scenario_index] 
shifted_field = circshift(field, (96, 0))
surface!(ax, nlongitude, latitude, shifted_field; colormap = cmap, colorrange = (0, 1), shading = NoShading)
Colorbar(fig[1:2,3], colormap=cmap, colorrange=(0, 1), height = Relative(2/4), label = "Temperature Error (K)", labelsize = legend_ls, ticklabelsize = legend_ls)
hidedecorations!(ax)
save(figure_directory * "figure_4_emulated_truth_truth_error_unweighted_scenarios_tas.png", fig) # this is the one currently being used in the paper
@info "Generated Figure 4."