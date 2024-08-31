using GeoMakie 

ts = 40
xls = 40 
yls = 40
tls = 40
legend_ls = 35
resolution = (1800, 600)
common_options = (; titlesize = ts, xlabelsize = xls, ylabelsize = yls, xticklabelsize = tls, yticklabelsize = tls)

if process_data 
    field_name = "tas"
    hfile = h5open(save_directory * field_name * "_basis.hdf5", "r")
    latitude = read(hfile["latitude"])
    longitude = read(hfile["longitude"])
    metric = read(hfile["metric"])
    close(hfile)
    sqrt_f_metric = sqrt.(reshape(metric, 192 * 96))

    include("emulator.jl")

    tas_fields = []
    hurs_fields = []
    temperatures = []
    scenarios = ["historical", "ssp585", "ssp119", "ssp245"]
    for scenario in scenarios
        historical_tas = common_array(scenario, "tas"; ensemble_members = 45)
        historical_temperatures = regression_variable(scenario)
        push!(tas_fields, historical_tas)
        push!(temperatures, historical_temperatures)
    end

    emulated_truth = zeros(Float32, 192, 96)
    truth = zeros(Float32, 192, 96)
    emulated_truth_truth_error = zeros(Float32, 192, 96, length(scenarios))
    emulated_truth_truth_error_weighted = zeros(Float32, 192, 96, length(scenarios))

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
                emulated_truth_truth_error_weighted[:, :, scenario_index] .+= sqrt.(((abs.(emulated_truth .- truth) ).^2) .* metric) / ( N * 12)
            end
        end
    end

end

fig = Figure(; resolution)
nlongitude = longitude .- 180
scenario_index = 1
ax = GeoAxis(fig[1,1]; title = "Historical", common_options...)
field = emulated_truth_truth_error[:, :, scenario_index]
shifted_field = circshift(field, (96, 0))
surface!(ax, nlongitude, latitude, shifted_field; colormap = :afmhot, colorrange = (0, 1), shading = NoShading)

scenario_index = 2 
ax = GeoAxis(fig[1,2]; title = "SSP5-8.5", common_options...)
field = emulated_truth_truth_error[:, :, scenario_index]
shifted_field = circshift(field, (96, 0))
surface!(ax, nlongitude, latitude, shifted_field; colormap = :afmhot, colorrange = (0, 1), shading = NoShading)
scenario_index = 3
ax = GeoAxis(fig[1,3]; title = "SSP1-1.9", common_options...)
field = emulated_truth_truth_error[:, :, scenario_index]
shifted_field = circshift(field, (96, 0))
surface!(ax, nlongitude, latitude, shifted_field; colormap = :afmhot, colorrange = (0, 1), shading = NoShading)
Colorbar(fig[1,4], colormap=:afmhot, colorrange=(0, 1), height = Relative(2/4), label = "Temperature Error (K)")
save(figure_directory * "emulated_truth_truth_error_unweighted_scenarios_tas.png", fig)