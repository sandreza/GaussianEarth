using CairoMakie
using GeoMakie

field_name = "tas"
hfile = h5open(save_directory * field_name * "_basis.hdf5", "r")
latitude = read(hfile["latitude"])
longitude = read(hfile["longitude"])
metric = read(hfile["metric"])
close(hfile)
sqrt_f_metric = sqrt.(reshape(metric, 192 * 96))

include("emulator.jl")
include("emulator_hurs.jl")

tas_fields = []
hurs_fields = []
temperatures = []
scenarios = ["historical", "ssp585", "ssp119", "ssp245"]
for scenario in scenarios
    historical_tas = common_array(scenario, "tas"; ensemble_members = 45)
    historical_hurs = common_array(scenario, "hurs"; ensemble_members = 29)
    historical_temperatures = regression_variable(scenario)
    push!(tas_fields, historical_tas)
    push!(hurs_fields, historical_hurs)
    push!(temperatures, historical_temperatures)
end


##
## Plotting

model = zeros(Float32, 192, 96)
truth = zeros(Float32, 192, 96)
model_truth_error = zeros(Float32, 192, 96, length(scenarios))
model_truth_error_weighted = zeros(Float32, 192, 96, length(scenarios))
model_truth_error_hurs = zeros(Float32, 192, 96, length(scenarios))
model_truth_error_weighted_hurs = zeros(Float32, 192, 96, length(scenarios))

scenario_index = 1
year = 1
for scenario_index in ProgressBar(eachindex(scenarios))
    for year in ProgressBar(eachindex(temperatures[scenario_index]))
        N = length(temperatures[scenario_index])
        for month in 1:12
            # tas 
            emulator.global_mean_temperature[1] = temperatures[scenario_index][year]
            emulator.month[1] = month
            model .= reshape(mean(emulator), 192, 96)
            truth .= mean(tas_fields[scenario_index][:, :, year, month, :], dims=3)[:,:, 1]
            model_truth_error[:, :, scenario_index] .+= abs.(model .- truth) / ( N * 12)
            model_truth_error_weighted[:, :, scenario_index] .+= sqrt.(((abs.(model .- truth) ).^2) .* metric) / ( N * 12)
            # hurs 
            emulator_hurs.global_mean_temperature[1] = temperatures[scenario_index][year]
            emulator_hurs.month[1] = month
            model .= reshape(mean(emulator_hurs), 192, 96)
            truth .= mean(hurs_fields[scenario_index][:, :, year, month, :], dims=3)[:,:, 1]
            model_truth_error_hurs[:, :, scenario_index] .+= abs.(model .- truth) / ( N * 12)
            model_truth_error_weighted_hurs[:, :, scenario_index] .+= sqrt.(((abs.(model .- truth) ).^2) .* metric) / ( N * 12)
        end
    end
end

##
scenario_index = 1
fig = Figure(resolution = (1000, 1000))
ax = GeoAxis(fig[1,1])
surface!(ax, longitude, latitude, model_truth_error_weighted[:, :, scenario_index]; colormap = :plasma, shading = NoShading)
save("model_truth_error_weighted.png", fig)
##
nlongitude = range(-180, 180, length = 192)
scenario_index = 1
fig = Figure(resolution = (2000, 1000))

ax = GeoAxis(fig[1,1]; title = "Historical")

field = model_truth_error[:, :, scenario_index]
shifted_field = circshift(field, (96, 0))
surface!(ax, nlongitude, latitude, shifted_field; colormap = :afmhot, colorrange = (0, 1), shading = NoShading)
Colorbar(fig[1,2], colormap=:afmhot, colorrange=(0, 1), height = Relative(2/4), label = "Temperature Error (K)")
ax2 = GeoAxis(fig[1,3]; title = "Historical")

field = model_truth_error_hurs[:, :, scenario_index]
shifted_field = circshift(field, (96, 0))
surface!(ax2, nlongitude, latitude, shifted_field; colormap = :plasma, colorrange = (0, 2), shading = NoShading)
Colorbar(fig[1,4], colormap=:plasma, colorrange=(0, 2), height = Relative(2/4), label = "Relative Humidity Error (%)")
save("model_truth_error_unweighted.png", fig)
##


fig = Figure(resolution = (3000, 2000))
scenario_index = 1
ax = GeoAxis(fig[1,1]; title = "Historical")
field = model_truth_error[:, :, scenario_index]
shifted_field = circshift(field, (96, 0))
surface!(ax, nlongitude, latitude, shifted_field; colormap = :afmhot, colorrange = (0, 1), shading = NoShading)
ax2 = GeoAxis(fig[2,1]; title = "Historical")
field = model_truth_error_hurs[:, :, scenario_index]
shifted_field = circshift(field, (96, 0))
surface!(ax2, nlongitude, latitude, shifted_field; colormap = :plasma, colorrange = (0, 2), shading = NoShading)

scenario_index = 2 
ax = GeoAxis(fig[1,2]; title = "SSP585")
field = model_truth_error[:, :, scenario_index]
shifted_field = circshift(field, (96, 0))
surface!(ax, nlongitude, latitude, shifted_field; colormap = :afmhot, colorrange = (0, 1), shading = NoShading)

ax2 = GeoAxis(fig[2,2]; title = "SSP585")
field = model_truth_error_hurs[:, :, scenario_index]
shifted_field = circshift(field, (96, 0))
surface!(ax2, nlongitude, latitude, shifted_field; colormap = :plasma, colorrange = (0, 2), shading = NoShading)

scenario_index = 3
ax = GeoAxis(fig[1,3]; title = "SSP119")
field = model_truth_error[:, :, scenario_index]
shifted_field = circshift(field, (96, 0))
surface!(ax, nlongitude, latitude, shifted_field; colormap = :afmhot, colorrange = (0, 1), shading = NoShading)
Colorbar(fig[1,4], colormap=:afmhot, colorrange=(0, 1), height = Relative(2/4), label = "Temperature Error (K)")

ax2 = GeoAxis(fig[2,3]; title = "SSP119")
field = model_truth_error_hurs[:, :, scenario_index]
shifted_field = circshift(field, (96, 0))
surface!(ax2, nlongitude, latitude, shifted_field; colormap = :plasma, colorrange = (0, 2), shading = NoShading)
Colorbar(fig[2,4], colormap=:plasma, colorrange=(0, 2), height = Relative(2/4), label = "Relative Humidity Error (%)")
save("model_truth_error_unweighted_scenarios.png", fig)

##
fig = Figure(resolution = (3000, 1000))
scenario_index = 1
ax = GeoAxis(fig[1,1]; title = "Historical")
field = model_truth_error[:, :, scenario_index]
shifted_field = circshift(field, (96, 0))
surface!(ax, nlongitude, latitude, shifted_field; colormap = :afmhot, colorrange = (0, 1), shading = NoShading)

scenario_index = 2 
ax = GeoAxis(fig[1,2]; title = "SSP585")
field = model_truth_error[:, :, scenario_index]
shifted_field = circshift(field, (96, 0))
surface!(ax, nlongitude, latitude, shifted_field; colormap = :afmhot, colorrange = (0, 1), shading = NoShading)
scenario_index = 3
ax = GeoAxis(fig[1,3]; title = "SSP119")
field = model_truth_error[:, :, scenario_index]
shifted_field = circshift(field, (96, 0))
surface!(ax, nlongitude, latitude, shifted_field; colormap = :afmhot, colorrange = (0, 1), shading = NoShading)
Colorbar(fig[1,4], colormap=:afmhot, colorrange=(0, 1), height = Relative(2/4), label = "Temperature Error (K)")
save("model_truth_error_unweighted_scenarios_tas.png", fig)
##
fig = Figure(resolution = (3000, 1000))
scenario_index = 1
ax = GeoAxis(fig[1,1]; title = "Historical")
field = model_truth_error_hurs[:, :, scenario_index]
shifted_field = circshift(field, (96, 0))
surface!(ax, nlongitude, latitude, shifted_field; colormap = :plasma, colorrange = (0, 2), shading = NoShading)
scenario_index = 2
ax = GeoAxis(fig[1,2]; title = "SSP585")
field = model_truth_error_hurs[:, :, scenario_index]
shifted_field = circshift(field, (96, 0))
surface!(ax, nlongitude, latitude, shifted_field; colormap = :plasma, colorrange = (0, 2), shading = NoShading)
scenario_index = 3
ax = GeoAxis(fig[1,3]; title = "SSP119")
field = model_truth_error_hurs[:, :, scenario_index]
shifted_field = circshift(field, (96, 0))
surface!(ax, nlongitude, latitude, shifted_field; colormap = :plasma, colorrange = (0, 2), shading = NoShading)
Colorbar(fig[1,4], colormap=:plasma, colorrange=(0, 2), height = Relative(2/4), label = "Relative Humidity Error (%)")
save("model_truth_error_unweighted_scenarios_hurs.png", fig)
