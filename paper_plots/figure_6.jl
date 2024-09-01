using GeoMakie 

ts = 40
xls = 40 
yls = 40
tls = 40
legend_ls = 35
xlp = 15 
ylp = 15 
resolution = (1800, 600) .* 2
common_options = (; titlesize = ts, xlabelsize = xls, ylabelsize = yls, xticklabelsize = tls, yticklabelsize = tls, xlabelpadding = xlp, ylabelpadding = ylp)

if process_data
    save_directory = "/net/fs06/d3/sandre/GaussianEarthData/"
    data_directory = "/net/fs06/d3/mgeo/CMIP6/interim/"
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
    scenarios = ["historical"]
    for scenario in scenarios
        historical_tas = common_array(scenario, "tas"; ensemble_members = 45)
        # historical_hurs = common_array(scenario, "hurs"; ensemble_members = 29)
        historical_temperatures = regression_variable(scenario)
        push!(tas_fields, historical_tas)
        # push!(hurs_fields, historical_hurs)
        push!(temperatures, historical_temperatures)
    end
end

##
colormap = :balance
Random.seed!(124)
year_index = 1
month_index = 1
emulator.month[1] = month_index
emulator.global_mean_temperature[1] = temperatures[1][year_index]
emulator_mean = mean(emulator)
crange = (-12, 12)
fig = Figure(; resolution)
scenario_index = 1
ax = GeoAxis(fig[1,1]; title = "1850 January (Data Realization)", common_options...)
data_realization = tas_fields[1][:, :, year_index, month_index, 1] - mean(tas_fields[1][:, :, year_index, month_index, :] , dims = 3)[:, :, 1]
nlongitude = range(-180, 180, length = 192)
ndata_realization = circshift(data_realization, (96, 0))
surface!(ax, nlongitude, latitude, ndata_realization; colormap = colormap, colorrange = crange, shading = NoShading)
ax = GeoAxis(fig[1,2]; title = "1850 January (Data Realization)", common_options...)
data_realization = tas_fields[1][:, :, year_index, month_index, 2]  - mean(tas_fields[1][:, :, year_index, month_index, :] , dims = 3)[:, :, 1]
ndata_realization = circshift(data_realization, (96, 0))
surface!(ax, nlongitude, latitude, ndata_realization; colormap = colormap, colorrange = crange, shading = NoShading)

ax2 = GeoAxis(fig[1,3]; title = "1850 January (Emulator Realization)", common_options...)
emulator_realization = rand(emulator) - emulator_mean
r_emulator_realization = reshape(emulator_realization, (192, 96))
nemulator_realization = circshift(r_emulator_realization, (96, 0))
surface!(ax2, nlongitude, latitude, nemulator_realization; colormap = colormap, colorrange = crange, shading = NoShading)
ax2 = GeoAxis(fig[1,4]; title = "1850 January (Emulator Realization)", common_options...)
emulator_realization = rand(emulator) - emulator_mean
r_emulator_realization = reshape(emulator_realization, (192, 96))
nemulator_realization = circshift(r_emulator_realization, (96, 0))
surface!(ax2, nlongitude, latitude, nemulator_realization; colormap = colormap, colorrange = crange, shading = NoShading)
Colorbar(fig[1,5], label = "Temperature Anomaly (K)", colorrange = crange, colormap = colormap, height = Relative(1/3), labelsize = legend_ls, ticklabelsize = legend_ls)
save(figure_directory * "tas_realization_comparison.png", fig)