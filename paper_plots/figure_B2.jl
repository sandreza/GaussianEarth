ts = 60
xls = 60 
yls = 60
tls = 60
legend_ls = 60
xlp = 20 
ylp = 20 
resolution = (1800, 1200) .* 2
common_options = (; titlesize = ts, xlabelsize = xls, ylabelsize = yls, xticklabelsize = tls, yticklabelsize = tls, xlabelpadding = xlp, ylabelpadding = ylp)
field_name = "tas"


# repeat for linear emulator
use_bundle = @isdefined(use_ground_truth_bundle) ? use_ground_truth_bundle : false
bundle_path = @isdefined(ground_truth_bundle_path) ? ground_truth_bundle_path : get_ground_truth_bundle_path()
scenarios = @isdefined(scenarios) ? scenarios : ["historical", "ssp585", "ssp119", "ssp245"]

if !@isdefined(temperatures)
    temperatures = [regression_variable(scenario) for scenario in scenarios]
end

if !@isdefined(latitude) || !@isdefined(longitude)
    hfile = h5open(save_directory * field_name * "_basis.hdf5", "r")
    latitude = read(hfile["latitude"])
    longitude = read(hfile["longitude"])
    close(hfile)
end

#read in linear emulator
#load in saved out tas emulator
hfile = h5open(save_directory * field_name * "_model.hdf5", "r")
μmodel = read(hfile["mean"])
Lmodel = read(hfile["L model"])
basis = read(hfile["basis"])
close(hfile)
linear_emulator = CovarEmulator(μmodel, Lmodel, basis)

linear_emulated_truth = zeros(Float32, 192, 96)
truth = zeros(Float32, 192, 96)
linear_emulated_truth_truth_error = zeros(Float32, 192, 96, length(scenarios))

tas_mean_monthly = nothing
if use_bundle && has_ground_truth_bundle(bundle_path)
    tas_mean_monthly = [
        read_ground_truth("tas/mean_monthly/$scenario"; bundle_path) for scenario in scenarios
    ]
elseif !@isdefined(tas_fields)
    error("tas_fields not defined and no ground-truth bundle found; set data_directory or provide bundle.")
end

for scenario_index in ProgressBar(eachindex(scenarios))
    for year in ProgressBar(eachindex(temperatures[scenario_index]))
        N = length(temperatures[scenario_index])
        for month in 1:12
            # tas 
            linear_emulator.global_mean_temperature[1] = temperatures[scenario_index][year]
            linear_emulator.month[1] = month
            linear_emulated_truth .= reshape(mean(linear_emulator), 192, 96)
            if tas_mean_monthly === nothing
                truth .= mean(tas_fields[scenario_index][:, :, year, month, :], dims=3)[:,:, 1]
            else
                truth .= tas_mean_monthly[scenario_index][:, :, year, month]
            end
            linear_emulated_truth_truth_error[:, :, scenario_index] .+= abs.(linear_emulated_truth .- truth) / ( N * 12)
        end
    end
end

fig = Figure(; resolution)
nlongitude = longitude .- 180
scenario_index = 1
cmap = Reverse(:linear_tritanopic_krjcw_5_98_c46_n256) #:plasma
ax = GeoAxis(fig[1,1]; title = "Historical", common_options...)
field =  linear_emulated_truth_truth_error[:, :, scenario_index]
shifted_field = circshift(field, (96, 0))
surface!(ax, nlongitude, latitude, shifted_field; colormap = cmap, colorrange = (0, 1), shading = false)
hidedecorations!(ax)
scenario_index = 2 
ax = GeoAxis(fig[1,2]; title = "SSP5-8.5", common_options...)
field =  linear_emulated_truth_truth_error[:, :, scenario_index]
shifted_field = circshift(field, (96, 0))
surface!(ax, nlongitude, latitude, shifted_field; colormap = cmap, colorrange = (0, 1), shading = false)
hidedecorations!(ax)
scenario_index = 3
ax = GeoAxis(fig[2,1]; title = "SSP1-1.9", common_options...)
field =  linear_emulated_truth_truth_error[:, :, scenario_index]
shifted_field = circshift(field, (96, 0))
surface!(ax, nlongitude, latitude, shifted_field; colormap = cmap, colorrange = (0, 1), shading = false)
hidedecorations!(ax)
scenario_index = 4
ax = GeoAxis(fig[2,2]; title = "SSP2-4.5", common_options...)
field =  linear_emulated_truth_truth_error[:, :, scenario_index]
shifted_field = circshift(field, (96, 0))
surface!(ax, nlongitude, latitude, shifted_field; colormap = cmap, colorrange = (0, 1), shading = false)
Colorbar(fig[1:2,3], colormap=cmap, colorrange=(0, 1), height = Relative(2/4), label = "Temperature Error (K)", labelsize = legend_ls, ticklabelsize = legend_ls)
hidedecorations!(ax)
save(figure_directory *  "figure_B2_emulated_truth_truth_error_unweighted_scenarios_tas_linear.png", fig)
# display(fig)
