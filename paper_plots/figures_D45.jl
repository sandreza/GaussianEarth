ts = 60
xls = 60 
yls = 60
tls = 60
legend_ls = 60
xlp = 20 
ylp = 20 
resolution = (1800 * 2, 1200 *3) 
common_options = (; titlesize = ts, xlabelsize = xls, ylabelsize = yls, xticklabelsize = tls, yticklabelsize = tls, xlabelpadding = xlp, ylabelpadding = ylp)

field_name = "tas"
hfile = h5open(save_directory * field_name * "_basis.hdf5", "r")
latitude = read(hfile["latitude"])
longitude = read(hfile["longitude"])
metric = read(hfile["metric"])
close(hfile)
sqrt_f_metric = sqrt.(reshape(metric, 192 * 96))

include("../emulator.jl")

## construct raw fields
use_bundle = @isdefined(use_ground_truth_bundle) ? use_ground_truth_bundle : false
bundle_path = @isdefined(ground_truth_bundle_path) ? ground_truth_bundle_path : get_ground_truth_bundle_path()

tas_fields = []
std_fields = []
temperatures = []
scenarios = ["historical", "ssp585", "ssp119", "ssp245"]
for scenario in scenarios
    if use_bundle && has_ground_truth_bundle(bundle_path)
        mean_monthly = read_ground_truth("tas/mean_monthly/$scenario"; bundle_path)
        scenario_tas = reshape(mean_monthly, size(mean_monthly)..., 1)
    else
        scenario_tas = common_array(scenario, "tas"; ensemble_members = 45)
    end
    push!(tas_fields, scenario_tas)
    scenario_temperatures = regression_variable(scenario)
    push!(temperatures, scenario_temperatures)
    true_std = ensemble_averaging(scenario, field_name; return_std=true)
    push!(std_fields, true_std)
end

## construct emulated fields
emulated_truth = zeros(Float32, 192, 96, length(scenarios))
emulated_truth_std = zeros(Float32, 192, 96, length(scenarios))
ground_truth = zeros(Float32, 192, 96, length(scenarios))
ground_truth_std = zeros(Float32, 192, 96, length(scenarios))
for scenario_index in ProgressBar(eachindex(scenarios))
    emulator.global_mean_temperature[1] = temperatures[scenario_index][end] #for 2100
    for month in 1:12
        emulator.month[1] = month
        emulated_truth[:,:,scenario_index] .+= reshape(mean(emulator), 192, 96)
        ground_truth[:,:,scenario_index] .+= mean(tas_fields[scenario_index][:, :, end, month, :], dims=3)[:,:, 1]
        if month == 1
            emulated_truth_std[:,:,scenario_index] .= reshape(sqrt.(variance(emulator)), 192, 96)
            ground_truth_std[:,:,scenario_index] .= std_fields[scenario_index][1][:, :, end, month]
        end
    end
    emulated_truth[:,:,scenario_index] ./= 12
    ground_truth[:,:,scenario_index] ./= 12
end

for plot_std in (false, true)

    fig = Figure(; resolution)
    nlongitude = longitude .- 180
    cmap = plot_std ? :viridis : :plasma
    crange = plot_std ? extrema(ground_truth_std[:,:,2]) : extrema(ground_truth[:,:,2])
    scenario_index = 2 
    ax = GeoAxis(fig[1,1]; title = "SSP5-8.5 (Emulator)", common_options...)
    field = plot_std ? emulated_truth_std[:,:, scenario_index] : emulated_truth[:, :, scenario_index]
    shifted_field = circshift(field, (96, 0))
    surface!(ax, nlongitude, latitude, shifted_field; colormap = cmap, colorrange = crange, shading = false) 
    hidedecorations!(ax)
    ax = GeoAxis(fig[1,2]; title = "SSP5-8.5 (MPI)", common_options...)
    field = plot_std ? ground_truth_std[:,:, scenario_index] : ground_truth[:, :, scenario_index]
    shifted_field = circshift(field, (96, 0))
    surface!(ax, nlongitude, latitude, shifted_field; colormap = cmap, colorrange = crange, shading = false) 
    hidedecorations!(ax)
    scenario_index = 3
    ax = GeoAxis(fig[3,1]; title = "SSP1-1.9 (Emulator)", common_options...)
    field = plot_std ? emulated_truth_std[:,:, scenario_index] : emulated_truth[:, :, scenario_index]
    shifted_field = circshift(field, (96, 0))
    surface!(ax, nlongitude, latitude, shifted_field; colormap = cmap, colorrange = crange, shading = false)
    hidedecorations!(ax)
    ax = GeoAxis(fig[3,2]; title = "SSP1-1.9 (MPI)", common_options...)
    field = plot_std ? ground_truth_std[:,:, scenario_index] : ground_truth[:, :, scenario_index]
    shifted_field = circshift(field, (96, 0))
    surface!(ax, nlongitude, latitude, shifted_field; colormap = cmap, colorrange = crange, shading = false) 
    hidedecorations!(ax)
    scenario_index = 4
    ax = GeoAxis(fig[2,1]; title = "SSP2-4.5 (Emulator)", common_options...)
    field = plot_std ? emulated_truth_std[:,:, scenario_index] : emulated_truth[:, :, scenario_index]
    shifted_field = circshift(field, (96, 0))
    surface!(ax, nlongitude, latitude, shifted_field; colormap = cmap, colorrange = crange, shading = false)
    hidedecorations!(ax)
    ax = GeoAxis(fig[2,2]; title = "SSP2-4.5 (MPI)", common_options...)
    field = plot_std ? ground_truth_std[:,:, scenario_index] : ground_truth[:, :, scenario_index]
    shifted_field = circshift(field, (96, 0))
    surface!(ax, nlongitude, latitude, shifted_field; colormap = cmap, colorrange = crange, shading = false)
    hidedecorations!(ax)
    cb_label = plot_std ? "Standard Deviation (K)" : "Temperature (K)"
    Colorbar(fig[4,1], colormap=cmap, colorrange=crange, width = Relative(3/4), height=Relative(1/30), vertical=false, label = cb_label, labelsize = legend_ls, ticklabelsize = legend_ls)
    # Colorbar(fig[4,2], colormap=cmap_std, colorrange=crange_std, width = Relative(3/4), height=Relative(1/30), vertical=false, label = "Standard Deviation (K)", labelsize = legend_ls, ticklabelsize = legend_ls)
    # display(fig)
    if plot_std
        save(figure_directory * "figure_D5_emulated_mpi_comparison_2100_std.png", fig)
    else
        save(figure_directory * "figure_D4_emulated_mpi_comparison_2100_mean.png", fig)
    end
end
