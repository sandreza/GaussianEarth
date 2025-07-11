#figure 5 but its a difference map

using GeoMakie 

ts = 60
xls = 60 
yls = 60
tls = 60
legend_ls = 60
xlp = 20 
ylp = 20 
resolution = (1800, 1200) .* 2
common_options = (; titlesize = ts, xlabelsize = xls, ylabelsize = yls, xticklabelsize = tls, yticklabelsize = tls, xlabelpadding = xlp, ylabelpadding = ylp)

# calculatem errors for quadratic emulator
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
        historical_temperatures = regression_variable(scenario)
        push!(tas_fields, historical_tas)
        push!(temperatures, historical_temperatures)
    end

    include("emulator.jl")

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
                # crucially the errors are calculated per month and then averaged
            end
        end
    end 
end

# repeat for INSTEAD the original mean trend difference
if process_data 

    truth = zeros(Float32, 192, 96, 12, 2, length(scenarios))
    truth_diff = zeros(Float32, 192, 96, 12, length(scenarios))

    for scenario_index in ProgressBar(eachindex(scenarios))
        N = length(temperatures[scenario_index])
        # for year in [1, N]

            truth[:,:,:,1, scenario_index] = mean(tas_fields[scenario_index][:, :, 1, :, :], dims=4)[:,:, :,1] # also even taking the ensemble average here vs later is a choice
            truth[:,:,:,2, scenario_index] = mean(tas_fields[scenario_index][:, :, N, :, :], dims=4)[:,:,:, 1] # this is relying on the assumption that the difference is maximal for the last year FIX !!
            truth_diff .= truth[:,:,:,2,:] .- truth[:,:,:,1,:]
            # truth .= mean(tas_fields[scenario_index][:, :, year, month, :], dims=3)[:,:, 1]
            # linear_emulated_truth_truth_error[:, :, scenario_index] .+= abs.(linear_emulated_truth .- truth) / ( N * 12) 
    end
end

avg_diff = mean(truth_diff, dims=3)[:,:,1, :] # average over the 12 months
# calculate the difference
relative =  emulated_truth_truth_error ./ avg_diff

plot_linear = false
if plot_linear
    fig = Figure(; resolution)
    nlongitude = longitude .- 180
    scenario_index = 1
    cmap = Reverse(:linear_tritanopic_krjcw_5_98_c46_n256) #:plasma
    ax = GeoAxis(fig[1,1]; title = "Historical", common_options...)
    field =  linear_emulated_truth_truth_error[:, :, scenario_index]
    shifted_field = circshift(field, (96, 0))
    surface!(ax, nlongitude, latitude, shifted_field; colormap = cmap, colorrange = (0, 1), shading = NoShading)
    hidedecorations!(ax)
    scenario_index = 2 
    ax = GeoAxis(fig[1,2]; title = "SSP5-8.5", common_options...)
    field =  linear_emulated_truth_truth_error[:, :, scenario_index]
    shifted_field = circshift(field, (96, 0))
    surface!(ax, nlongitude, latitude, shifted_field; colormap = cmap, colorrange = (0, 1), shading = NoShading)
    hidedecorations!(ax)
    scenario_index = 3
    ax = GeoAxis(fig[2,1]; title = "SSP1-1.9", common_options...)
    field =  linear_emulated_truth_truth_error[:, :, scenario_index]
    shifted_field = circshift(field, (96, 0))
    surface!(ax, nlongitude, latitude, shifted_field; colormap = cmap, colorrange = (0, 1), shading = NoShading)
    hidedecorations!(ax)
    scenario_index = 4
    ax = GeoAxis(fig[2,2]; title = "SSP2-4.5", common_options...)
    field =  linear_emulated_truth_truth_error[:, :, scenario_index]
    shifted_field = circshift(field, (96, 0))
    surface!(ax, nlongitude, latitude, shifted_field; colormap = cmap, colorrange = (0, 1), shading = NoShading)
    Colorbar(fig[1:2,3], colormap=cmap, colorrange=(0, 1), height = Relative(2/4), label = "Temperature Error (K)", labelsize = legend_ls, ticklabelsize = legend_ls)
    hidedecorations!(ax)
    # save(figure_directory *  "emulated_truth_truth_error_unweighted_scenarios_tas_linear.png", fig)
    display(fig)
else
    fig = Figure(; resolution)
    nlongitude = longitude .- 180
    scenario_index = 1
    # cmap = Reverse(:linear_tritanopic_krjcw_5_98_c46_n256) #:plasma
    cmap = :balance
    ax = GeoAxis(fig[1,1]; title = "Historical", common_options...)
    field = relative[:, :, scenario_index]
    shifted_field = circshift(field, (96, 0))
    surface!(ax, nlongitude, latitude, shifted_field; colormap = cmap, colorrange = (-1, 1), shading = NoShading)
    hidedecorations!(ax)
    scenario_index = 2 
    ax = GeoAxis(fig[1,2]; title = "SSP5-8.5", common_options...)
    field = relative[:, :, scenario_index]
    shifted_field = circshift(field, (96, 0))
    surface!(ax, nlongitude, latitude, shifted_field; colormap = cmap, colorrange = (-1, 1), shading = NoShading)
    hidedecorations!(ax)
    scenario_index = 3
    ax = GeoAxis(fig[2,1]; title = "SSP1-1.9", common_options...)
    field = relative[:, :, scenario_index]
    shifted_field = circshift(field, (96, 0))
    surface!(ax, nlongitude, latitude, shifted_field; colormap = cmap, colorrange = (-1, 1), shading = NoShading)
    hidedecorations!(ax)
    scenario_index = 4
    ax = GeoAxis(fig[2,2]; title = "SSP2-4.5", common_options...)
    field = relative[:, :, scenario_index]
    shifted_field = circshift(field, (96, 0))
    surface!(ax, nlongitude, latitude, shifted_field; colormap = cmap, colorrange = (-1, 1), shading = NoShading)
    Colorbar(fig[1:2,3], colormap=cmap, colorrange=(-1, 1), height = Relative(2/4), label = "Difference in Temperature Error (K)", labelsize = legend_ls, ticklabelsize = legend_ls)
    hidedecorations!(ax)
    # save(figure_directory * "emulated_truth_difference_errors_tas.png", fig)
    display(fig)
end