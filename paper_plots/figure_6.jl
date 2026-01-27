ts = 40
xls = 40 
yls = 40
tls = 40
legend_ls = 35
resolution = (300, 140) .* 9
common_options = (; titlesize = ts, xlabelsize = xls, ylabelsize = yls, xticklabelsize = tls, yticklabelsize = tls)
ts = 40
xls = 40 
yls = 40
tls = 40
legend_ls = 35
common_options_2 = (; titlesize = ts, xlabelsize = xls, ylabelsize = yls, xticklabelsize = tls, yticklabelsize = tls)
if process_data
    use_bundle = @isdefined(use_ground_truth_bundle) ? use_ground_truth_bundle : false
    bundle_path = @isdefined(ground_truth_bundle_path) ? ground_truth_bundle_path : get_ground_truth_bundle_path()

    month = 1
    field = "tas"

    hfile = h5open(save_directory * field * "_basis.hdf5", "r")
    latitude = read(hfile["latitude"])
    longitude = read(hfile["longitude"])
    metric = read(hfile["metric"])
    close(hfile)

    if use_bundle && has_ground_truth_bundle(bundle_path)
        latitude_mean = read_ground_truth("stats/figure_6/latitude_mean"; bundle_path)
        latitude_std = read_ground_truth("stats/figure_6/latitude_std"; bundle_path)
        lat_samples = Dict(i => read_ground_truth("samples/figure_6/lat_$i"; bundle_path) for i in [1, 24, 48, 72, 96])
    else
        eof_mode, temperature = concatenate_regression(field, ["historical"])
        eofs = eof_mode[:,:, 1:45] # [1:48..., 50:50...]]

        historical_field = common_array("historical", field)
        hfile = h5open(save_directory * field * "_basis.hdf5", "r")
        latitude = read(hfile["latitude"])
        longitude = read(hfile["longitude"])
        metric = read(hfile["metric"])
        close(hfile)

        year_inds = 1:(argmin(temperature)-1) # volcano year
        acceptible_inds_lower = (temperature .> minimum(temperature[year_inds]))
        acceptible_inds_upper = (temperature .< maximum(temperature[year_inds]))
        acceptible_inds = acceptible_inds_lower .& acceptible_inds_upper
        year_inds = collect(eachindex(temperature))[acceptible_inds]

        global_mean_field = mean(historical_field[:, :, year_inds, month, :], dims = 1)[1, :, :, :]
        latitude_mean = mean(global_mean_field, dims = (2, 3))[:]
        latitude_std = std(global_mean_field, dims = (2, 3))[:]
    end

    include("../emulator.jl")

    rmetric = reshape(metric, (192, 96, 1, 1, 1))
    fmetric = reshape(metric, (192 * 96, 1))

    global_mean_basis = mean(reshape(emulator.basis, (192, 96, 1980)), dims = 1)[1, :, 1:1000]
    Σ = emulator_variance(emulator)
    mean_modes = mode_mean(emulator; modes = 1000)
    σs = [sqrt(global_mean_basis[i, :]' * (Σ * global_mean_basis[i,:])) for i in ProgressBar(1:96)]
    μs = [mean_modes' * global_mean_basis[i, :] for i in 1:96]
end

##
fig = Figure(; resolution) 
ga = fig[1, 1] = GridLayout()
ax = Axis(ga[1,1]; title = "Zonal Average (Data)", ylabel = "Temperature (K)", xlabel = "Latitude", common_options...)
lines!(ax, latitude[1:96], latitude_mean, color = :purple, label = "data")
band!(ax, latitude[1:96], latitude_mean .- 3 * latitude_std, latitude_mean .+ 3 * latitude_std, color = (:purple, 0.2))
ylims!(ax, 220, 310)
ax = Axis(ga[1,2]; title = "Zonal Average (Emulator)", xlabel = "Latitude", common_options...)
lines!(ax, latitude[1:96], μs, color = :blue, label = "data")
band!(ax, latitude[1:96], μs .- 3 * σs, μs .+ 3 * σs, color = (:blue, 0.2))
ylims!(ax, 220, 310)
gb = fig[2, 1] = GridLayout()
for (i, lat_index) in enumerate([1, 24, 48, 72, 96])
    latitude_string = @sprintf("Latitude %.2f", latitude[lat_index])
    if i == 1
        ax = Axis(gb[1,i]; xlabel = "Temperature (K)", ylabel = "Probability Density", title = latitude_string, common_options_2...)
    else
        ax = Axis(gb[1,i]; xlabel = "Temperature (K)", title = latitude_string, common_options_2...)
    end
    σ = σs[lat_index]
    μ = μs[lat_index]
    x = range(μ - 4σ, μ + 4σ, length = 100)
    y = gaussian.(x, μ, σ)
    data_samples = (use_bundle && has_ground_truth_bundle(bundle_path)) ? lat_samples[lat_index] : global_mean_field[lat_index, :, :][:]
    hist!(ax, data_samples, bins = 25, color = (:purple, 0.5), normalization = :pdf, label = "Data")
    lines!(ax, x, y, color = :blue, label = "Emulator")
    if i == 1
        axislegend(ax, position = :lt, labelsize = legend_ls)
    end
end
save(figure_directory * "figure_6_latitude_mean_model_emulator_locations_together.png", fig)
@info "Generated Figure 6"
