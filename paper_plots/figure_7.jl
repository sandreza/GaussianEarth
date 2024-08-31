ts = 40
xls = 40 
yls = 40
tls = 40
legend_ls = 35
resolution = (1800, 600)
common_options = (; titlesize = ts, xlabelsize = xls, ylabelsize = yls, xticklabelsize = tls, yticklabelsize = tls)

if process_data
    month = 1
    field = "tas"

    hfile = h5open(save_directory * field * "_basis.hdf5", "r")
    latitude = read(hfile["latitude"])
    longitude = read(hfile["longitude"])
    metric = read(hfile["metric"])
    close(hfile)

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

    include("emulator.jl")

    global_mean_field = mean(historical_field[:, :, year_inds, month, :], dims = 1)[1, :, :, :]
    latitude_mean = mean(global_mean_field, dims = (2, 3))[:]
    latitude_std = std(global_mean_field, dims = (2, 3))[:]

    rmetric = reshape(metric, (192, 96, 1, 1, 1))
    fmetric = reshape(metric, (192 * 96, 1))

    global_mean_basis = mean(reshape(emulator.basis, (192, 96, 1980)), dims = 1)[1, :, 1:1000]
    Σ = emulator_variance(emulator)
    mean_modes = mode_mean(emulator; modes = 1000)
    σs = [sqrt(global_mean_basis[i, :]' * (Σ * global_mean_basis[i,:])) for i in ProgressBar(1:96)]
    μs = [mean_modes' * global_mean_basis[i, :] for i in 1:96]
end

fig = Figure(; resolution) 
ax = Axis(fig[1,1]; title = "Latitude Mean and Variance (Data)", xlabel = "Latitude", ylabel = "Temperature (K)", common_options...)
lines!(ax, latitude[1:96], latitude_mean, color = :purple, label = "data")
band!(ax, latitude[1:96], latitude_mean .- 3 * latitude_std, latitude_mean .+ 3 * latitude_std, color = (:purple, 0.2))
ax = Axis(fig[1,2]; title = "Latitude Mean and Variance (Emulator)", xlabel = "Latitude", ylabel = "Temperature (K)", common_options...)
lines!(ax, latitude[1:96], μs, color = :blue, label = "data")
band!(ax, latitude[1:96], μs .- 3 * σs, μs .+ 3 * σs, color = (:blue, 0.2))
display(fig)
save(figure_directory * "latitude_mean_model_emulator.png", fig)