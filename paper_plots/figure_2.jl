ts = 40
xls = 40 
yls = 40
tls = 40
legend_ls = 35
resolution = (3000, 1000)
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
    modes = 1000 
    mean_field = mean(emulator; modes)
    variance_field = variance(emulator; modes)
    mean_modes = mode_mean(emulator; modes)
    variance_modes = mode_variance(emulator; modes)
end

modes = 1000 
mean_field = mean(emulator; modes)
variance_field = variance(emulator; modes)
mean_modes = mode_mean(emulator; modes)
variance_modes = mode_variance(emulator; modes)

κs = Float64[]
ρs = Float64[]
for mode_number in 1:modes
    month_eof = eofs[mode_number, month:12:end, :]
    month_eof = month_eof[year_inds, :][:]
    μ = mean(month_eof)
    σ = std(month_eof)
    κ = kurtosis(month_eof)
    ρ = skewness(month_eof)
    push!(κs, κ)
    push!(ρs, ρ)
end


kurtosis_max = argmax(κs)
skewness_max = argmax(ρs)
kurtosis_min = argmin(κs)
skewness_min = argmin(ρs)
gaussian_max = argmin(abs.(κs) + abs.(ρs))

fig = Figure(; resolution)
inds = 1:45
lower_order_statistics = Vector{Float64}[]
for (j, mode_number) in enumerate(sort([gaussian_max, skewness_min, skewness_max, kurtosis_min, kurtosis_max]))
    ax = Axis(fig[1, j]; title = "Mode $mode_number", xlabel = "Amplitude", ylabel = "Probability Density", common_options...)

    month_eof = eofs[mode_number, month:12:end, :]
    month_eof = month_eof[year_inds, :][:]
    μ = mean(month_eof)
    σ = std(month_eof)
    κ = kurtosis(month_eof)
    ρ = skewness(month_eof)
    push!(lower_order_statistics, [μ, σ, κ, ρ])
    x = range(μ - 4σ, μ + 4σ, length = 100)
    y = gaussian.(x, μ, σ)
    hist!(ax, month_eof, bins = 25, color = (:purple, 0.5), normalization = :pdf, label = "Data")
    # lines!(ax, x, y, color = :red)

    μ = mean_modes[mode_number]
    σ = sqrt(variance_modes[mode_number])
    y = gaussian.(x, μ, σ)
    lines!(ax, x, y, color = :blue, label = "Emulator")
    xlims!(ax, μ - 4σ, μ + 4σ)
    if j == 1 
        axislegend(ax, position = :lt, labelsize = legend_ls)
    end
end
display(fig)
# save(figure_directory * field * "_eof_gaussian_with_model.png", fig)

κs = Float64[]
ρs = Float64[]
for j in ProgressBar(1:96)
    for i in 1:192
        month_field = historical_field[i, j, year_inds, month, :][:]
        μ = mean(month_field)
        σ = std(month_field)
        κ = kurtosis(month_field)
        ρ = skewness(month_field)
        push!(κs, κ)
        push!(ρs, ρ)
    end
end

kurtosis_max = argmax(κs)
skewness_max = argmax(ρs)
kurtosis_min = argmin(κs)
skewness_min = argmin(ρs)
gaussian_max = argmin(abs.(κs) + abs.(ρs))


for (j, mode_number) in enumerate(sort([gaussian_max, kurtosis_min, skewness_min, kurtosis_max, skewness_max]))
    ii = (mode_number-1)%192 + 1
    jj = (mode_number-1)÷192 + 1

    lat = latitude[jj]
    lon = longitude[ii]
    titlestring = @sprintf("%.2fᵒ, %.2fᵒ", lon, lat)
    ax = Axis(fig[2, j]; title = titlestring, xlabel = "Temperature (K)", ylabel = "Probability Density", common_options...)
    month_field = historical_field[ii, jj, year_inds, month, :][:]
    μ = mean(month_field)
    σ = std(month_field)
    κ = kurtosis(month_field)
    ρ = skewness(month_field)
    push!(lower_order_statistics, [μ, σ, κ, ρ])
    x = range(μ - 4σ, μ + 4σ, length = 100)
    y = gaussian.(x, μ, σ)
    hist!(ax, month_field, bins = 25, color = (:purple, 0.5), normalization = :pdf, label = "Data")
    # lines!(ax, x, y, color = :red, label = "fit")
    y = gaussian.(x, mean_field[mode_number], sqrt(variance_field[mode_number]))
    lines!(ax, x, y, color = :blue, label = "Emulator")
    xlims!(ax, μ - 4σ, μ + 4σ)
end

save(figure_directory * "check_model_fit_modes_$modes.png", fig)
@info "Figure 2 generated"