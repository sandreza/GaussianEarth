using CairoMakie, Statistics, ProgressBars, Printf

include("utils.jl")

hfile = h5open(save_directory * field_name * "_basis.hdf5", "r")
latitude = read(hfile["latitude"])
longitude = read(hfile["longitude"])
metric = read(hfile["metric"])
close(hfile)

month = 1
field = "tas"
eof_mode, temperature = concatenate_regression(field, ["historical"])
eofs = eof_mode[:,:, 1:45] # [1:48..., 50:50...]]

historical_field = common_array("historical", field)

year_inds = 1:(argmin(temperature)-1) # volcano year
acceptible_inds_lower = (temperature .> minimum(temperature[year_inds]))
acceptible_inds_upper = (temperature .< maximum(temperature[year_inds]))
acceptible_inds = acceptible_inds_lower .& acceptible_inds_upper
year_inds = collect(eachindex(temperature))[acceptible_inds]
##
κs = Float64[]
ρs = Float64[]
for mode_number in 1:1000
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

fig = Figure(resolution = (2000, 500))
inds = 1:45
lower_order_statistics = Vector{Float64}[]
for (j, mode_number) in enumerate(sort([gaussian_max, skewness_min, skewness_max, kurtosis_min, kurtosis_max]))
    ax = Axis(fig[1, j]; title = "Mode $mode_number")

    month_eof = eofs[mode_number, month:12:end, :]
    month_eof = month_eof[year_inds, :][:]
    μ = mean(month_eof)
    σ = std(month_eof)
    κ = kurtosis(month_eof)
    ρ = skewness(month_eof)
    push!(lower_order_statistics, [μ, σ, κ, ρ])
    x = range(μ - 4σ, μ + 4σ, length = 100)
    y = gaussian.(x, μ, σ)
    hist!(ax, month_eof, bins = 25, color = (:blue, 0.5), normalization = :pdf)
    lines!(ax, x, y, color = :red)
    xlims!(ax, μ - 4σ, μ + 4σ)
end
display(fig)
save(field * "_eof_gaussian.png", fig)

κs = Float64[]
ρs = Float64[]
for i in ProgressBar(1:192)
    for j in 1:96
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
    ii = (mode_number-1)÷96 + 1
    jj = (mode_number-1)%96 + 1

    lat = latitude[ii]
    lon = longitude[jj]
    title_string = @sprintf("Location (%.2f, %.2f)", lat, lon)
    ax = Axis(fig[2, j]; title = title_string)
    month_field = historical_field[ii, jj, year_inds, month, :][:]
    μ = mean(month_field)
    σ = std(month_field)
    κ = kurtosis(month_field)
    ρ = skewness(month_field)
    push!(lower_order_statistics, [μ, σ, κ, ρ])
    x = range(μ - 4σ, μ + 4σ, length = 100)
    y = gaussian.(x, μ, σ)
    hist!(ax, month_field, bins = 25, color = (:blue, 0.5), normalization = :pdf)
    lines!(ax, x, y, color = :red)
    xlims!(ax, μ - 4σ, μ + 4σ)
end
save(field * "_eof_points_gaussian.png", fig)
##

#=
month = 1
field = "hurs"
eof_mode, temperature = concatenate_regression(field, ["historical"])
eofs = eof_mode[:,:, :] # [1:48..., 50:50...]]

historical_field = common_array("historical", field; ensemble_members = size(eofs)[end])

year_inds = 1:(argmin(temperature)-1) # volcano year
acceptible_inds_lower = (temperature .> minimum(temperature[year_inds]))
acceptible_inds_upper = (temperature .< maximum(temperature[year_inds]))
acceptible_inds = acceptible_inds_lower .& acceptible_inds_upper
year_inds = collect(eachindex(temperature))[acceptible_inds]
##
κs = Float64[]
ρs = Float64[]
for mode_number in 1:1000
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

fig = Figure(resolution = (2000, 500))
inds = 1:45
lower_order_statistics = Vector{Float64}[]
for (j, mode_number) in enumerate(sort([gaussian_max, kurtosis_min, skewness_min, kurtosis_max, skewness_max]))
    ax = Axis(fig[1, j]; title = "Mode $mode_number")

    month_eof = eofs[mode_number, month:12:end, :]
    month_eof = month_eof[year_inds, :][:]
    μ = mean(month_eof)
    σ = std(month_eof)
    κ = kurtosis(month_eof)
    ρ = skewness(month_eof)
    push!(lower_order_statistics, [μ, σ, κ, ρ])
    x = range(μ - 4σ, μ + 4σ, length = 100)
    y = gaussian.(x, μ, σ)
    hist!(ax, month_eof, bins = 25, color = (:blue, 0.5), normalization = :pdf)
    lines!(ax, x, y, color = :red)
    xlims!(ax, μ - 4σ, μ + 4σ)
end
display(fig)
save(field * "_eof_gaussian.png", fig)

κs = Float64[]
ρs = Float64[]
for i in ProgressBar(1:192)
    for j in 1:96
        month_field = historical_field[i, j, year_inds, month, :][:]
        μ = mean(month_field)
        σ = std(month_field)
        κ = kurtosis(month_field)
        ρ = skewness(month_field)
        push!(κs, κ)
        push!(ρs, ρ)
    end
end

kurtosis_max = sortperm(κs)[2]
skewness_max = argmax(ρs)
kurtosis_min = argmin(κs)
skewness_min = argmin(ρs)
gaussian_max = argmin(abs.(κs) + abs.(ρs))


for (j, mode_number) in enumerate(sort([gaussian_max, kurtosis_min, skewness_min, kurtosis_max, skewness_max]))
    ii = (mode_number-1)÷96 + 1
    jj = (mode_number-1)%96 + 1

    ax = Axis(fig[2, j]; title = "Location ($ii, $jj)")
    month_field = historical_field[ii, jj, year_inds, month, :][:]
    μ = mean(month_field)
    σ = std(month_field)
    κ = kurtosis(month_field)
    ρ = skewness(month_field)
    push!(lower_order_statistics, [μ, σ, κ, ρ])
    x = range(μ - 4σ, μ + 4σ, length = 100)
    y = gaussian.(x, μ, σ)
    hist!(ax, month_field, bins = 25, color = (:blue, 0.5), normalization = :pdf)
    lines!(ax, x, y, color = :red)
    xlims!(ax, μ - 4σ, μ + 4σ)
end
save(field * "_eof_points_gaussian.png", fig)

=#