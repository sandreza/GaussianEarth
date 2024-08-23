using CairoMakie, Statistics, ProgressBars
using Printf

include("utils.jl")
include("emulator.jl")
##
month = 1
field = "tas"

hfile = h5open(save_directory * field * "_basis.hdf5", "r")
latitude = read(hfile["latitude"])
longitude = read(hfile["longitude"])
metric = read(hfile["metric"])
close(hfile)
fmetric = reshape(metric, (192*96, 1))

##
global_zonal_average_basis = mean(reshape(emulator.basis, (192, 96, 1980)), dims = 1)[1, :, 1:1000]
global_mean_basis = sum(emulator.basis .* fmetric, dims = 1)[1:1000]
min_temp = 1.0544952f0
max_temp = 1.0720718f0
emulator.global_mean_temperature[1] = min_temp 
emulator.month[1] = month

Σ = emulator_variance(emulator)
mean_modes = mode_mean(emulator; modes = 1000)
σs = [sqrt(global_zonal_average_basis[i, :]' * (Σ * global_zonal_average_basis[i,:])) for i in ProgressBar(1:96)]
μs = [mean_modes' * global_zonal_average_basis[i, :] for i in 1:96]

fig = Figure()
ax = Axis(fig[1, 1]; title = "Zonal Average", xlabel = "Latitude", ylabel = "Temperature (K)")
lines!(ax, latitude[1:96], μs, color = :blue, label = "Coolest")
band!(ax, latitude[1:96], μs .- 3 * σs, μs .+ 3 * σs, color = (:blue, 0.2))

emulator.global_mean_temperature[1] = max_temp 
emulator.month[1] = month
Σ = emulator_variance(emulator)
mean_modes = mode_mean(emulator; modes = 1000)
σs = [sqrt(global_zonal_average_basis[i, :]' * (Σ * global_zonal_average_basis[i,:])) for i in ProgressBar(1:96)]
μs = [mean_modes' * global_zonal_average_basis[i, :] for i in 1:96]

lines!(ax, latitude[1:96], μs, color = :orange, label = "Warmest")
band!(ax, latitude[1:96], μs .- 3 * σs, μs .+ 3 * σs, color = (:orange, 0.2))

axislegend(ax, position = :rt)
save("climage_change_shifts.png", fig)



ax = Axis(fig[1, 2]; title = "Global Average", xlabel = "Time (years)", ylabel = "PDF")
emulator.global_mean_temperature[1] = min_temp 
emulator.month[1] = month
Σ = emulator_variance(emulator)
mean_modes = mode_mean(emulator; modes = 1000)
σ = sqrt(global_mean_basis' * (Σ * global_mean_basis)) 
μ = mean_modes' * global_mean_basis
x = range(μ - 4σ, μ + 4σ, length = 100)
y = gaussian.(x, μ, σ)
lines!(ax, x, y, color = :blue, label = "Coolest")

emulator.global_mean_temperature[1] = max_temp
emulator.month[1] = month
Σ = emulator_variance(emulator)
mean_modes = mode_mean(emulator; modes = 1000)
σ = sqrt(global_mean_basis' * (Σ * global_mean_basis))
μ = mean_modes' * global_mean_basis
x = range(μ - 4σ, μ + 4σ, length = 100)
y = gaussian.(x, μ, σ)
lines!(ax, x, y, color = :orange, label = "Warmest")

save("climage_change_shifts_with_global_mean.png", fig)

##

cold_mus = Float32[]
cold_sigmas = Float32[]
cold_mus_lat = Vector{Float32}[]
cold_sigmas_lat = Vector{Float32}[]
for month in 1:12
    emulator.global_mean_temperature[1] = min_temp 
    emulator.month[1] = month
    Σ = emulator_variance(emulator)
    mean_modes = mode_mean(emulator; modes = 1000)
    σ = sqrt(global_mean_basis' * (Σ * global_mean_basis)) 
    μ = mean_modes' * global_mean_basis
    push!(cold_mus, μ)
    push!(cold_sigmas, σ)
    σs_lat = [sqrt(global_zonal_average_basis[i, :]' * (Σ * global_zonal_average_basis[i,:])) for i in ProgressBar(1:96)]
    μs_lat = [mean_modes' * global_zonal_average_basis[i, :] for i in 1:96]
    push!(cold_mus_lat, μs_lat)
    push!(cold_sigmas_lat, σs_lat)
end

hot_mus = Float32[]
hot_sigmas = Float32[]
hot_mus_lat = Vector{Float32}[]
hot_sigmas_lat = Vector{Float32}[]
for month in 1:12
    emulator.global_mean_temperature[1] = max_temp 
    emulator.month[1] = month
    Σ = emulator_variance(emulator)
    mean_modes = mode_mean(emulator; modes = 1000)
    σ = sqrt(global_mean_basis' * (Σ * global_mean_basis)) 
    μ = mean_modes' * global_mean_basis
    push!(hot_mus, μ)
    push!(hot_sigmas, σ)
    σs_lat = [sqrt(global_zonal_average_basis[i, :]' * (Σ * global_zonal_average_basis[i,:])) for i in ProgressBar(1:96)]
    μs_lat = [mean_modes' * global_zonal_average_basis[i, :] for i in 1:96]
    push!(hot_mus_lat, μs_lat)
    push!(hot_sigmas_lat, σs_lat)
end

month_string = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
xticks = (1:12, month_string)
fig = Figure()
ax = Axis(fig[1, 1]; title = "Global Average (Emulator)", xlabel = "Month", ylabel = "Temperature (K)", xticks)
lines!(ax, 1:12, cold_mus, color = :blue, label = "Coolest")
band!(ax, 1:12, cold_mus .- 3 * cold_sigmas, cold_mus .+ 3 * cold_sigmas, color = (:blue, 0.2))
lines!(ax, 1:12, hot_mus, color = :orange, label = "Warmest")
band!(ax, 1:12, hot_mus .- 3 * hot_sigmas, hot_mus .+ 3 * hot_sigmas, color = (:orange, 0.2))
xlims!(ax, (0.6, 12.4))
axislegend(ax, position = :rt)
save("climage_change_shifts_monthly.png", fig)


fig = Figure(resolution = (2000, 400))
ax = Axis(fig[1, 1]; title = "Global Average (Emulator)", xlabel = "Month", ylabel = "Temperature (K)", xticks)
lines!(ax, 1:12, cold_mus, color = :blue, label = "Coolest")
band!(ax, 1:12, cold_mus .- 3 * cold_sigmas, cold_mus .+ 3 * cold_sigmas, color = (:blue, 0.2))
lines!(ax, 1:12, hot_mus, color = :orange, label = "Warmest")
band!(ax, 1:12, hot_mus .- 3 * hot_sigmas, hot_mus .+ 3 * hot_sigmas, color = (:orange, 0.2))
xlims!(ax, (0.6, 12.4))
axislegend(ax, position = :lt)
for (jj, month) in enumerate([1, 4, 7, 10])
    ax = Axis(fig[1, jj + 1], xlabel = "Latitude", ylabel = "Temperature (K)", title = "Zonal Average Month " * month_string[month] * " (Emulator)")
    lines!(ax, latitude[1:96], cold_mus_lat[month], color = :blue, label = "Coolest")
    band!(ax, latitude[1:96], cold_mus_lat[month] .- 3 * cold_sigmas_lat[month], cold_mus_lat[month] .+ 3 * cold_sigmas_lat[month], color = (:blue, 0.2))
    lines!(ax, latitude[1:96], hot_mus_lat[month], color = :orange, label = "Warmest")
    band!(ax, latitude[1:96], hot_mus_lat[month] .- 3 * hot_sigmas_lat[month], hot_mus_lat[month] .+ 3 * hot_sigmas_lat[month], color = (:orange, 0.2))
end
save("climage_change_shifts_monthly_with_lat.png", fig)




