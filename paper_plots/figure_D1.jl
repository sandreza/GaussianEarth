ts = 30
xls = 40 
yls = 40
tls = 30
legend_ls = 35
resolution = (400, 100) .* 6 
common_options = (; titlesize = ts, xlabelsize = xls, ylabelsize = yls, xticklabelsize = tls, yticklabelsize = tls)

##
month = 1
field = "tas"

hfile = h5open(save_directory * field * "_basis.hdf5", "r")
latitude = read(hfile["latitude"])
longitude = read(hfile["longitude"])
metric = read(hfile["metric"])
close(hfile)
fmetric = reshape(metric, (192*96, 1))

include("../emulator.jl")

##
global_zonal_average_basis = mean(reshape(emulator.basis, (192, 96, 1980)), dims = 1)[1, :, 1:1000]
global_mean_basis = sum(emulator.basis .* fmetric, dims = 1)[1:1000]
global_mean_basis_upper = sum(reshape(emulator.basis .* fmetric, (192, 96, 1980))[:, 49:end, :], dims = (1, 2))[1:1000] * 2
global_mean_basis_lower = sum(reshape(emulator.basis .* fmetric, (192, 96, 1980))[:, 1:48, :], dims = (1, 2))[1:1000] * 2
min_temp = 1.0544952f0
max_temp = 1.0720718f0
emulator.global_mean_temperature[1] = min_temp 
emulator.month[1] = month

Σ = emulator_variance(emulator)
mean_modes = mode_mean(emulator; modes = 1000)
σs = [sqrt(global_zonal_average_basis[i, :]' * (Σ * global_zonal_average_basis[i,:])) for i in ProgressBar(1:96)]
μs = [mean_modes' * global_zonal_average_basis[i, :] for i in 1:96]

cold_mus = Float32[]
cold_sigmas = Float32[]
cold_mus_lat = Vector{Float32}[]
cold_sigmas_lat = Vector{Float32}[]
cold_mus_upper = Float32[]
cold_sigmas_upper = Float32[]
cold_mus_lower = Float32[]
cold_sigmas_lower = Float32[]
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
    σ_lower = sqrt(global_mean_basis_lower' * (Σ * global_mean_basis_lower))
    μ_lower = mean_modes' * global_mean_basis_lower
    push!(cold_mus_lower, μ_lower)
    push!(cold_sigmas_lower, σ_lower)
    σ_upper = sqrt(global_mean_basis_upper' * (Σ * global_mean_basis_upper))
    μ_upper = mean_modes' * global_mean_basis_upper
    push!(cold_mus_upper, μ_upper)
    push!(cold_sigmas_upper, σ_upper)
end

hot_mus = Float32[]
hot_sigmas = Float32[]
hot_mus_lat = Vector{Float32}[]
hot_sigmas_lat = Vector{Float32}[]
hot_mus_upper = Float32[]
hot_sigmas_upper = Float32[]
hot_mus_lower = Float32[]
hot_sigmas_lower = Float32[]
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
    σ_lower = sqrt(global_mean_basis_lower' * (Σ * global_mean_basis_lower))
    μ_lower = mean_modes' * global_mean_basis_lower
    push!(hot_mus_lower, μ_lower)
    push!(hot_sigmas_lower, σ_lower)
    σ_upper = sqrt(global_mean_basis_upper' * (Σ * global_mean_basis_upper))
    μ_upper = mean_modes' * global_mean_basis_upper
    push!(hot_mus_upper, μ_upper)
    push!(hot_sigmas_upper, σ_upper)
end
 
month_string = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
xticks = (1:12, month_string)

fig = Figure(; resolution )
ax = Axis(fig[1, 1]; title = "Global Average (Emulator)", xlabel = "Month", ylabel = "Temperature (K)", xticks, common_options...)
lines!(ax, 1:12, cold_mus, color = :blue, label = "Coolest")
band!(ax, 1:12, cold_mus .- 3 * cold_sigmas, cold_mus .+ 3 * cold_sigmas, color = (:blue, 0.2))
lines!(ax, 1:12, hot_mus, color = :orange, label = "Warmest")
band!(ax, 1:12, hot_mus .- 3 * hot_sigmas, hot_mus .+ 3 * hot_sigmas, color = (:orange, 0.2))
xlims!(ax, (0.6, 12.4))
ylims!(ax, (280, 302))
axislegend(ax, position = :lt, labelsize = legend_ls)
maximum(cold_mus) - minimum(cold_mus)
maximum(hot_mus) - minimum(hot_mus)

ax = Axis(fig[1, 2]; title = "Northern Hemisphere Average (Emulator)", xlabel = "Month", ylabel = "Temperature (K)", xticks, common_options...)
lines!(ax, 1:12, cold_mus_upper, color = :blue, label = "Coolest")
band!(ax, 1:12, cold_mus_upper .- 3 * cold_sigmas_upper, cold_mus_upper .+ 3 * cold_sigmas_upper, color = (:blue, 0.2))
lines!(ax, 1:12, hot_mus_upper, color = :orange, label = "Warmest")
band!(ax, 1:12, hot_mus_upper .- 3 * hot_sigmas_upper, hot_mus_upper .+ 3 * hot_sigmas_upper, color = (:orange, 0.2))
hideydecorations!(ax, grid = false)
xlims!(ax, (0.6, 12.4))
ylims!(ax, (280, 302))
maximum(cold_mus_upper) - minimum(cold_mus_upper)
maximum(hot_mus_upper) - minimum(hot_mus_upper)

ax = Axis(fig[1, 3]; title = "Southern Hemisphere Average (Emulator)", xlabel = "Month", ylabel = "Temperature (K)", xticks, common_options...)
lines!(ax, 1:12, cold_mus_lower, color = :blue, label = "Coolest")
band!(ax, 1:12, cold_mus_lower .- 3 * cold_sigmas_lower, cold_mus_lower .+ 3 * cold_sigmas_lower, color = (:blue, 0.2))
lines!(ax, 1:12, hot_mus_lower, color = :orange, label = "Warmest")
band!(ax, 1:12, hot_mus_lower .- 3 * hot_sigmas_lower, hot_mus_lower .+ 3 * hot_sigmas_lower, color = (:orange, 0.2))
hideydecorations!(ax, grid = false)
xlims!(ax, (0.6, 12.4))
ylims!(ax, (280, 302))
maximum(cold_mus_lower) - minimum(cold_mus_lower)
maximum(hot_mus_lower) - minimum(hot_mus_lower)

save(figure_directory * "figure_D1_climate_change_shifts_monthly_with_upper_lower.png", fig)
