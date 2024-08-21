using CairoMakie, Statistics, ProgressBars

include("utils.jl")
include("emulator.jl")

#=
month = 1
field = "tas"
eof_mode, temperature = concatenate_regression(field, ["historical"])
eofs = eof_mode[:,:, 1:45] # [1:48..., 50:50...]]

historical_field = common_array("historical", field)
hfile = h5open(save_directory * field * "_basis.hdf5", "r")
metric = read(hfile["metric"])
close(hfile)

year_inds = 1:(argmin(temperature)-1) # volcano year
acceptible_inds_lower = (temperature .> minimum(temperature[year_inds]))
acceptible_inds_upper = (temperature .< maximum(temperature[year_inds]))
acceptible_inds = acceptible_inds_lower .& acceptible_inds_upper
year_inds = collect(eachindex(temperature))[acceptible_inds]
=#
##
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
    hist!(ax, month_eof, bins = 25, color = (:purple, 0.5), normalization = :pdf)
    lines!(ax, x, y, color = :red)

    μ = mean_modes[mode_number]
    σ = sqrt(variance_modes[mode_number])
    y = gaussian.(x, μ, σ)
    lines!(ax, x, y, color = :blue)
    xlims!(ax, μ - 4σ, μ + 4σ)
end
display(fig)
save(field * "_eof_gaussian_with_model.png", fig)

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

    ax = Axis(fig[2, j]; title = "Location ($ii, $jj)")
    month_field = historical_field[ii, jj, year_inds, month, :][:]
    μ = mean(month_field)
    σ = std(month_field)
    κ = kurtosis(month_field)
    ρ = skewness(month_field)
    push!(lower_order_statistics, [μ, σ, κ, ρ])
    x = range(μ - 4σ, μ + 4σ, length = 100)
    y = gaussian.(x, μ, σ)
    hist!(ax, month_field, bins = 25, color = (:purple, 0.5), normalization = :pdf)
    lines!(ax, x, y, color = :red, label = "fit")
    y = gaussian.(x, mean_field[mode_number], sqrt(variance_field[mode_number]))
    lines!(ax, x, y, color = :blue, label = "model")
    xlims!(ax, μ - 4σ, μ + 4σ)
end

save("check_model_fit_modes_$modes.png", fig)

##   
modes = 10
mean_field = mean(emulator; modes)
variance_field = variance(emulator; modes)
mean_modes = mode_mean(emulator; modes)
variance_modes = mode_variance(emulator; modes)

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
    hist!(ax, month_eof, bins = 25, color = (:purple, 0.5), normalization = :pdf)
    lines!(ax, x, y, color = :red)

    μ = mean_modes[mode_number]
    σ = sqrt(variance_modes[mode_number])
    y = gaussian.(x, μ, σ)
    lines!(ax, x, y, color = :blue)
    xlims!(ax, μ - 4σ, μ + 4σ)
end
display(fig)
save(field * "_eof_gaussian_with_model.png", fig)

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

    ax = Axis(fig[2, j]; title = "Location ($ii, $jj)")
    month_field = historical_field[ii, jj, year_inds, month, :][:]
    μ = mean(month_field)
    σ = std(month_field)
    κ = kurtosis(month_field)
    ρ = skewness(month_field)
    push!(lower_order_statistics, [μ, σ, κ, ρ])
    x = range(μ - 4σ, μ + 4σ, length = 100)
    y = gaussian.(x, μ, σ)
    hist!(ax, month_field, bins = 25, color = (:purple, 0.5), normalization = :pdf)
    lines!(ax, x, y, color = :red, label = "fit")
    y = gaussian.(x, mean_field[mode_number], sqrt(variance_field[mode_number]))
    lines!(ax, x, y, color = :blue, label = "model")
    xlims!(ax, μ - 4σ, μ + 4σ)
end

save("check_model_fit_modes_$modes.png", fig)

##
modes = 100
mean_field = mean(emulator; modes)
variance_field = variance(emulator; modes)
mean_modes = mode_mean(emulator; modes)
variance_modes = mode_variance(emulator; modes)

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
    hist!(ax, month_eof, bins = 25, color = (:purple, 0.5), normalization = :pdf)
    lines!(ax, x, y, color = :red)

    μ = mean_modes[mode_number]
    σ = sqrt(variance_modes[mode_number])
    y = gaussian.(x, μ, σ)
    lines!(ax, x, y, color = :blue)
    xlims!(ax, μ - 4σ, μ + 4σ)
end
display(fig)
save(field * "_eof_gaussian_with_model.png", fig)

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

    ax = Axis(fig[2, j]; title = "Location ($ii, $jj)")
    month_field = historical_field[ii, jj, year_inds, month, :][:]
    μ = mean(month_field)
    σ = std(month_field)
    κ = kurtosis(month_field)
    ρ = skewness(month_field)
    push!(lower_order_statistics, [μ, σ, κ, ρ])
    x = range(μ - 4σ, μ + 4σ, length = 100)
    y = gaussian.(x, μ, σ)
    hist!(ax, month_field, bins = 25, color = (:purple, 0.5), normalization = :pdf)
    lines!(ax, x, y, color = :red, label = "fit")
    y = gaussian.(x, mean_field[mode_number], sqrt(variance_field[mode_number]))
    lines!(ax, x, y, color = :blue, label = "model")
    xlims!(ax, μ - 4σ, μ + 4σ)
end

save("check_model_fit_modes_$modes.png", fig)

##
fig = Figure() 
ax = Axis(fig[1,1]; title = "Global Mean (no metric)")
rmetric = reshape(metric, (192, 96, 1, 1, 1))
fmetric = reshape(metric, (192 * 96, 1))
global_mean_field = mean(historical_field[:, :, year_inds, month, :], dims = (1, 2))[:]
hackstimate = mean(std(historical_field[:,:, year_inds, month, :], dims = (3,4)))
global_mean_basis = mean(emulator.basis, dims = 1)[1:1000]
Σ = emulator_variance(emulator)
mean_modes = mode_mean(emulator; modes = 1000)
σ = sqrt(global_mean_basis' * (Σ * global_mean_basis))
μ = mean_modes' * global_mean_basis
std(global_mean_field )
x = range(μ - 4σ, μ + 4σ, length = 100)
y = gaussian.(x, μ, σ)
hist!(ax, global_mean_field, bins = 25, color = (:purple, 0.5), normalization = :pdf)
lines!(ax, x, y, color = :blue, label = "emulator")
display(fig)
save("global_mean_model_emulator.png", fig)
##
fig = Figure() 
ax = Axis(fig[1,1]; title = "Global Mean (metric)")
rmetric = reshape(metric, (192, 96, 1, 1, 1))
fmetric = reshape(metric, (192 * 96, 1))
global_mean_field = sum(historical_field[:, :, year_inds, month, :] .* rmetric, dims = (1, 2))[:]
hackstimate = sum(std(historical_field[:,:, year_inds, month, :] , dims = (3,4)) .* metric)
global_mean_basis = sum(emulator.basis .* fmetric, dims = 1)[1:1000]
Σ = emulator_variance(emulator)
mean_modes = mode_mean(emulator; modes = 1000)
σ = sqrt(global_mean_basis' * (Σ * global_mean_basis))
μ = mean_modes' * global_mean_basis
std(global_mean_field )
x = range(μ - 4σ, μ + 4σ, length = 100)
y = gaussian.(x, μ, σ)
hist!(ax, global_mean_field, bins = 25, color = (:purple, 0.5), normalization = :pdf)
lines!(ax, x, y, color = :blue, label = "emulator")
display(fig)
save("global_mean_model_emulator_2.png", fig)