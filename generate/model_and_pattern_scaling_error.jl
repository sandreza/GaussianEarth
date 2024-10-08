using NCDatasets, LinearAlgebra, Statistics, HDF5, ProgressBars, CairoMakie
include("utils.jl")
##
save_directory = "/net/fs06/d3/sandre/GaussianEarthData/"
data_directory = "/net/fs06/d3/mgeo/CMIP6/interim/"
scenario_directories = readdir(data_directory)
current_path = joinpath(data_directory, scenario_directories[1])
variable_directories = readdir(current_path)
##
field_name = "tas" # ["tas", "hurs"]
colors = [:red4, :red, :indigo, :magenta3]
hfile = h5open(save_directory * field_name * "_basis.hdf5", "r")
latitude = read(hfile["latitude"])
longitude = read(hfile["longitude"])
metric = read(hfile["metric"])
close(hfile)
sqrt_f_metric = sqrt.(reshape(metric, 192 * 96))

hfile = h5open(save_directory * field_name * "_mean_regression.hdf5", "r")
regression_coefficients = read(hfile["regression_coefficients 1"])
linear_coefficients = zeros(Float32, size(regression_coefficients)..., 12)
for month in ProgressBar(1:12)
    linear_coefficients[:, :, month] .= read(hfile["regression_coefficients $month"])
end
close(hfile)

Φ = eof_basis(field_name)

averaged_coefficients = mean(linear_coefficients, dims = 3)[:, :, 1]

eof_model10 = zeros(Float32, size(Φ)[1]..., 2)
eof_model100 = zeros(Float32, size(Φ)[1]..., 2)
eof_model1000 = zeros(Float32, size(Φ)[1]..., 2)

for i in ProgressBar(1:10)
    eof_model10[:, 1] .+= Φ[:, i] * averaged_coefficients[i, 1]
    eof_model10[:, 2] .+= Φ[:, i] * averaged_coefficients[i, 2]
end
for i in ProgressBar(1:100)
    eof_model100[:, 1] .+= Φ[:, i] * averaged_coefficients[i, 1]
    eof_model100[:, 2] .+= Φ[:, i] * averaged_coefficients[i, 2]
end
for i in ProgressBar(1:1000)
    eof_model1000[:, 1] .+= Φ[:, i] * averaged_coefficients[i, 1]
    eof_model1000[:, 2] .+= Φ[:, i] * averaged_coefficients[i, 2]
end
##
filename  = field_name * "_pattern_scaling.hdf5"
hfile = h5open(save_directory * filename, "r")
linear_fit_yearly = read(hfile["linear fit yearly"])
close(hfile)
scenarios = ["historical", "ssp585", "ssp119", "ssp245"]
temperatures = []
fields = []
for scenario in scenarios
    temperature = regression_variable(scenario)
    a, b = ensemble_averaging(scenario, field_name; ensemble_members = 29)
    push!(temperatures, temperature)
    push!(fields, b)
end

##
i = 2
Ts = temperatures[i]
field = fields[i]
error_list10 = [norm((reshape(field[:, :, j], size(Φ)[1]) .- (eof_model10[:, 1] .+ Ts[j] * eof_model10[:, 2]))[:] .* sqrt_f_metric) for j in eachindex(Ts)]
ymax = maximum(error_list10)

yrange = (0, ymax)
fig_tas = Figure(resolution = (2000, 400)) 
ts = [1850:2014, 2015:2100, 2015:2100, 2015:2100]
ax = Axis(fig_tas[1,1]; title = "Pattern Scaling", xlabel = "Year", ylabel = "Temperature RMSE (K)")
for (i, field) in enumerate(fields)
    Ts = temperatures[i]
    error_list = [norm((field[:, :, j] .- (linear_fit_yearly[:, :, 1] .+ Ts[j] * linear_fit_yearly[:, :, 2]))[:] .* sqrt_f_metric) for j in eachindex(Ts)]
    lines!(ax, ts[i], error_list, label = scenarios[i], color = colors[i])
    ylims!(ax, yrange...)
end
axislegend(ax; position = :lt)

ax = Axis(fig_tas[1,4]; title = "10 Modes", xlabel = "Year")
for (i, field) in enumerate(fields)
    Ts = temperatures[i]
    error_list10 = [norm((reshape(field[:, :, j], size(Φ)[1]) .- (eof_model10[:, 1] .+ Ts[j] * eof_model10[:, 2]) )[:] .* sqrt_f_metric) for j in eachindex(Ts)]
    lines!(ax, ts[i], error_list10, label = scenarios[i] * " 10 modes", color = colors[i])
    ylims!(ax, yrange...)
end
hideydecorations!(ax, grid = false)

ax = Axis(fig_tas[1,3]; title = "100 Modes", xlabel = "Year")
for (i, field) in enumerate(fields)
    Ts = temperatures[i]
    error_list100 = [norm((reshape(field[:, :, j], size(Φ)[1]) .- (eof_model100[:, 1] .+ Ts[j] * eof_model100[:, 2]) )[:] .* sqrt_f_metric) for j in eachindex(Ts)]
    lines!(ax, ts[i], error_list100, label = scenarios[i] * " 100 modes", color = colors[i])
    ylims!(ax, yrange...)
end
hideydecorations!(ax, grid = false)

ax = Axis(fig_tas[1,2]; title = "1000 Modes", xlabel = "Year")
for (i, field) in enumerate(fields)
    Ts = temperatures[i]
    error_list1000 = [norm((reshape(field[:, :, j], size(Φ)[1]) .- (eof_model1000[:, 1] .+ Ts[j] * eof_model1000[:, 2]) )[:] .* sqrt_f_metric) for j in eachindex(Ts)]
    lines!(ax, ts[i], error_list1000, label = scenarios[i] * " 1000 modes", color = colors[i])
    ylims!(ax, yrange...)
end
hideydecorations!(ax, grid = false)
save(field_name * "_pattern_scaling_and_model_errors_3_metric.png", fig_tas)

##
field_name = "hurs" # ["tas", "hurs"]

hfile = h5open(save_directory * field_name * "_basis.hdf5", "r")
latitude = read(hfile["latitude"])
longitude = read(hfile["longitude"])
metric = read(hfile["metric"])
close(hfile)
sqrt_f_metric = sqrt.(reshape(metric, 192 * 96))

hfile = h5open(save_directory * field_name * "_mean_regression.hdf5", "r")
regression_coefficients = read(hfile["regression_coefficients 1"])
linear_coefficients = zeros(Float32, size(regression_coefficients)..., 12)
for month in ProgressBar(1:12)
    linear_coefficients[:, :, month] .= read(hfile["regression_coefficients $month"])
end
close(hfile)

Φ = eof_basis(field_name)

averaged_coefficients = mean(linear_coefficients, dims = 3)[:, :, 1]

eof_model10 = zeros(Float32, size(Φ)[1]..., 2)
eof_model100 = zeros(Float32, size(Φ)[1]..., 2)
eof_model1000 = zeros(Float32, size(Φ)[1]..., 2)

for i in ProgressBar(1:10)
    eof_model10[:, 1] .+= Φ[:, i] * averaged_coefficients[i, 1]
    eof_model10[:, 2] .+= Φ[:, i] * averaged_coefficients[i, 2]
end
for i in ProgressBar(1:100)
    eof_model100[:, 1] .+= Φ[:, i] * averaged_coefficients[i, 1]
    eof_model100[:, 2] .+= Φ[:, i] * averaged_coefficients[i, 2]
end
for i in ProgressBar(1:1000)
    eof_model1000[:, 1] .+= Φ[:, i] * averaged_coefficients[i, 1]
    eof_model1000[:, 2] .+= Φ[:, i] * averaged_coefficients[i, 2]
end
##
filename  = field_name * "_pattern_scaling.hdf5"
hfile = h5open(save_directory * filename, "r")
linear_fit_yearly = read(hfile["linear fit yearly"])
close(hfile)
scenarios = ["historical", "ssp585", "ssp119", "ssp245"]
temperatures = []
fields = []
for scenario in scenarios
    temperature = regression_variable(scenario)
    a, b = ensemble_averaging(scenario, field_name; ensemble_members = 29)
    push!(temperatures, temperature)
    push!(fields, b)
end

##
##  Dict("historical" => :red4, "ssp585" => :red, "ssp245" => :magenta3, "ssp119" => :indigo)
colors = [:red4, :red, :indigo, :magenta3]
i = 2
Ts = temperatures[i]
field = fields[i]
error_list10 = [norm((reshape(field[:, :, j], size(Φ)[1]) .- (eof_model10[:, 1] .+ Ts[j] * eof_model10[:, 2]))[:] .* sqrt_f_metric) for j in eachindex(Ts)]
ymax = maximum(error_list10)

yrange = (0, ymax)
fig = Figure(resolution = (2000, 400)) 
ts = [1850:2014, 2015:2100, 2015:2100, 2015:2100]
ax = Axis(fig[1,1]; title = "Pattern Scaling", xlabel = "Year", ylabel = "Relative Humidity RMSE (Percent)")
for (i, field) in enumerate(fields)
    Ts = temperatures[i]
    error_list = [norm((field[:, :, j] .- (linear_fit_yearly[:, :, 1] .+ Ts[j] * linear_fit_yearly[:, :, 2]))[:] .* sqrt_f_metric) for j in eachindex(Ts)]
    lines!(ax, ts[i], error_list, label = scenarios[i], color = colors[i])
    ylims!(ax, yrange...)
end
axislegend(ax; position = :lt)

ax = Axis(fig[1,4]; title = "10 Modes", xlabel = "Year")
for (i, field) in enumerate(fields)
    Ts = temperatures[i]
    error_list10 = [norm((reshape(field[:, :, j], size(Φ)[1]) .- (eof_model10[:, 1] .+ Ts[j] * eof_model10[:, 2]) )[:] .* sqrt_f_metric) for j in eachindex(Ts)]
    lines!(ax, ts[i], error_list10, label = scenarios[i] * " 10 modes", color = colors[i])
    ylims!(ax, yrange...)
end
hideydecorations!(ax, grid = false)

ax = Axis(fig[1,3]; title = "100 Modes", xlabel = "Year")
for (i, field) in enumerate(fields)
    Ts = temperatures[i]
    error_list100 = [norm((reshape(field[:, :, j], size(Φ)[1]) .- (eof_model100[:, 1] .+ Ts[j] * eof_model100[:, 2]) )[:] .* sqrt_f_metric) for j in eachindex(Ts)]
    lines!(ax, ts[i], error_list100, label = scenarios[i] * " 100 modes", color = colors[i])
    ylims!(ax, yrange...)
end
hideydecorations!(ax, grid = false)

ax = Axis(fig[1,2]; title = "1000 Modes", xlabel = "Year")
for (i, field) in enumerate(fields)
    Ts = temperatures[i]
    error_list1000 = [norm((reshape(field[:, :, j], size(Φ)[1]) .- (eof_model1000[:, 1] .+ Ts[j] * eof_model1000[:, 2]) )[:] .* sqrt_f_metric) for j in eachindex(Ts)]
    lines!(ax, ts[i], error_list1000, label = scenarios[i] * " 1000 modes", color = colors[i])
    ylims!(ax, yrange...)
end
hideydecorations!(ax, grid = false)
save(field_name * "_pattern_scaling_and_model_errors_3_metric.png", fig)

##
fig_hurs_tas = Figure(resolution = (2000, 800))
ax = Axis(fig_hurs_tas[2,1]; title = "Pattern Scaling", xlabel = "Year", ylabel = "Relative Humidity RMSE (Percent)")
for (i, field) in enumerate(fields)
    Ts = temperatures[i]
    error_list = [norm((field[:, :, j] .- (linear_fit_yearly[:, :, 1] .+ Ts[j] * linear_fit_yearly[:, :, 2]))[:] .* sqrt_f_metric) for j in eachindex(Ts)]
    lines!(ax, ts[i], error_list, label = scenarios[i], color = colors[i])
    ylims!(ax, yrange...)
end
# axislegend(ax; position = :lt)

ax = Axis(fig_hurs_tas[2,4]; title = "10 Modes", xlabel = "Year")
for (i, field) in enumerate(fields)
    Ts = temperatures[i]
    error_list10 = [norm((reshape(field[:, :, j], size(Φ)[1]) .- (eof_model10[:, 1] .+ Ts[j] * eof_model10[:, 2]) )[:] .* sqrt_f_metric) for j in eachindex(Ts)]
    lines!(ax, ts[i], error_list10, label = scenarios[i] * " 10 modes", color = colors[i])
    ylims!(ax, yrange...)
end
hideydecorations!(ax, grid = false)

ax = Axis(fig_hurs_tas[2,3]; title = "100 Modes", xlabel = "Year")
for (i, field) in enumerate(fields)
    Ts = temperatures[i]
    error_list100 = [norm((reshape(field[:, :, j], size(Φ)[1]) .- (eof_model100[:, 1] .+ Ts[j] * eof_model100[:, 2]) )[:] .* sqrt_f_metric) for j in eachindex(Ts)]
    lines!(ax, ts[i], error_list100, label = scenarios[i] * " 100 modes", color = colors[i])
    ylims!(ax, yrange...)
end
hideydecorations!(ax, grid = false)

ax = Axis(fig_hurs_tas[2,2]; title = "1000 Modes", xlabel = "Year")
for (i, field) in enumerate(fields)
    Ts = temperatures[i]
    error_list1000 = [norm((reshape(field[:, :, j], size(Φ)[1]) .- (eof_model1000[:, 1] .+ Ts[j] * eof_model1000[:, 2]) )[:] .* sqrt_f_metric) for j in eachindex(Ts)]
    lines!(ax, ts[i], error_list1000, label = scenarios[i] * " 1000 modes", color = colors[i])
    ylims!(ax, yrange...)
end
hideydecorations!(ax, grid = false)

field_name = "tas" # ["tas", "hurs"]

hfile = h5open(save_directory * field_name * "_basis.hdf5", "r")
latitude = read(hfile["latitude"])
longitude = read(hfile["longitude"])
metric = read(hfile["metric"])
close(hfile)
sqrt_f_metric = sqrt.(reshape(metric, 192 * 96))

hfile = h5open(save_directory * field_name * "_mean_regression.hdf5", "r")
regression_coefficients = read(hfile["regression_coefficients 1"])
linear_coefficients = zeros(Float32, size(regression_coefficients)..., 12)
for month in ProgressBar(1:12)
    linear_coefficients[:, :, month] .= read(hfile["regression_coefficients $month"])
end
close(hfile)

Φ = eof_basis(field_name)

averaged_coefficients = mean(linear_coefficients, dims = 3)[:, :, 1]

eof_model10 = zeros(Float32, size(Φ)[1]..., 2)
eof_model100 = zeros(Float32, size(Φ)[1]..., 2)
eof_model1000 = zeros(Float32, size(Φ)[1]..., 2)

for i in ProgressBar(1:10)
    eof_model10[:, 1] .+= Φ[:, i] * averaged_coefficients[i, 1]
    eof_model10[:, 2] .+= Φ[:, i] * averaged_coefficients[i, 2]
end
for i in ProgressBar(1:100)
    eof_model100[:, 1] .+= Φ[:, i] * averaged_coefficients[i, 1]
    eof_model100[:, 2] .+= Φ[:, i] * averaged_coefficients[i, 2]
end
for i in ProgressBar(1:1000)
    eof_model1000[:, 1] .+= Φ[:, i] * averaged_coefficients[i, 1]
    eof_model1000[:, 2] .+= Φ[:, i] * averaged_coefficients[i, 2]
end
##
filename  = field_name * "_pattern_scaling.hdf5"
hfile = h5open(save_directory * filename, "r")
linear_fit_yearly = read(hfile["linear fit yearly"])
close(hfile)
scenarios = ["historical", "ssp585", "ssp119", "ssp245"]
temperatures = []
fields = []
for scenario in scenarios
    temperature = regression_variable(scenario)
    a, b = ensemble_averaging(scenario, field_name; ensemble_members = 29)
    push!(temperatures, temperature)
    push!(fields, b)
end

##
i = 2
Ts = temperatures[i]
field = fields[i]
error_list10 = [norm((reshape(field[:, :, j], size(Φ)[1]) .- (eof_model10[:, 1] .+ Ts[j] * eof_model10[:, 2]))[:] .* sqrt_f_metric) for j in eachindex(Ts)]
ymax = maximum(error_list10)

yrange = (0, ymax)
ts = [1850:2014, 2015:2100, 2015:2100, 2015:2100]
ax = Axis(fig_hurs_tas[1,1]; title = "Pattern Scaling", xlabel = "Year", ylabel = "Temperature RMSE (K)")
for (i, field) in enumerate(fields)
    Ts = temperatures[i]
    error_list = [norm((field[:, :, j] .- (linear_fit_yearly[:, :, 1] .+ Ts[j] * linear_fit_yearly[:, :, 2]))[:] .* sqrt_f_metric) for j in eachindex(Ts)]
    lines!(ax, ts[i], error_list, label = scenarios[i], color = colors[i])
    ylims!(ax, yrange...)
end
axislegend(ax; position = :lt)

ax = Axis(fig_hurs_tas[1,4]; title = "10 Modes", xlabel = "Year")
for (i, field) in enumerate(fields)
    Ts = temperatures[i]
    error_list10 = [norm((reshape(field[:, :, j], size(Φ)[1]) .- (eof_model10[:, 1] .+ Ts[j] * eof_model10[:, 2]) )[:] .* sqrt_f_metric) for j in eachindex(Ts)]
    lines!(ax, ts[i], error_list10, label = scenarios[i] * " 10 modes", color = colors[i])
    ylims!(ax, yrange...)
end
hideydecorations!(ax, grid = false)

ax = Axis(fig_hurs_tas[1,3]; title = "100 Modes", xlabel = "Year")
for (i, field) in enumerate(fields)
    Ts = temperatures[i]
    error_list100 = [norm((reshape(field[:, :, j], size(Φ)[1]) .- (eof_model100[:, 1] .+ Ts[j] * eof_model100[:, 2]) )[:] .* sqrt_f_metric) for j in eachindex(Ts)]
    lines!(ax, ts[i], error_list100, label = scenarios[i] * " 100 modes", color = colors[i])
    ylims!(ax, yrange...)
end
hideydecorations!(ax, grid = false)

ax = Axis(fig_hurs_tas[1,2]; title = "1000 Modes", xlabel = "Year")
for (i, field) in enumerate(fields)
    Ts = temperatures[i]
    error_list1000 = [norm((reshape(field[:, :, j], size(Φ)[1]) .- (eof_model1000[:, 1] .+ Ts[j] * eof_model1000[:, 2]) )[:] .* sqrt_f_metric) for j in eachindex(Ts)]
    lines!(ax, ts[i], error_list1000, label = scenarios[i] * " 1000 modes", color = colors[i])
    ylims!(ax, yrange...)
end
hideydecorations!(ax, grid = false)

save("hurs_tas_pattern_scaling_and_model_errors_3_metric.png", fig_hurs_tas)