using NCDatasets, LinearAlgebra, Statistics, HDF5, ProgressBars
include("utils.jl")
##
save_directory = "/net/fs06/d3/sandre/GaussianEarthData/"
data_directory = "/net/fs06/d3/mgeo/CMIP6/interim/"
scenario_directories = readdir(data_directory)
current_path = joinpath(data_directory, scenario_directories[1])
variable_directories = readdir(current_path)
##
field = "tas"
hfile = h5open(save_directory * field * "_mean_regression.hdf5", "r")
regression_coefficients = read(hfile["regression_coefficients 1"])
linear_coefficients = zeros(Float32, size(regression_coefficients)..., 12)
for month in ProgressBar(1:12)
    linear_coefficients[:, :, month] .= read(hfile["regression_coefficients $month"])
end
close(hfile)

Φ = eof_basis("tas")

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
include("pattern_scaling_error.jl")
##

fig = Figure() 
ax = Axis(fig[1,1])
ts = [1850:2014, 2015:2100, 2015:2100, 2015:2100]
for (i, field) in enumerate(fields)
    Ts = temperatures[i]
    error_list = [norm((field[:, :, j] .- (linear_fit_yearly[:, :, 1] .+ Ts[j] * linear_fit_yearly[:, :, 2]))[:]) for j in eachindex(Ts)]
    error_list10 = [norm(reshape(field[:, :, j], size(Φ)[1]) .- (eof_model10[:, 1] .+ Ts[j] * eof_model10[:, 2])) for j in eachindex(Ts)]
    error_list100 = [norm(reshape(field[:, :, j], size(Φ)[1]) .- (eof_model100[:, 1] .+ Ts[j] * eof_model100[:, 2])) for j in eachindex(Ts)]
    error_list1000 = [norm(reshape(field[:, :, j], size(Φ)[1]) .- (eof_model1000[:, 1] .+ Ts[j] * eof_model1000[:, 2])) for j in eachindex(Ts)]
    lines!(ax, ts[i], error_list, label = scenarios[i])
    lines!(ax, ts[i], error_list10, label = scenarios[i] * " 10 modes")
    lines!(ax, ts[i], error_list100, label = scenarios[i] * " 100 modes")
    lines!(ax, ts[i], error_list1000, label = scenarios[i] * " 1000 modes")
end
axislegend(ax; position = :lt)
save("pattern_scaling_and_model_errors.png", fig)

##
yrange = (0, 150)
fig = Figure(resolution = (2000, 400)) 
ts = [1850:2014, 2015:2100, 2015:2100, 2015:2100]
ax = Axis(fig[1,1]; title = "Pattern Scaling", xlabel = "Year", ylabel = "Error")
for (i, field) in enumerate(fields)
    Ts = temperatures[i]
    error_list = [norm((field[:, :, j] .- (linear_fit_yearly[:, :, 1] .+ Ts[j] * linear_fit_yearly[:, :, 2]))[:]) for j in eachindex(Ts)]
    lines!(ax, ts[i], error_list, label = scenarios[i])
    ylims!(ax, yrange...)
end
axislegend(ax; position = :lt)

ax = Axis(fig[1,4]; title = "10 Modes", xlabel = "Year")
for (i, field) in enumerate(fields)
    Ts = temperatures[i]
    error_list10 = [norm(reshape(field[:, :, j], size(Φ)[1]) .- (eof_model10[:, 1] .+ Ts[j] * eof_model10[:, 2])) for j in eachindex(Ts)]
    lines!(ax, ts[i], error_list10, label = scenarios[i] * " 10 modes")
    ylims!(ax, yrange...)
end
hideydecorations!(ax, grid = false)

ax = Axis(fig[1,3]; title = "100 Modes", xlabel = "Year")
for (i, field) in enumerate(fields)
    Ts = temperatures[i]
    error_list100 = [norm(reshape(field[:, :, j], size(Φ)[1]) .- (eof_model100[:, 1] .+ Ts[j] * eof_model100[:, 2])) for j in eachindex(Ts)]
    lines!(ax, ts[i], error_list100, label = scenarios[i] * " 100 modes")
    ylims!(ax, yrange...)
end
hideydecorations!(ax, grid = false)

ax = Axis(fig[1,2]; title = "1000 Modes", xlabel = "Year")
for (i, field) in enumerate(fields)
    Ts = temperatures[i]
    error_list1000 = [norm(reshape(field[:, :, j], size(Φ)[1]) .- (eof_model1000[:, 1] .+ Ts[j] * eof_model1000[:, 2])) for j in eachindex(Ts)]
    lines!(ax, ts[i], error_list1000, label = scenarios[i] * " 1000 modes")
    ylims!(ax, yrange...)
end
hideydecorations!(ax, grid = false)
save("pattern_scaling_and_model_errors_3.png", fig)