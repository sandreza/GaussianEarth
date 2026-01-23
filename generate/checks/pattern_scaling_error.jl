using NCDatasets, LinearAlgebra, Statistics, HDF5, ProgressBars, CairoMakie

save_directory = "/net/fs06/d3/sandre/GaussianEarthData/"

include("utils.jl")
##
filename  = "tas_pattern_scaling.hdf5"
hfile = h5open(save_directory * filename, "r")
linear_fit_yearly = read(hfile["linear fit yearly"])
close(hfile)
##
scenarios = ["historical", "ssp585", "ssp119", "ssp245"]
temperatures = []
fields = []
for scenario in scenarios
    temperature = regression_variable(scenario)
    a, b = ensemble_averaging(scenario, "tas")
    push!(temperatures, temperature)
    push!(fields, b)
end
##
fig = Figure() 
ax = Axis(fig[1,1])
ts = [1850:2014, 2015:2100, 2015:2100, 2015:2100]
for (i, field) in enumerate(fields)
    Ts = temperatures[i]
    error_list = [norm(field[:, :, j] .- (linear_fit_yearly[:, :, 1] .+ Ts[j] * linear_fit_yearly[:, :, 2])) for j in eachindex(Ts)]
    lines!(ax, ts[i], error_list, label = scenarios[i])
end
axislegend(ax)
save("pattern_scaling_errors.png", fig)
