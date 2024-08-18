using NCDatasets, LinearAlgebra, Statistics, HDF5, ProgressBars, CairoMakie

save_directory = "/net/fs06/d3/sandre/GaussianEarthData/"

include("utils.jl")

function ensemble_averaging(scenario, variable; data_directory = "/net/fs06/d3/mgeo/CMIP6/interim/", ensemble_members = 30)
    current_path = joinpath(data_directory, scenario)
    local_current_path = joinpath(current_path, variable)
    file_names = readdir(local_current_path)
    fields = Array{Float32, 3}[]
    if length(file_names) > 0 # sometimes the directory is empty
        file_name = file_names[1] # pick the first file for obtaining varaibles
        file_path = joinpath(local_current_path, file_name)
        for (i, file_name) in ProgressBar(enumerate(file_names[1:ensemble_members]))
            file_path = joinpath(local_current_path, file_name)
            ds = Dataset(file_path)
            field = Float32.(ds[variable][:,:,:])
            push!(fields, field)
        end
    end
    time1 = size(fields[1])[end]
    monthtime1 = time1 ÷ 12
    all_together = zeros(Float32, size(fields[1])[1], size(fields[1])[2], monthtime1, 12, ensemble_members);
    for ω in ProgressBar(1:ensemble_members)
        for month in 1:12
            @inbounds all_together[:,:, 1:monthtime1, month, ω] .= fields[ω][:, :, month:12:end]
        end
    end
    return mean(all_together, dims = 5)[:,:,:,:,1], mean(all_together, dims = (4, 5))[:,:,:,1, 1]
end
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
