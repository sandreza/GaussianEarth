#### start by essentially repeating emulator.jl from paper_plots
## !!! will need to standardize how guassian_emulator is saved to separately save different orders of fit # probably just add a \mumodel_2 variable

field = "tas"

# create quadratic emulator
rC = []
hfile = h5open("./figures/" * field * "_mean_regression_quadratic.hdf5", "r") 
for month in 1:12
    regression_coefficients = read(hfile["regression_coefficients $month"])
    push!(rC, regression_coefficients)
end
close(hfile)

hfile = h5open(save_directory * field * "_basis.hdf5", "r")
basis = read(hfile["basis"] )
close(hfile)

covsave = h5open(save_directory * field * "_covariances_model.hdf5", "r")
Ls = []
for month in 1:12
    L = read(covsave["L" * "$month"])
    push!(Ls, L)
end
close(covsave)

hfile = h5open("./figures/" * field * "_gaussian_model_quadratic.hdf5", "w") # CHANGE THIS to be the save_directory
Lmodel = zeros(Float32, size(Ls[1])..., 12)
[Lmodel[:, :, :, i] .= Ls[i] for i in 1:12]
μmodel = zeros(Float32, size(rC[1])..., 12)
[μmodel[:, :, i] .= rC[i] for i in 1:12]
hfile["mean"] = μmodel
hfile["L model"] = Lmodel
hfile["basis"] = basis
hfile["scale factor"] = 273 
close(hfile)

##### 
# load in saved out quadratic emulator
hfile = h5open("./figures/" * field * "_gaussian_model_quadratic.hdf5", "r")
μmodel = read(hfile["mean"])
Lmodel = read(hfile["L model"])
basis = read(hfile["basis"])
close(hfile)
quadratic_emulator = GaussianEmulator(μmodel, Lmodel, basis)

#load in saved out linear emulator
hfile = h5open(save_directory * field * "_gaussian_model.hdf5", "r")
μmodel = read(hfile["mean"])
Lmodel = read(hfile["L model"])
basis = read(hfile["basis"])
close(hfile)
emulator = GaussianEmulator(μmodel, Lmodel, basis)


# now we have emulator and quadratic emulator, can do the error comparison
include("utils.jl")
##
save_directory = "/net/fs06/d3/sandre/GaussianEarthData/"
data_directory = "/net/fs06/d3/mgeo/CMIP6/interim/"
scenario_directories = readdir(data_directory)
current_path = joinpath(data_directory, scenario_directories[1])
variable_directories = readdir(current_path)
##
field_name = "tas" 
colors = [:red4, :red, :indigo, :magenta3]

# get metric correction terms
if process_data
    hfile = h5open(save_directory * field_name * "_basis.hdf5", "r")
    latitude = read(hfile["latitude"])
    longitude = read(hfile["longitude"])
    metric = read(hfile["metric"])
    close(hfile)
    sqrt_f_metric = sqrt.(reshape(metric, 192 * 96))
end

# get eof models for 10, 100, and 1000 modes 
if process_data
    hfile = h5open(save_directory * field_name * "_mean_regression.hdf5", "r")
    regression_coefficients = read(hfile["regression_coefficients 1"])
    linear_coefficients = zeros(Float32, size(regression_coefficients)..., 12)
    for month in ProgressBar(1:12)
        linear_coefficients[:, :, month] .= read(hfile["regression_coefficients $month"])
    end
    close(hfile)

    Φ = eof_basis(field_name) #this is the only thing we're using here
    averaged_coefficients = mean(linear_coefficients, dims = 3)[:, :, 1] #he was using the annual-avg coefficients

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
end

# repeat for quadratic
if process_data
    hfile = h5open("./figures/" * field_name * "_mean_regression_quadratic.hdf5", "r")
    regression_coefficients = read(hfile["regression_coefficients 1"])
    quadratic_coefficients = zeros(Float32, size(regression_coefficients)..., 12)
    for month in ProgressBar(1:12)
        quadratic_coefficients[:, :, month] .= read(hfile["regression_coefficients $month"])
    end
    close(hfile)

    # Φ = eof_basis(field_name) #this is the only thing we're using here
    averaged_quad_coefficients = mean(quadratic_coefficients, dims = 3)[:, :, 1] #he was using the annual-avg coefficients

    quad_eof_model10 = zeros(Float32, size(Φ)[1]..., 3)
    quad_eof_model100 = zeros(Float32, size(Φ)[1]..., 3)
    quad_eof_model1000 = zeros(Float32, size(Φ)[1]..., 3)

    for i in ProgressBar(1:10)
        quad_eof_model10[:, 1] .+= Φ[:, i] * averaged_quad_coefficients[i, 1]
        quad_eof_model10[:, 2] .+= Φ[:, i] * averaged_quad_coefficients[i, 2]
        quad_eof_model10[:, 3] .+= Φ[:, i] * averaged_quad_coefficients[i, 3]
    end
    for i in ProgressBar(1:100)
        quad_eof_model100[:, 1] .+= Φ[:, i] * averaged_quad_coefficients[i, 1]
        quad_eof_model100[:, 2] .+= Φ[:, i] * averaged_quad_coefficients[i, 2]
        quad_eof_model100[:, 3] .+= Φ[:, i] * averaged_quad_coefficients[i, 3]
    end
    for i in ProgressBar(1:1000)
        quad_eof_model1000[:, 1] .+= Φ[:, i] * averaged_quad_coefficients[i, 1]
        quad_eof_model1000[:, 2] .+= Φ[:, i] * averaged_quad_coefficients[i, 2]
        quad_eof_model1000[:, 3] .+= Φ[:, i] * averaged_quad_coefficients[i, 3]
    end
end

## get true data 
if process_data
    scenarios = ["historical", "ssp585", "ssp119", "ssp245"]
    temperatures = [] 
    fields = [] #list of four arrays of time x space data values for given variable
    for scenario in scenarios
        temperature = regression_variable(scenario) #this gets the list of temps to regress onto
        a, b = ensemble_averaging(scenario, field_name; ensemble_members = 29, return_std=false)
        push!(temperatures, temperature)
        push!(fields, b)
    end
end


## make the plot
ts = 40
xls = 40 
yls = 40
tls = 40
legend_ls = 35
resolution = (1350, 600)
common_options = (; titlesize = ts, xlabelsize = xls, ylabelsize = yls, xticklabelsize = tls, yticklabelsize = tls)
linewidth = 5

scenario_names = ["Historical", "SSP5-8.5", "SSP1-1.9", "SSP2-4.5"]
i = 2
Ts = temperatures[i]
field = fields[i]
error_list10 = [norm((reshape(field[:, :, j], size(Φ)[1]) .- (eof_model10[:, 1] .+ Ts[j] * eof_model10[:, 2]))[:] .* sqrt_f_metric) for j in eachindex(Ts)]
quad_error_list10 = [norm((reshape(field[:, :, j], size(Φ)[1]) .- (quad_eof_model10[:, 1] .+ Ts[j] * quad_eof_model10[:, 2] .+ Ts[j]^2 * quad_eof_model10[:, 3]))[:] .* sqrt_f_metric) for j in eachindex(Ts)]
diff = error_list10 .- quad_error_list10
ymax = maximum(error_list10)
maxdiff = maximum(diff) 

xrange = (1850, 2140)
# yrange = (0, ymax)
# yrange = (-0.5, 0.5)
yrange = (-(maxdiff + 0.1), maxdiff + 0.1)
fig_tas = Figure(; resolution) 
ts = [1850:2014, 2015:2100, 2015:2100, 2015:2100]

ax = Axis(fig_tas[1,3]; title = "10 Modes", xlabel = "Year", common_options...)
for (i, field) in enumerate(fields)
    Ts = temperatures[i]
    error_list10 = [norm((reshape(field[:, :, j], size(Φ)[1]) .- (eof_model10[:, 1] .+ Ts[j] * eof_model10[:, 2]) )[:] .* sqrt_f_metric) for j in eachindex(Ts)]
    quad_error_list10 = [norm((reshape(field[:, :, j], size(Φ)[1]) .- (quad_eof_model10[:, 1] .+ Ts[j] * quad_eof_model10[:, 2] .+ Ts[j]^2 * quad_eof_model10[:, 3]))[:] .* sqrt_f_metric) for j in eachindex(Ts)]
    diff = error_list10 .- quad_error_list10
    # lines!(ax, ts[i], error_list10, label = scenario_names[i], color = colors[i], linestyle = :dash, linewidth = linewidth)
    # lines!(ax, ts[i], quad_error_list10,  color = colors[i], linestyle = :dot, linewidth = linewidth)
    lines!(ax, ts[i], diff, color = colors[i], label = scenario_names[i], linestyle = :solid, linewidth = linewidth)
    xlims!(ax, xrange...)
    ylims!(ax, yrange...)
end
hideydecorations!(ax, grid = false)
axislegend(ax; position = :lt, labelsize = legend_ls)

ax = Axis(fig_tas[1,2]; title = "100 Modes", xlabel = "Year", common_options...)
for (i, field) in enumerate(fields)
    Ts = temperatures[i]
    error_list100 = [norm((reshape(field[:, :, j], size(Φ)[1]) .- (eof_model100[:, 1] .+ Ts[j] * eof_model100[:, 2]) )[:] .* sqrt_f_metric) for j in eachindex(Ts)]
    quad_error_list100 = [norm((reshape(field[:, :, j], size(Φ)[1]) .- (quad_eof_model100[:, 1] .+ Ts[j] * quad_eof_model100[:, 2] .+ Ts[j]^2 * quad_eof_model100[:, 3]))[:] .* sqrt_f_metric) for j in eachindex(Ts)]
    diff = error_list100 .- quad_error_list100
    # lines!(ax, ts[i], error_list100, label = scenario_names[i], color = colors[i], linestyle = :dash, linewidth = linewidth)
    # lines!(ax, ts[i], quad_error_list100, color = colors[i], linestyle = :dot, linewidth = linewidth)
    lines!(ax, ts[i], diff, color = colors[i], linestyle = :solid, linewidth = linewidth)
    xlims!(ax, xrange...)
    ylims!(ax, yrange...)
end
hideydecorations!(ax, grid = false)

ax = Axis(fig_tas[1,1]; title = "1000 Modes", xlabel = "Year", common_options...)
for (i, field) in enumerate(fields)
    Ts = temperatures[i]
    error_list1000 = [norm((reshape(field[:, :, j], size(Φ)[1]) .- (eof_model1000[:, 1] .+ Ts[j] * eof_model1000[:, 2]) )[:] .* sqrt_f_metric) for j in eachindex(Ts)]
    quad_error_list1000 = [norm((reshape(field[:, :, j], size(Φ)[1]) .- (quad_eof_model1000[:, 1] .+ Ts[j] * quad_eof_model1000[:, 2] .+ Ts[j]^2 * quad_eof_model1000[:, 3]))[:] .* sqrt_f_metric) for j in eachindex(Ts)]
    diff = error_list1000 .- quad_error_list1000
    # lines!(ax, ts[i], error_list1000, label = scenario_names[i], color = colors[i], linestyle = :dash, linewidth = linewidth)
    # lines!(ax, ts[i], quad_error_list1000,  color = colors[i], linestyle = :dot, linewidth = linewidth)
    lines!(ax, ts[i], diff, color = colors[i], linestyle = :solid, linewidth = linewidth)
    xlims!(ax, xrange...)
    ylims!(ax, yrange...)
end
# display(fig_tas)

save(figure_directory * field_name * "_errors_diff.png", fig_tas)