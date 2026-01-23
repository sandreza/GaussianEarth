ts = 40
xls = 40 
yls = 40
tls = 40
legend_ls = 35
resolution = (1350, 600)
common_options = (; titlesize = ts, xlabelsize = xls, ylabelsize = yls, xticklabelsize = tls, yticklabelsize = tls)
linewidth = 5

##
field_name = "tas" 
colors = [:red4, :red, :indigo, :magenta3]

# get basis and metric correction terms
if process_data
    hfile = h5open(save_directory * field_name * "_basis.hdf5", "r")
    latitude = read(hfile["latitude"])
    longitude = read(hfile["longitude"])
    metric = read(hfile["metric"])
    close(hfile)
    sqrt_f_metric = sqrt.(reshape(metric, 192 * 96))

    Φ = eof_basis(field_name) 
end

## get true data for comparison
if process_data
    scenarios = ["historical", "ssp585", "ssp119", "ssp245"]
    temperatures = [] 
    fields = [] #list of four arrays of time x space data values for given variable #NEEDS to be standard deviation !!
    for scenario in scenarios
        temperature = regression_variable(scenario) #this gets the list of temps to regress onto
        a, b = ensemble_averaging(scenario, field_name; ensemble_members = 29, return_std=true) #gets the ensemble avg for that var? #changed this 
        push!(temperatures, temperature)
        push!(fields, a[:,:,:,:]) # stds for all months
    end
end

## get the errors

error_historical = zeros(Float32, 1, 3, 165, 12)
error_future = zeros(Float32, 3, 3, 86, 12)

for (e, d) in enumerate([10, 100, 1000])
    println("working on $d")
    basis = Φ[:, 1:d]
    for (i, field) in enumerate(fields)
        Ts = temperatures[i]
        for (j, t) in ProgressBar(enumerate(Ts))
            for month in 1:12
                Σ = emulator_variance(emulator; modes=d, month=month, global_mean_temperature=t)
                σ = sqrt.(diag((basis * Σ) * basis'))
                error = norm((reshape(field[:, :, j, month], size(Φ)[1]) .- σ)[:] .* sqrt_f_metric)
                if i == 1
                    error_historical[1, e, j, month] = error
                else
                    error_future[i-1, e, j, month] = error
                end
            end
        end
    end
end

hfile = h5open(save_directory * field_name * "_full_errors.hdf5", "w") 
write(hfile, "error_historical", error_historical)
write(hfile, "error_future", error_future)
close(hfile)

### load in pre-calculated error (that's saved out above) #EDIT TO REGENERATE IF NONEXISTENT
# hfile = h5open(save_directory * field_name * "_errors.hdf5", "r")
# error_historical = read(hfile["error_historical"])
# error_future = read(hfile["error_future"])
# close(hfile)

## average the errors
avg_error_historical = mean(error_historical, dims = 4)[:, :, :, 1]
avg_error_future = mean(error_future, dims = 4)[:, :, :, 1]

## set plotting parameters

scenario_names = ["Historical", "SSP5-8.5", "SSP1-1.9", "SSP2-4.5"]
xrange = (1850, 2140)
yrange = (0,1)
fig_tas = Figure(; resolution) 
ts = [1850:2014, 2015:2100, 2015:2100, 2015:2100]

ax = Axis(fig_tas[1,3]; title = "10 Modes", xlabel = "Year", common_options...)
for (i, field) in ProgressBar(enumerate(fields))
    if i == 1
        error_list10 = avg_error_historical[1, 1, :]
    else
        error_list10 = avg_error_future[i-1, 1, :]
    end
    lines!(ax, ts[i], error_list10, label = scenario_names[i] * " 10 modes", color = colors[i], linewidth = linewidth) 
    xlims!(ax, xrange...)
    ylims!(ax, yrange...)
end
hideydecorations!(ax, grid = false)

ax = Axis(fig_tas[1,2]; title = "100 Modes", xlabel = "Year", common_options...)
for (i, field) in ProgressBar(enumerate(fields))
    if i == 1
        error_list100 = avg_error_historical[1, 2, :]
    else
        error_list100 = avg_error_future[i-1, 2, :]
    end
    lines!(ax, ts[i], error_list100, label = scenario_names[i] * " 100 modes", color = colors[i], linewidth = linewidth)
    xlims!(ax, xrange...)
    ylims!(ax, yrange...)
end
hideydecorations!(ax, grid = false)

ax = Axis(fig_tas[1,1]; title = "1000 Modes", xlabel = "Year", common_options...)
for (i, field) in ProgressBar(enumerate(fields))
    if i == 1
        error_list1000 = avg_error_historical[1, 3, :]
    else
        error_list1000 = avg_error_future[i-1, 3, :]
    end
    lines!(ax, ts[i], error_list1000, label = scenario_names[i], color = colors[i], linewidth = linewidth)
    xlims!(ax, xrange...)
    ylims!(ax, yrange...)
end
axislegend(ax; position = :lt, labelsize = legend_ls)

# display(fig_tas)
save(figure_directory * field_name * "_model_errors_3_metric_std.png", fig_tas)
