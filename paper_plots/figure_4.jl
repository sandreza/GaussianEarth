ts = 40
xls = 40 
yls = 40
tls = 40
legend_ls = 35
resolution = (1800, 600)
common_options = (; titlesize = ts, xlabelsize = xls, ylabelsize = yls, xticklabelsize = tls, yticklabelsize = tls)
linewidth = 5

##
field_name = "tas" 
colors = [:red4, :red, :indigo, :magenta3]
if process_data
    hfile = h5open(save_directory * field_name * "_basis.hdf5", "r")
    latitude = read(hfile["latitude"])
    longitude = read(hfile["longitude"])
    metric = read(hfile["metric"])
    close(hfile)
    sqrt_f_metric = sqrt.(reshape(metric, 192 * 96))
end

if process_data
    # hfile = h5open(save_directory * field_name * "_mean_regression.hdf5", "r")
    # regression_coefficients = read(hfile["regression_coefficients 1"])
    # linear_coefficients = zeros(Float32, size(regression_coefficients)..., 12)
    # for month in ProgressBar(1:12)
    #     linear_coefficients[:, :, month] .= read(hfile["regression_coefficients $month"])
    # end
    # close(hfile)
    hfile = h5open(save_directory * field * "_mean_regression_quadratic.hdf5", "r")
    regression_coefficients = read(hfile["regression_coefficients 1"])
    quadratic_coefficients = zeros(Float32, size(regression_coefficients)..., 12)
    for month in ProgressBar(1:12)
        quadratic_coefficients[:, :, month] .= read(hfile["regression_coefficients $month"])
    end
    close(hfile)

    Φ = eof_basis(field_name)
    # averaged_coefficients = mean(linear_coefficients, dims = 3)[:, :, 1]
    averaged_coefficients = mean(quadratic_coefficients, dims = 3)[:, :, 1]

    eof_model10 = zeros(Float32, size(Φ)[1]..., 3)
    eof_model100 = zeros(Float32, size(Φ)[1]..., 3)
    eof_model1000 = zeros(Float32, size(Φ)[1]..., 3)

    for i in ProgressBar(1:10)
        eof_model10[:, 1] .+= Φ[:, i] * averaged_coefficients[i, 1]
        eof_model10[:, 2] .+= Φ[:, i] * averaged_coefficients[i, 2]
        eof_model10[:, 3] .+= Φ[:, i] * averaged_coefficients[i, 3]
    end
    for i in ProgressBar(1:100)
        eof_model100[:, 1] .+= Φ[:, i] * averaged_coefficients[i, 1]
        eof_model100[:, 2] .+= Φ[:, i] * averaged_coefficients[i, 2]
        eof_model100[:, 3] .+= Φ[:, i] * averaged_coefficients[i, 3]
    end
    for i in ProgressBar(1:1000)
        eof_model1000[:, 1] .+= Φ[:, i] * averaged_coefficients[i, 1]
        eof_model1000[:, 2] .+= Φ[:, i] * averaged_coefficients[i, 2]
        eof_model1000[:, 3] .+= Φ[:, i] * averaged_coefficients[i, 3]
    end
end

##
if process_data
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
end
##
scenario_names = ["Historical", "SSP5-8.5", "SSP1-1.9", "SSP2-4.5"]
i = 2
Ts = temperatures[i]
field = fields[i]
# error_list10 = [norm((reshape(field[:, :, j], size(Φ)[1]) .- (eof_model10[:, 1] .+ Ts[j] * eof_model10[:, 2]))[:] .* sqrt_f_metric) for j in eachindex(Ts)]
error_list10 = [norm((reshape(field[:, :, j], size(Φ)[1]) .- (eof_model10[:, 1] .+ Ts[j] * eof_model10[:, 2] .+ Ts[j]^2 * eof_model10[:, 3]))[:] .* sqrt_f_metric) for j in eachindex(Ts)]
ymax = maximum(error_list10)

xrange = (1850, 2140)
yrange = (0, ymax)
fig_tas = Figure(; resolution) 
ts = [1850:2014, 2015:2100, 2015:2100, 2015:2100]
ax = Axis(fig_tas[1,1]; title = "Pattern Scaling", xlabel = "Year", ylabel = "Temperature RMSE (K)", common_options...)
for (i, field) in enumerate(fields)
    Ts = temperatures[i]
    error_list = [norm((field[:, :, j] .- (linear_fit_yearly[:, :, 1] .+ Ts[j] * linear_fit_yearly[:, :, 2]))[:] .* sqrt_f_metric) for j in eachindex(Ts)]
    lines!(ax, ts[i], error_list, label = scenario_names[i], color = colors[i], linewidth = linewidth)
    xlims!(ax, xrange...)
    ylims!(ax, yrange...)
end
axislegend(ax; position = :lt, labelsize = legend_ls)

ax = Axis(fig_tas[1,4]; title = "10 Modes", xlabel = "Year", common_options...)
for (i, field) in enumerate(fields)
    Ts = temperatures[i]
    # error_list10 = [norm((reshape(field[:, :, j], size(Φ)[1]) .- (eof_model10[:, 1] .+ Ts[j] * eof_model10[:, 2]) )[:] .* sqrt_f_metric) for j in eachindex(Ts)]
    error_list10 = [norm((reshape(field[:, :, j], size(Φ)[1]) .- (eof_model10[:, 1] .+ Ts[j] * eof_model10[:, 2] .+ Ts[j]^2 * eof_model10[:, 3]))[:] .* sqrt_f_metric) for j in eachindex(Ts)]
    lines!(ax, ts[i], error_list10, label = scenario_names[i] * " 10 modes", color = colors[i], linewidth = linewidth)
    xlims!(ax, xrange...)
    ylims!(ax, yrange...)
end
hideydecorations!(ax, grid = false)

ax = Axis(fig_tas[1,3]; title = "100 Modes", xlabel = "Year", common_options...)
for (i, field) in enumerate(fields)
    Ts = temperatures[i]
    # error_list100 = [norm((reshape(field[:, :, j], size(Φ)[1]) .- (eof_model100[:, 1] .+ Ts[j] * eof_model100[:, 2]) )[:] .* sqrt_f_metric) for j in eachindex(Ts)]
    error_list100 = [norm((reshape(field[:, :, j], size(Φ)[1]) .- (eof_model100[:, 1] .+ Ts[j] * eof_model100[:, 2] .+ Ts[j]^2 * eof_model100[:, 3]))[:] .* sqrt_f_metric) for j in eachindex(Ts)]
    lines!(ax, ts[i], error_list100, label = scenario_names[i] * " 100 modes", color = colors[i], linewidth = linewidth)
    xlims!(ax, xrange...)
    ylims!(ax, yrange...)
end
hideydecorations!(ax, grid = false)

ax = Axis(fig_tas[1,2]; title = "1000 Modes", xlabel = "Year", common_options...)
for (i, field) in enumerate(fields)
    Ts = temperatures[i]
    # error_list1000 = [norm((reshape(field[:, :, j], size(Φ)[1]) .- (eof_model1000[:, 1] .+ Ts[j] * eof_model1000[:, 2]) )[:] .* sqrt_f_metric) for j in eachindex(Ts)]
    error_list1000 = [norm((reshape(field[:, :, j], size(Φ)[1]) .- (eof_model1000[:, 1] .+ Ts[j] * eof_model1000[:, 2] .+ Ts[j]^2 * eof_model1000[:, 3]))[:] .* sqrt_f_metric) for j in eachindex(Ts)]
    lines!(ax, ts[i], error_list1000, label = scenario_names[i] * " 1000 modes", color = colors[i], linewidth = linewidth)
    xlims!(ax, xrange...)
    ylims!(ax, yrange...)
end
hideydecorations!(ax, grid = false)
save(figure_directory * field_name * "_pattern_scaling_and_model_errors_3_metric.png", fig_tas)