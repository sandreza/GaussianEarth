ts = 60
xls = 70 
yls = 70
tls = 70
legend_ls = 60
resolution = (400, 150) .* 12
lw = 7
global_common_options = (; titlesize = ts, xlabelsize = xls, ylabelsize = yls, xticklabelsize = tls, yticklabelsize = tls)
##
hfile = h5open("new_scenarios.hdf5", "r")
scenario_1 = read(hfile["scenario_1"])
scenario_2 = read(hfile["scenario_2"])
close(hfile)

##
field = "tas" 
field_name = "tas"

#load in saved out linear emulator
# include("emulator.jl")
# include("emulator_hurs.jl")
d = 1000

# load in arbitrary scenario
min_temp_1 = scenario_1[56]/273 #2030
max_temp_1 = scenario_1[end]/273 #2070

min_temp_2 = scenario_2[56]/273 #2030
max_temp_2 = scenario_2[end]/273 #2070

hfile = h5open(save_directory * field * "_basis.hdf5", "r")
latitude = read(hfile["latitude"])
longitude = read(hfile["longitude"])
metric = read(hfile["metric"])
close(hfile)
fmetric = reshape(metric, (192*96, 1))

index_1 = [157, 46]
index_2 = [54, 48]
index_3 = [19, 87]

function losangeles(field)
    rfield = reshape(field, (192, 96))
    return rfield[130, 67]
end

function sahara(field)
    rfield = reshape(field, (192, 96))
    return rfield[11, 58] #northern chad
end

function norcal(field)
    rfield = reshape(field, (192, 96))
    return rfield[128, 70]
end

function manaus(field)
    rfield = reshape(field, (192, 96))
    return rfield[161, 46]
end

function chicago(field)
    rfield = reshape(field, (192, 96))
    return rfield[146, 70]
end

# function nw_india(field)
#     rfield = reshape(field, (192, 96))
#     return rfield[148, 61]
# end

function nw_india(field) #CORRECT ONE
    rfield = reshape(field, (192, 96))
    return rfield[46, 61]
end
    

observables = [nw_india, losangeles]

##
fig = Figure(; resolution)
observable_names_base = ["Observable 1 ", "Observable 2 "]
observable_names = observable_names_base .* ["(Jul, 2100)"]
month_string = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
xticks = (1:12, month_string)
num_stds = 2
for i in eachindex(observables)
    month = 7
    emulator.month[1] = month
    emulator.global_mean_temperature[1] = max_temp_1
    emulator_hurs.month[1] = month
    emulator_hurs.global_mean_temperature[1] = max_temp_1
    tasmean, tasstd = emulator_mean_variance_linear_functionals(observables, emulator)
    hursmean, hursstd = emulator_mean_variance_linear_functionals(observables, emulator_hurs)

    common_options_1 = (; xlabel = "Temperature (K)")
    common_options_2 = (; xlabel = "Relative Humidity (%)")

    ax1 = Axis(fig[1, 2*i]; title = observable_names[i], common_options_1..., global_common_options..., ylabel = "PDF")
    xs, ys = gaussian_grid(tasmean[i], tasstd[i])
    lines!(ax1, xs, ys, color = :blue, linewidth = lw, label = "Scenario 1")

    ax2 = Axis(fig[2, 2*i]; title = observable_names[i], common_options_2..., global_common_options..., ylabel = "PDF")
    xs, ys = gaussian_grid(hursmean[i], hursstd[i])
    lines!(ax2, xs, ys, color = :blue, linewidth = lw)

    emulator.month[1] = month
    emulator.global_mean_temperature[1] = max_temp_2
    emulator_hurs.month[1] = month
    emulator_hurs.global_mean_temperature[1] = max_temp_2
    tasmean, tasstd = emulator_mean_variance_linear_functionals(observables, emulator)
    hursmean, hursstd = emulator_mean_variance_linear_functionals(observables, emulator_hurs)

    xs, ys = gaussian_grid(tasmean[i], tasstd[i])
    lines!(ax1, xs, ys, color = :orangered, linewidth = lw, label = "Scenario 2")
    xs, ys = gaussian_grid(hursmean[i], hursstd[i])
    lines!(ax2, xs, ys, color = :orangered, linewidth = lw)

    if i == 1
        axislegend(ax1, position = :lt, labelsize = legend_ls)
    end

    ax3 = Axis(fig[1, 2*i-1]; title = observable_names_base[i] * "(2100)", xticks, common_options_1..., global_common_options..., ylabel = "Temperature (K)")
    ax4 = Axis(fig[2, 2*i-1]; title = observable_names_base[i] * "(2100)", xticks, common_options_2..., global_common_options..., ylabel = "Relative Humidity (%)")
    tasmeanlist = zeros(Float32, 12)
    tasstdlist = zeros(Float32, 12)
    hursmeanlist = zeros(Float32, 12)
    hursstdlist = zeros(Float32, 12)
    for month in 1:12
        emulator.month[1] = month
        emulator.global_mean_temperature[1] = max_temp_1
        tasmean, tasstd = emulator_mean_variance_linear_functionals([observables[i]], emulator)
        emulator_hurs.month[1] = month
        emulator_hurs.global_mean_temperature[1] = max_temp_1
        hursmean, hursstd = emulator_mean_variance_linear_functionals([observables[i]], emulator_hurs)
        
        tasmeanlist[month] = tasmean[1]
        hursmeanlist[month] = hursmean[1]
        tasstdlist[month] = tasstd[1]
        hursstdlist[month] = hursstd[1]
    end
    lines!(ax3, 1:12, tasmeanlist, color = :blue, linewidth = lw)
    band!(ax3, 1:12, tasmeanlist .- num_stds .* tasstdlist, tasmeanlist .+ num_stds .* tasstdlist, color = (:blue, 0.2))
    lines!(ax4, 1:12, hursmeanlist, color = :blue, linewidth = lw)
    band!(ax4, 1:12, hursmeanlist .- num_stds .* hursstdlist, hursmeanlist .+ num_stds .* hursstdlist, color = (:blue, 0.2))

    
    tasmeanlist = zeros(Float32, 12)
    tasstdlist = zeros(Float32, 12)
    hursmeanlist = zeros(Float32, 12)
    hursstdlist = zeros(Float32, 12)
    for month in 1:12
        emulator.month[1] = month
        emulator.global_mean_temperature[1] = max_temp_2
        tasmean, tasstd = emulator_mean_variance_linear_functionals([observables[i]], emulator)
        emulator_hurs.month[1] = month
        emulator_hurs.global_mean_temperature[1] = max_temp_2
        hursmean, hursstd = emulator_mean_variance_linear_functionals([observables[i]], emulator_hurs)
        
        tasmeanlist[month] = tasmean[1]
        hursmeanlist[month] = hursmean[1]
        tasstdlist[month] = tasstd[1]
        hursstdlist[month] = hursstd[1]
    end
    lines!(ax3, 1:12, tasmeanlist, color = :orangered, linewidth = lw)
    band!(ax3, 1:12, tasmeanlist .- num_stds .* tasstdlist, tasmeanlist .+ num_stds .* tasstdlist, color = (:orangered, 0.2))
    lines!(ax4, 1:12, hursmeanlist, color = :orangered, linewidth = lw)
    band!(ax4, 1:12, hursmeanlist .- num_stds .* hursstdlist, hursmeanlist .+ num_stds .* hursstdlist, color = (:orangered, 0.2))

end
display(fig)
# save(figure_directory * "case_study.png", fig)