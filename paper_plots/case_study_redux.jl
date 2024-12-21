ts = 60
xls = 70 
yls = 70
tls = 70
legend_ls = 60
resolution = (400, 150) .* 12
lw = 7
global_common_options = (; titlesize = ts, xlabelsize = xls, ylabelsize = yls, xticklabelsize = tls, yticklabelsize = tls)
##
scenarios = ["historical", "ssp119", "ssp245", "ssp585"]
time_history = 1850:2014
time_future = 2015:2100
scenario_colors = Dict("historical" => :red4, "ssp585" => :red, "ssp245" => :magenta3, "ssp119" => :indigo)
scenario_labels = ["Historical","SSP1-1.9",  "SSP2-4.5", "SSP5-8.5"]

#
# arbitrary scenario
using DelimitedFiles
# structure: Year,Baseline,Current Scenario
csvdata = readdlm("./arbitrary-scenario-temp.csv", ',', Float64)
baseline = regression_variable("historical")[1]*273
arbitrary = csvdata[16:end, 3] .+ baseline

hfile = h5open("arbitrary_scenario.hdf5", "w")
write(hfile, "arbitrary", arbitrary)
close(hfile)

##
field = "tas" 
field_name = "tas"
# field = "hurs" 
# field_name = "hurs"

#load in saved out linear emulator
include("emulator.jl")
include("emulator_hurs.jl")
d = 1000

# load in arbitrary scenario
hfile = h5open("arbitrary_scenario.hdf5", "r")
temperatures = read(hfile["arbitrary"]) 
close(hfile)
min_temp = temperatures[1]/273 #2015
max_temp = temperatures[end]/273 #2100

temps585 = regression_variable("ssp585")
max_temp585 = temps585[end]

hfile = h5open(save_directory * field * "_basis.hdf5", "r")
latitude = read(hfile["latitude"])
longitude = read(hfile["longitude"])
metric = read(hfile["metric"])
close(hfile)
fmetric = reshape(metric, (192*96, 1))

index_1 = [157, 46]
index_2 = [54, 48]
index_3 = [19, 87]

function location_1(field)
    rfield = reshape(field, (192, 96))
    return rfield[157, 46]
end

function location_2(field)
    rfield = reshape(field, (192, 96))
    return rfield[54, 48]
end

function location_3(field)
    rfield = reshape(field, (192, 96))
    return rfield[19, 87]
end

observables = [location_1, location_3]

##
fig = Figure(; resolution)
observable_names_base = [@sprintf("%.2fᵒ, %.2fᵒ ", longitude[tmp[1]], latitude[tmp[2]]) for tmp in [index_1, index_2, index_3]]
observable_names = observable_names_base .* ["(Jul)"]
month_string = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
xticks = (1:12, month_string)

for i in eachindex(observables)
    month = 7
    emulator.month[1] = month
    emulator.global_mean_temperature[1] = max_temp585
    emulator_hurs.month[1] = month
    emulator_hurs.global_mean_temperature[1] = max_temp585
    tasmean, tasstd = emulator_mean_variance_linear_functionals(observables, emulator)
    hursmean, hursstd = emulator_mean_variance_linear_functionals(observables, emulator_hurs)

    common_options_1 = (; xlabel = "Temperature (K)")
    common_options_2 = (; xlabel = "Relative Humidity (%)")

    ax1 = Axis(fig[1, 2*i]; title = observable_names[i], common_options_1..., global_common_options..., ylabel = "PDF")
    xs, ys = gaussian_grid(tasmean[i], tasstd[i])
    lines!(ax1, xs, ys, color = :blue, linewidth = lw, label = "SSP5-8.5")

    ax2 = Axis(fig[2, 2*i]; title = observable_names[i], common_options_2..., global_common_options..., ylabel = "PDF")
    xs, ys = gaussian_grid(hursmean[i], hursstd[i])
    lines!(ax2, xs, ys, color = :blue, linewidth = lw)

    emulator.month[1] = month
    emulator.global_mean_temperature[1] = max_temp
    emulator_hurs.month[1] = month
    emulator_hurs.global_mean_temperature[1] = max_temp
    tasmean, tasstd = emulator_mean_variance_linear_functionals(observables, emulator)
    hursmean, hursstd = emulator_mean_variance_linear_functionals(observables, emulator_hurs)

    xs, ys = gaussian_grid(tasmean[i], tasstd[i])
    lines!(ax1, xs, ys, color = :orange, linewidth = lw, label = "New Scenario")
    xs, ys = gaussian_grid(hursmean[i], hursstd[i])
    lines!(ax2, xs, ys, color = :orange, linewidth = lw)

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
        emulator.global_mean_temperature[1] = max_temp
        tasmean, tasstd = emulator_mean_variance_linear_functionals([observables[i]], emulator)
        emulator_hurs.month[1] = month
        emulator_hurs.global_mean_temperature[1] = max_temp
        hursmean, hursstd = emulator_mean_variance_linear_functionals([observables[i]], emulator_hurs)
        
        tasmeanlist[month] = tasmean[1]
        hursmeanlist[month] = hursmean[1]
        tasstdlist[month] = tasstd[1]
        hursstdlist[month] = hursstd[1]
    end
    lines!(ax3, 1:12, tasmeanlist, color = :orange, linewidth = lw)
    band!(ax3, 1:12, tasmeanlist .- 3 .* tasstdlist, tasmeanlist .+ 3 .* tasstdlist, color = (:orange, 0.2))
    lines!(ax4, 1:12, hursmeanlist, color = :orange, linewidth = lw)
    band!(ax4, 1:12, hursmeanlist .- 3 .* hursstdlist, hursmeanlist .+ 3 .* hursstdlist, color = (:orange, 0.2))

end
display(fig)
save(figure_directory * "case_study.png", fig)