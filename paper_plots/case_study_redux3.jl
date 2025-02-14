#the actual code for the case study figure

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
new_scenario = read(hfile["new_scenario"])
close(hfile)
new_scenario_temperatures = new_scenario ./273

ssp5 = regression_variable("ssp585").*273   
ssp5_temperatures = regression_variable("ssp585")
ssp2 = regression_variable("ssp245").*273
ssp2_temperatures = regression_variable("ssp245")

min_temp_1 = new_scenario[1]/273 #2030
max_temp_1 = new_scenario[end]/273 #2070

min_temp_2 = ssp5[1]/273 #2030
max_temp_2 = ssp5[end]/273 #2070

##
field = "tas" 
field_name = "tas"

#load in saved out linear emulator
include("emulator.jl")
include("emulator_hurs.jl")
d = 1000


hfile = h5open(save_directory * field * "_basis.hdf5", "r")
latitude = read(hfile["latitude"])
longitude = read(hfile["longitude"])
metric = read(hfile["metric"])
close(hfile)
fmetric = reshape(metric, (192*96, 1))

# define observables
begin
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

# function nw_india(field) #this one is wrong, it's in the bahamas lol
#     rfield = reshape(field, (192, 96))
#     return rfield[148, 61] #75W and 23N
# end

    function nw_india(field) #CORRECT ONE
        rfield = reshape(field, (192, 96))
        return rfield[46, 61]
    end
        

end
observables = [nw_india, chicago]


#get divergences
divs_tas = zeros(Float32, (length(observables), length(new_scenario)))
divs_hurs = zeros(Float32, (length(observables), length(new_scenario)))

month = 1
for (i, year) in ProgressBar(enumerate(2015:2100))
    emulator.month[1] = month
    emulator_hurs.month[1] = month
    temp = new_scenario_temperatures[i]
    emulator.global_mean_temperature[1] = temp
    emulator_hurs.global_mean_temperature[1] = temp
    tasmean, tasstd = emulator_mean_variance_linear_functionals(observables, emulator; show_progress = false)
    hursmean, hursstd = emulator_mean_variance_linear_functionals(observables, emulator_hurs; show_progress = false)

    temp = ssp5_temperatures[i]
    emulator.global_mean_temperature[1] = temp
    emulator_hurs.global_mean_temperature[1] = temp
    tasmean_ssp5, tasstd_ssp5 = emulator_mean_variance_linear_functionals(observables, emulator; show_progress = false)
    hursmean_ssp5, hursstd_ssp5 = emulator_mean_variance_linear_functionals(observables, emulator_hurs; show_progress = false)

    for j in eachindex(observables)
        # divergence of new scenario from ssp5
        divs_tas[j, i] = kl_div(tasmean[j], tasstd[j], tasmean_ssp5[j], tasstd_ssp5[j])
        divs_hurs[j, i] = kl_div(hursmean[j], hursstd[j], hursmean_ssp5[j], hursstd_ssp5[j])
    end
end


fig = Figure(; resolution)
observable_names_base = ["NW India", "Midwest US"]
month_string = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
observable_names = observable_names_base  .* [" (" * month_string[month] * ")"]
xticks = (1:12, month_string)
num_stds = 2
common_options_1 = (;)# xlabel = "Temperature (K)")
common_options_2 = (;)# xlabel = "Relative Humidity (%)")
ax0 = Axis(fig[1,1]; common_options_1..., global_common_options..., xlabel="Year", ylabel="GMT (K)", title="Global Mean Temperature") # xticks=1850:50:2100
lines!(ax0, 2015:2100, new_scenario, color = :blue, linewidth = lw, label = "New Scenario")
lines!(ax0, 2015:2100, ssp5, color = :red, linewidth = lw, label = "SSP5-8.5")
lines!(ax0, 2015:2100, ssp2, color = :magenta3, linewidth = lw,linestyle=:dash, label = "SSP2-4.5")
axislegend(ax0, position = :lt, labelsize = legend_ls)

for i in eachindex(observables)
    if i == 2
        ax1 = Axis(fig[1, 2*i-1]; title = month_string[month], common_options_1..., global_common_options..., ylabel = "D")
        ax2 = Axis(fig[2, 2*i-1]; title = month_string[month], common_options_2..., global_common_options..., ylabel = "D")
        
        # plot divergences
        lines!(ax1, 2015:2100, divs_tas[1, :], color = :green, linewidth = lw, label = observable_names_base[1])
        lines!(ax1, 2015:2100, divs_tas[2, :], color = :purple, linewidth = lw, label =observable_names_base[2])
        lines!(ax2, 2015:2100, divs_hurs[1, :], color = :green, linewidth = lw, label = observable_names_base[1])
        lines!(ax2, 2015:2100, divs_hurs[2, :], color = :purple, linewidth = lw, label = observable_names_base[2])
    end

    if i == 2
        axislegend(ax1, position = :lt, labelsize = legend_ls)
    end

    ax3 = Axis(fig[1, 2*i]; title = observable_names_base[i] * " (2100)", xticks,xticklabelrotation=45.0, common_options_1..., global_common_options..., ylabel = "Temperature (K)")
    ax4 = Axis(fig[2, 2*i]; title = observable_names_base[i] * " (2100)", xticks,xticklabelrotation=45.0, common_options_2..., global_common_options..., ylabel = "Relative Humidity (%)")
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
    band!(ax3, 1:12, tasmeanlist .- num_stds .* tasstdlist, tasmeanlist .+ num_stds .* tasstdlist, color = (:blue, 0.2), label = "New Scenario")
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
    band!(ax3, 1:12, tasmeanlist .- num_stds .* tasstdlist, tasmeanlist .+ num_stds .* tasstdlist, color = (:orangered, 0.2), label = "SSP5-8.5")
    lines!(ax4, 1:12, hursmeanlist, color = :orangered, linewidth = lw)
    band!(ax4, 1:12, hursmeanlist .- num_stds .* hursstdlist, hursmeanlist .+ num_stds .* hursstdlist, color = (:orangered, 0.2))

    if i == 1
        axislegend(ax3, position = :rt, labelsize = legend_ls)
    end

end
display(fig)
save(figure_directory * "case_study_divergence.png", fig)