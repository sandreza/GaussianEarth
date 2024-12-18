ts = 60
xls = 70 
yls = 70
tls = 70
legend_ls = 60
resolution = (400, 150) .* 12
lw = 7
global_common_options = (; titlesize = ts, xlabelsize = xls, ylabelsize = yls, xticklabelsize = tls, yticklabelsize = tls)

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
month_names = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
for i in eachindex(observables)
    month = 7
    emulator.month[1] = month
    emulator.global_mean_temperature[1] = min_temp
    emulator_hurs.month[1] = month
    emulator_hurs.global_mean_temperature[1] = min_temp
    tasmean, tasstd = emulator_mean_variance_linear_functionals(observables, emulator)
    hursmean, hursstd = emulator_mean_variance_linear_functionals(observables, emulator_hurs)


    common_options_1 = (; xlabel = "Temperature (K)")
    common_options_2 = (; xlabel = "Relative Humidity (%)")

    ax1 = Axis(fig[1, 2*i]; title = observable_names[i], common_options_1..., global_common_options..., ylabel = "PDF")
    xs, ys = gaussian_grid(tasmean[i], tasstd[i])
    lines!(ax1, xs, ys, color = :blue, linewidth = lw, label = "2015")

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
    lines!(ax1, xs, ys, color = :orange, linewidth = lw, label = "2100")
    xs, ys = gaussian_grid(hursmean[i], hursstd[i])
    lines!(ax2, xs, ys, color = :orange, linewidth = lw)

    if i == 1
        axislegend(ax1, position = :lt, labelsize = legend_ls)
    end

    ax3 = Axis(fig[1, 2*i-1]; title = observable_names_base[i] * "(2100)", common_options_1..., global_common_options..., ylabel = "Temperature (K)")
    ax4 = Axis(fig[2, 2*i-1]; title = observable_names_base[i] * "(2100)", common_options_2..., global_common_options..., ylabel = "Relative Humidity (%)")
    for n in 1:20
        taslist = zeros(12)
        hurslist = zeros(12)
        for month in 1:12
            emulator.month[1] = month
            emulator.global_mean_temperature[1] = max_temp
            emulator_hurs.month[1] = month
            emulator_hurs.global_mean_temperature[1] = max_temp
            tas_instance = rand(emulator)
            hurs_instance = rand(emulator_hurs)
            
            taslist[month] = observables[i](tas_instance)
            hurslist[month] = observables[i](hurs_instance)
        end
        lines!(ax3, month_names, taslist, alpha = 0.7)
        lines!(ax4, month_names, hurslist, alpha = 0.7)
    end

end
display(fig)
save(figure_directory * "case_study.png", fig)