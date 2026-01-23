ts = 60
xls = 60 
yls = 60
tls = 60
legend_ls = 60
xlp = 20 
ylp = 20 
resolution = (1800, 1200) .* 2
common_options = (; titlesize = ts, xlabelsize = xls, ylabelsize = yls, xticklabelsize = tls, yticklabelsize = tls, xlabelpadding = xlp, ylabelpadding = ylp)

# calculate errors for quadratic emulator
field_name = "tas"
hfile = h5open(save_directory * field_name * "_basis.hdf5", "r")
latitude = read(hfile["latitude"])
longitude = read(hfile["longitude"])
metric = read(hfile["metric"])
close(hfile)
sqrt_f_metric = sqrt.(reshape(metric, 192 * 96))

#### now for stds
## get true std data for comparison
if process_data
    scenarios = ["historical", "ssp585", "ssp119", "ssp245"]
    temperatures = [] 
    fields = [] #list of four arrays of time x space data values for given variable
    for scenario in scenarios
        temperature = regression_variable(scenario) #this gets the list of temps to regress onto
        a, b = ensemble_averaging(scenario, field_name; ensemble_members = 45, return_std=true) #gets the ensemble avg for that var? #changed this 
        push!(temperatures, temperature)
        push!(fields, a[:,:,:,:]) # stds for all months
    end
end

mean_stds = zeros(Float32, 192, 96, 12, length(scenarios))
for scenario_index in 1:4
    for month in 1:12
        mean_stds[:, :, month, scenario_index] = mean(fields[scenario_index][:, :, :, month], dims=3)
    end
end

year_mean_stds = mean(mean_stds, dims=3)[:,:,1,:] # average over the 12 months


#######

fig = Figure(; resolution)
nlongitude = longitude .- 180
scenario_index = 1
cmap = Reverse(:linear_tritanopic_krjcw_5_98_c46_n256) #:plasma
# cmap = :balance
ax = GeoAxis(fig[1,1]; title = "Historical", common_options...)
field = year_mean_stds[:, :,  scenario_index]
shifted_field = circshift(field, (96, 0))
crange = extrema(year_mean_stds)
surface!(ax, nlongitude, latitude, shifted_field; colormap = cmap, colorrange = crange, shading = NoShading)
hidedecorations!(ax)
scenario_index = 2 
ax = GeoAxis(fig[1,2]; title = "SSP5-8.5", common_options...)
field = year_mean_stds[:, :,  scenario_index]
shifted_field = circshift(field, (96, 0))
surface!(ax, nlongitude, latitude, shifted_field; colormap = cmap, colorrange = crange, shading = NoShading)
hidedecorations!(ax)
scenario_index = 3
ax = GeoAxis(fig[2,1]; title = "SSP1-1.9", common_options...)
field = year_mean_stds[:, :,  scenario_index]
shifted_field = circshift(field, (96, 0))
surface!(ax, nlongitude, latitude, shifted_field; colormap = cmap, colorrange = crange, shading = NoShading)
hidedecorations!(ax)
scenario_index = 4
ax = GeoAxis(fig[2,2]; title = "SSP2-4.5", common_options...)
field = year_mean_stds[:, :,  scenario_index]
shifted_field = circshift(field, (96, 0))
surface!(ax, nlongitude, latitude, shifted_field; colormap = cmap, colorrange = crange, shading = NoShading)
Colorbar(fig[1:2,3], colormap=cmap, colorrange=crange, height = Relative(2/4), label = "temperature standard deviation (K)", labelsize = legend_ls, ticklabelsize = legend_ls)
hidedecorations!(ax)
save(figure_directory * "figure_5_tas_averaged_std.png", fig)
# display(fig)