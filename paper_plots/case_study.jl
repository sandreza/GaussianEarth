
##
field = "tas" 
field_name = "tas"
# field = "hurs" 
# field_name = "hurs"

#load in saved out linear emulator
include(:"emulator.jl")
d = 1000

# load in arbitrary scenario
hfile = h5open("arbitrary_scenario.hdf5", "r")
temperatures = read(hfile["arbitrary"]) 
close(hfile)


# get metric correction terms
if process_data
    hfile = h5open(save_directory * field_name * "_basis.hdf5", "r")
    latitude = read(hfile["latitude"])
    longitude = read(hfile["longitude"])
    metric = read(hfile["metric"])
    close(hfile)
    sqrt_f_metric = sqrt.(reshape(metric, 192 * 96))
end

month = 1
# get eof model
if process_data
    hfile = h5open(save_directory * field_name * "_mean_regression.hdf5", "r")
    regression_coefficients = read(hfile["regression_coefficients 1"])
    linear_coefficients = zeros(Float32, size(regression_coefficients)..., 12)
    for month in ProgressBar(1:12)
        linear_coefficients[:, :, month] .= read(hfile["regression_coefficients $month"])
    end
    close(hfile)

    Φ = eof_basis(field_name) 
    basis = Φ[:, 1:d]

    # averaged_coefficients = mean(linear_coefficients, dims = 3)[:, :, 1] #using the annual-avg coefficients !!!!!
    coeffs = linear_coefficients[:, :, month]

    eof_model = zeros(Float32, size(Φ)[1]..., 2)
    for i in ProgressBar(1:d)
        # eof_model[:, 1] .+= Φ[:, i] * averaged_coefficients[i, 1]
        # eof_model[:, 2] .+= Φ[:, i] * averaged_coefficients[i, 2]
        eof_model[:, 1] .+= Φ[:, i] * coeffs[i, 1]
        eof_model[:, 2] .+= Φ[:, i] * coeffs[i, 2]
    end
end

### apply emulator to temperatures


# means = zeros(Float32, size(Φ)[1], length(temperatures), 12)
# stds = zeros(Float32, size(Φ)[1], length(temperatures), 12)

# for (i, t) in ProgressBar(enumerate(temperatures))
#     month = 1
#     # for month in 1:12
#         μ = eof_model[:, 1] .+ t * eof_model[:, 2]
#         Σ = emulator_variance(emulator; modes=d, month=month, global_mean_temperature=t)
#         σ = sqrt.(diag((basis * Σ) * basis'))
#     # end
# end

# get plots of interest
scenario = "ssp245"
temperatures = regression_variable(scenario)
# year = 2100
t = temperatures[end]
# these are already in temperatures
μ = eof_model[:, 1] .+ t * eof_model[:, 2]
Σ = emulator_variance(emulator; modes=d, month=month, global_mean_temperature=t)
σ = sqrt.(diag((basis * Σ) * basis'))

field = reshape(μ, (192, 96))
field_upper = reshape(μ .+ 2 .* σ, (192, 96))
field_lower = reshape(μ .- 2 .* σ, (192, 96))




### plot
ts = 60
xls = 60 
yls = 60
tls = 60
legend_ls = 60
xlp = 20 
ylp = 20 
resolution = (1800, 1200) .* 2
common_options = (; titlesize = ts, xlabelsize = xls, ylabelsize = yls, xticklabelsize = tls, yticklabelsize = tls, xlabelpadding = xlp, ylabelpadding = ylp)


##
if field == "tas"
    colormap = :matter
    crange = (minimum(field_lower), maximum(field_upper))
    fig = Figure(; resolution)
    ax = GeoAxis(fig[1,1]; title = "2100 January Temperature -2σ", common_options...)
    nlongitude = range(-180, 180, length = 192)
    ndata_realization = circshift(field_lower, (96, 0))
    surface!(ax, nlongitude, latitude, ndata_realization; colormap = colormap, colorrange = crange,  shading = NoShading) 
    hidedecorations!(ax)

    ax = GeoAxis(fig[1,2]; title = "2100 January Temperature +2σ", common_options...)
    nlongitude = range(-180, 180, length = 192)
    ndata_realization = circshift(field_upper, (96, 0))
    surface!(ax, nlongitude, latitude, ndata_realization; colormap = colormap, colorrange = crange,  shading = NoShading) 
    hidedecorations!(ax)

    Colorbar(fig[1,3], label = "Temperature (K)", colorrange = crange, colormap = colormap, height = Relative(2/4), labelsize = legend_ls, ticklabelsize = legend_ls)
    display(fig)
    save(figure_directory * "tas_map_comparison.png", fig)

elseif field == "hurs"
    colormap = Reverse(:pink)
    crange = (minimum(field_lower), maximum(field_upper))
    fig = Figure(; resolution)
    ax = GeoAxis(fig[1,1]; title = "2100 January Rel. Humidity -2σ", common_options...)
    nlongitude = range(-180, 180, length = 192)
    ndata_realization = circshift(field_lower, (96, 0))
    surface!(ax, nlongitude, latitude, ndata_realization; colormap = colormap, colorrange = crange,  shading = NoShading) 
    hidedecorations!(ax)

    ax = GeoAxis(fig[1,2]; title = "2100 January Rel. Humidity +2σ", common_options...)
    nlongitude = range(-180, 180, length = 192)
    ndata_realization = circshift(field_upper, (96, 0))
    surface!(ax, nlongitude, latitude, ndata_realization; colormap = colormap, colorrange = crange,  shading = NoShading) 
    hidedecorations!(ax)

    Colorbar(fig[1,3], label = "Rel. Humidity (%)", colorrange = crange, colormap = colormap, height = Relative(2/4), labelsize = legend_ls, ticklabelsize = legend_ls)
    display(fig)
    save(figure_directory * "hurs_map_comparison.png", fig)
end



# ax = GeoAxis(fig[1,2]; title = "2100 January Rel. Humidity -2σ", common_options...)
# nlongitude = range(-180, 180, length = 192)
# ndata_realization = circshift(field_lower, (96, 0))
# surface!(ax, nlongitude, latitude, ndata_realization; colormap = colormap, colorrange = crange,  shading = NoShading) 
# hidedecorations!(ax)

# ax = GeoAxis(fig[2,2]; title = "2100 January Rel. Humidity +2σ", common_options...)
# nlongitude = range(-180, 180, length = 192)
# ndata_realization = circshift(field_upper, (96, 0))
# surface!(ax, nlongitude, latitude, ndata_realization; colormap = colormap, colorrange = crange,  shading = NoShading) 
# hidedecorations!(ax)



# save(figure_directory * "tas_realization_comparison.png", fig)