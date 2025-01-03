using NCDatasets, LinearAlgebra, Statistics, HDF5, ProgressBars, CairoMakie, Interpolations, ProgressBars

lon = -71.0589
lat = 42.3601
location_string = "Boston"


location_string = "Capita_Bermudez"
lat = -32.82
lon = -60.72

save_directory = "/net/fs06/d3/sandre/GaussianEarthData/"
data_directory = "/net/fs06/d3/mgeo/CMIP6/interim/"

scenario_directories = readdir(data_directory)
scenario_directory = scenario_directories[1]

current_path = joinpath(data_directory, scenario_directory)
variable_directories = readdir(current_path)
variable_directory = variable_directories[end]

local_current_path = joinpath(current_path, variable_directory)
file_names = readdir(local_current_path)

save_directory = "/net/fs06/d3/sandre/GaussianEarthData/"
hfile = h5open(save_directory * "tas" * "_basis.hdf5", "r")
metric = read(hfile["metric"] )
close(hfile)


location_temp = zeros(1980, length(file_names))
gmt = zeros(1980, length(file_names))
for i in ProgressBar(eachindex(file_names)[1:48])
    file_name = file_names[i] # pick the first file for computing a basis
    file_path = joinpath(local_current_path, file_name)

    ds = Dataset(file_path)
    field_name = keys(ds)[end] # always the last key for variable of interest
    field = Float32.(ds[field_name][:,:,:])

    # yep
    tmp = ds["tas"][:, :, :]

    # Data
    M, N, T = size(tmp)
    shifted_tmp = circshift(tmp, (-M ÷ 2, 0, 0))
    x = range(-180, 180, length=M)
    y = range(-90, 90, length=N)
    if i == 1
        fig = Figure() 
        ax = Axis(fig[1, 1])
        heatmap!(ax, x, y, shifted_tmp[:, :, 1])
        scatter!(ax, [lon], [lat], color = :red)
        save(location_string * "_check.png", fig)
    end

    # Interpolation object (caches coefficients and such)
    for j in ProgressBar(1:T)
        field= shifted_tmp[:, :, j]
        itp = cubic_spline_interpolation((x, y), field)
        location_temp[j, i] = itp(lon, lat) 
        gmt[j, i] = sum(field .* metric)
    end

end

##
hfile = h5open(location_string * "_historical.hdf5", "w")
hfile["temp"] = location_temp[:, 1:48]
hfile["gmt"] = gmt[:, 1:48]
close(hfile)

##
hfile = h5open(location_string * "_historical.hdf5", "r")
location_temp = read(hfile["temp"] )
gmt = read(hfile["gmt"] )
close(hfile)

fig = Figure(resolution = (2000, 2000)) 
for i in 1:4
    rr = rand(1:48)
    ax = Axis(fig[i, 1]; title = location_string * " temperature")
    lines!(ax, location_temp[:, rr])
    ax = Axis(fig[i, 2]; title = "GMT")
    lines!(ax, gmt[:, rr])
end
save(location_string * "_temp_gmt.png", fig)
