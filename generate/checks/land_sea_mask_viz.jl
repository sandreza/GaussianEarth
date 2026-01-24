using NCDatasets, CairoMakie, HDF5

# Load the land-sea mask
filename = "/net/fs06/d3/sandre/GaussianEarthData/IMERG_land_sea_mask.nc"
ds = Dataset(filename)
mask = ds["landseamask"][:, :]

target_resolution = (192, 96)
current_resolution = size(mask)

ratios = current_resolution ./ target_resolution
lonindices = range(0.5, target_resolution[1]-0.5, length = target_resolution[1])
latindices = range(0.5, target_resolution[2]-0.5, length = target_resolution[2])

subsampled_lon = floor.(Int, ratios[1] .* lonindices)
subsampled_lat = floor.(Int, ratios[2] .* latindices)
subsampled_mask = mask[subsampled_lon, subsampled_lat]

fig = Figure()
ax = Axis(fig[1, 1])
heatmap!(ax, mask)
ax = Axis(fig[1, 2])
heatmap!(ax, subsampled_mask .> 90)
save("land_sea_mask.png", fig)
display(fig)

mpi_mask = subsampled_mask .> 90

##
#=
filename = "/net/fs06/d3/sandre/GaussianEarthData/sftlf050_fx_MPI-ESM1-2-LR_piControl_r1i1p1f1_native.nc"
ds = Dataset(filename)
mask = ds["mask"][:, :]
lon = ds["lon"][:]
lat = ds["lat"][:]
close(ds)
=#

hfile = h5open("/net/fs06/d3/sandre/GaussianEarthData/land_sea_mask.hdf5", "w")
hfile["mask"] = Float32.(mask)
close(hfile)


index_1 = [157, 46]
index_2 = [54, 48]
index_3 = [19, 87]

fig = Figure()
ax = Axis(fig[1, 1])
heatmap!(ax, lon, lat, mask * 1.0, colormap = :viridis)
colors = (:red, :blue, :orange)
for (i, index) in enumerate([index_1, index_2, index_3])
    vlines!(ax, lon[index[1]], color = colors[i])
    hlines!(ax, lat[index[2]], color = colors[i])
end
save("mpi_land_sea_mask.png", fig)
display(fig)

