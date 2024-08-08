using NCDatasets, LinearAlgebra, Statistics, HDF5, ProgressBars

save_directory = "/net/fs06/d3/sandre/GaussianEarthData/"
data_directory = "/net/fs06/d3/mgeo/CMIP6/interim/"

scenario_directories = readdir(data_directory)
field_name = "tas"

gmt = []
regression_file = h5open(save_directory * field_name * "_ensemble_yearly_average.hdf5", "w")
for scenario in ProgressBar(scenario_directories)
    hfile = h5open(save_directory * field_name * "_" * scenario * "_projection.hdf5", "r")
    global_mean = read(hfile["global mean"])
    close(hfile)
    ensemble_average = mean(global_mean, dims = 2)[:]
    months = size(ensemble_average)[1] ÷ 12
    yearly_average = mean(reshape(ensemble_average, (12, months)), dims = 1)[:]
    regression_file[scenario] = yearly_average
end
close(regression_file)