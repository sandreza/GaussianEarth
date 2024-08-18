using CairoMakie
using NCDatasets, LinearAlgebra, Statistics, HDF5, ProgressBars
using LinearAlgebra

include("utils.jl")
##
save_directory = "/net/fs06/d3/sandre/GaussianEarthData/"
data_directory = "/net/fs06/d3/mgeo/CMIP6/interim/"
scenario_directories = readdir(data_directory)
current_path = joinpath(data_directory, scenario_directories[1])
variable_directories = readdir(current_path)

##
field = "tas"
for month in ProgressBar(1:12)
covsave = h5open(save_directory * field * "_covariances.hdf5", "r")
cov1 = covsave["covariances $month"][:,:,:]
temperature = covsave["temperature"][:]
scale = read(covsave["scale"])
close(covsave)
covsave = h5open(save_directory * field * "_covariances_model.hdf5", "r")
L = covsave["L" * "$month"][:,:, :]
close(covsave)

Σ⁰ = L[:,:,1]' * L[:,:,1]
Σ¹ = L[:,:,1]' * L[:,:,2] + L[:,:,2]' * L[:,:,1]
Σ² = L[:,:,2]' * L[:,:,2]

function model(L, temperature)
    L⁰ = L[:,:,1]
    L¹ = L[:,:,2]
    L = L⁰ + L¹ * temperature
    return L' * L 
end
function model(Σ⁰, Σ¹, Σ², temperature)
    Σ = Σ⁰ + Σ¹ * temperature + Σ² * temperature^2
    return Σ 
end

model(L, temperature[1])

##
indlist = [[1,1], [2, 2], [3, 3], [4, 4], [1, 3], [10, 10], [50, 50], [100, 100], [1000, 1000]]
fig = Figure() 
i = 1
for i in 1:9
    ii = (i - 1) ÷ 3 + 1
    jj = (i - 1) % 3 + 1
    a, b = indlist[i]
    ax = Axis(fig[ii, jj]; title = "C($a, $b)")
    covs = [cov1[a, b, i] for i in eachindex(temperature)]
    covs_model = [model(Σ⁰[a, b], Σ¹[a, b], Σ²[a, b], t) for t in temperature]
    scatter!(ax, temperature * scale, covs, color = (:blue, 0.5))
    lines!(ax, temperature * scale, covs_model, color = :red)
end
display(fig)
save("covariance_model_fit_month_$month.png", fig)
end
