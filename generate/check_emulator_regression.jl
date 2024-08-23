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
month = 1
covsave = h5open(save_directory * field * "_covariances.hdf5", "r")
cov1 = covsave["covariances $month"][:,:,:]
temperature_cov = covsave["temperature"][:]
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
    covs = [cov1[a, b, i] for i in eachindex(temperature_cov)]
    covs_model = [model(Σ⁰[a, b], Σ¹[a, b], Σ²[a, b], t) for t in temperature_cov]
    scatter!(ax, temperature_cov * scale, covs, color = (:purple, 0.1))
    lines!(ax, temperature_cov * scale, covs_model, color = :blue)
end
display(fig)
save("covariance_model_fit_month_$month.png", fig)

##
scale = 273
month = 1
eof_mode, temperature = concatenate_regression(field, ["historical", "ssp585"])
eofs = eof_mode[:,month:12:end, 1:45] 

hfile = h5open(save_directory * field * "_mean_regression.hdf5", "r")
regression_coefficients = read(hfile["regression_coefficients 1"])
linear_coefficients = zeros(Float32, size(regression_coefficients)..., 12)
for month in ProgressBar(1:12)
    linear_coefficients[:, :, month] .= read(hfile["regression_coefficients $month"])
end
close(hfile)

##
fig = Figure(resolution = (1500, 600))
eof_index = 1
for (jj, eof_index) in enumerate([1, 10, 100, 1000])
    if jj == 1
        ax = Axis(fig[1,jj]; title = "Mode $eof_index", xlabel = "Temperature (K)", ylabel = "Amplitude")
    else
        ax = Axis(fig[1,jj]; title = "Mode $eof_index", xlabel = "Temperature (K)")
    end
    month_eof = eofs[eof_index, :, :]
    line_fit = linear_coefficients[eof_index, 1, month] .+ linear_coefficients[eof_index, 2, month] * temperature
    restructured = [month_eof[:, i] for i in 1:45]
    if jj == 1
        factor = -1
    else
        factor = 1
    end
    scatter!(ax, temperature[1:2] * 273, factor * restructured[1][1:2], color = (:purple, 0.5), label = "Data")
    for i in 1:45
        scatter!(ax, temperature * scale, factor * restructured[i], color = (:purple, 0.03))
    end
    lines!(ax, temperature * scale, factor * line_fit, color = :blue, label = "Emulator")
    if jj == 1
        axislegend(ax, position = :lt)
    end
end

for (jj, ab) in enumerate([[1, 1], [10, 10], [100, 100], [2, 15]])
    a = ab[1]
    b = ab[2]
    if jj == 1
        ax = Axis(fig[2, jj]; title = "Cov(Mode $a, Mode $b)", xlabel = "Temperature (K)", ylabel = "Covariance")
    else
        ax = Axis(fig[2, jj]; title = "Cov(Mode $a, Mode $b)", xlabel = "Temperature (K)")
    end
    covs = [cov1[a, b, i] for i in eachindex(temperature_cov)]
    covs_model = [model(Σ⁰[a, b], Σ¹[a, b], Σ²[a, b], t) for t in temperature]
    scatter!(ax, temperature_cov * scale, covs, color = (:purple, 0.3))
    lines!(ax, temperature * scale, covs_model, color = :blue)
end

save("regression_check.png", fig)

##
Cstart = model(L, temperature[1])
Cend = model(L, temperature[end])

change_inds = (Cstart .> 0) .&& (Cend .< 0)
argmax(change_inds)
change_inds

reverse(sortperm(change_inds[:]))