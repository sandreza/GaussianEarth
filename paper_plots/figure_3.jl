ts = 40
xls = 40 
yls = 40
tls = 40
legend_ls = 35
resolution = (3000, 1000)
common_options = (; titlesize = ts, xlabelsize = xls, ylabelsize = yls, xticklabelsize = tls, yticklabelsize = tls)

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

##
scale = 273
month = 1
eof_mode, temperature = concatenate_regression(field, ["historical", "ssp585"])
eofs = eof_mode[:,month:12:end, 1:45] 

hfile = h5open(save_directory * field * "_mean_regression_quadratic.hdf5", "r")
regression_coefficients = read(hfile["regression_coefficients 1"])
quadratic_coefficients = zeros(Float32, size(regression_coefficients)..., 12)
for month in ProgressBar(1:12)
    quadratic_coefficients[:, :, month] .= read(hfile["regression_coefficients $month"])
end
close(hfile)

##
fig = Figure(; resolution)
eof_index = 1
for (jj, eof_index) in enumerate([1, 10, 100, 1000])
    if jj == 1
        ax = Axis(fig[1,jj]; title = "Mode $eof_index", xlabel = "Temperature (K)", ylabel = "Amplitude", common_options...)
    else
        ax = Axis(fig[1,jj]; title = "Mode $eof_index", xlabel = "Temperature (K)", common_options...)
    end
    month_eof = eofs[eof_index, :, :]
    means = mean(month_eof, dims = 2)
    # line_fit = linear_coefficients[eof_index, 1, month] .+ linear_coefficients[eof_index, 2, month] * temperature
    quad_fit = quadratic_coefficients[eof_index, 1, month] .+ quadratic_coefficients[eof_index, 2, month] * temperature .+ quadratic_coefficients[eof_index, 3, month] * temperature.^2
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
    scatter!(ax, temperature * scale, factor * means[:], color = (:black, 0.5))
    # lines!(ax, temperature * scale, factor * line_fit, color = :blue, label = "Emulator")
    lines!(ax, temperature * scale, factor * quad_fit, color = :blue, label = "Emulator")
    if jj == 1
        axislegend(ax, position = :lt, labelsize = legend_ls)
    end
end

for (jj, ab) in enumerate([[1, 1], [10, 10], [100, 100], [2, 15]])
    a = ab[1]
    b = ab[2]
    if jj == 1
        ax = Axis(fig[2, jj]; title = "Cov(Mode $a, Mode $b)", xlabel = "Temperature (K)", ylabel = "Covariance", common_options...)
    else
        ax = Axis(fig[2, jj]; title = "Cov(Mode $a, Mode $b)", xlabel = "Temperature (K)", common_options...)
    end
    covs = [cov1[a, b, i] for i in eachindex(temperature_cov)]
    covs_model = [model(Σ⁰[a, b], Σ¹[a, b], Σ²[a, b], t) for t in temperature]
    scatter!(ax, temperature_cov * scale, covs, color = (:purple, 0.3))
    lines!(ax, temperature * scale, covs_model, color = :blue)
end

save(figure_directory * "figure_3_regression_check.png", fig)
 
@info "Figure 3 generated"