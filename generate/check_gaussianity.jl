using CairoMakie, Statistics

field = "pr"
eof_mode, temperature = concatenate_regression(field, ["historical"])
eofs = eof_mode[:,:, [1:48..., 50:50...]]
##
function kurtosis(x)
    n = length(x)
    μ = mean(x)
    σ = std(x)
    return sum((x .- μ).^4) / (n * σ^4) - 3
end
function skewness(x)
    n = length(x)
    μ = mean(x)
    σ = std(x)
    return sum((x .- μ).^3) / (n * σ^3)
end
##
fig = Figure(resolution = (1500, 1500))
inds = 1:30
N = 5
temp = temperature[inds] * 273
gaussian(x, μ, σ) = exp(-0.5 * ((x .- μ) ./ σ).^2) ./ sqrt(2π * σ^2)
lower_order_statistics = Vector{Float64}[]
for i in 1:N^2
    ii = div(i-1, N) + 1
    jj = mod(i-1, N) + 1
    ax = Axis(fig[ii, jj]; title = "Mode $(i)")
    month_eof = eofs[i, 1:12:end, :]
    μ = mean(month_eof[inds, :])
    σ = std(month_eof[inds, :])
    κ = kurtosis(month_eof[inds, :])
    ρ = skewness(month_eof[inds, :])
    push!(lower_order_statistics, [μ, σ, κ, ρ])
    x = range(extrema(month_eof[inds, :])..., length = 100)
    y = gaussian.(x, μ, σ)
    hist!(ax, month_eof[inds, :][:], bins = 20, color = (:blue, 0.5), normalization = :pdf)
    lines!(ax, x, y, color = :red)
end
display(fig)
save(field * "_eof_.png", fig)
##
κ = [lower_order_statistic[3] for lower_order_statistic in lower_order_statistics]
ρ = [lower_order_statistic[4] for lower_order_statistic in lower_order_statistics]

fig = Figure(resolution = (1500, 300))
ax = Axis(fig[1, 1]; title = "Kurtosis", xlabel = "Mode", ylabel = "Kurtosis")
scatter!(ax, κ, color = :blue)
ylims!(ax, (-1,1))
ax = Axis(fig[1, 2]; title = "Skewness", xlabel = "Mode", ylabel = "Skewness")
scatter!(ax, ρ, color = :red)
ylims!(ax, (-1,1))
display(fig)

save(field * "_skew_kurt_.png", fig)