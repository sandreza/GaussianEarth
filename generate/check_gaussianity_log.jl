using CairoMakie, Statistics, ProgressBars

include("utils.jl")

inds = 1:30
fields =  ["hurs", "huss", "pr", "tas"]
for field in fields
    eof_mode, temperature = concatenate_log_regression(field, ["historical"])
    eofs = eof_mode[:,:, inds] # [1:48..., 50:50...]]


    fig = Figure(resolution = (1500, 1500))
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
    save(field * "_log_eof.png", fig)

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

    save(field * "_log_skew_kurt.png", fig)
end

##
modes = 100
fig = Figure(resolution = (1500, 1500))
field = fields[4]
for (i, field) in enumerate(fields)
    eof_mode, temperature = concatenate_log_regression(field, ["historical"])
    if i !== 2
        eofs = eof_mode[:,:, 1:30] 
    else
        eofs = eof_mode[:,:, 1:47] 
    end
    lower_order_statistics = Vector{Float64}[]
    for i in 1:1980
        month_eof = eofs[i, 1:12:end, :]
        μ = mean(month_eof[inds, :])
        σ = std(month_eof[inds, :])
        κ = kurtosis(month_eof[inds, :])
        ρ = skewness(month_eof[inds, :])
        push!(lower_order_statistics, [μ, σ, κ, ρ])
    end

    κ = [lower_order_statistic[3] for lower_order_statistic in lower_order_statistics]
    ρ = [lower_order_statistic[4] for lower_order_statistic in lower_order_statistics]


    ax = Axis(fig[i, 1]; title = "log " * field, xlabel = "Mode", ylabel = "Kurtosis")
    scatter!(ax, κ[1:modes], color = :blue)
    ylims!(ax, (-1,1))
    ax = Axis(fig[i, 2]; title = "log " * field, xlabel = "Mode", ylabel = "Skewness")
    scatter!(ax, ρ[1:modes], color = :red)
    ylims!(ax, (-1,1))
end
display(fig)
save("all_fields_skew_kurt_log.png", fig)