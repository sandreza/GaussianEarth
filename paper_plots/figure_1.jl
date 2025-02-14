scenarios = ["historical", "ssp119", "ssp245", "ssp585"]
time_history = 1850:2014
time_future = 2015:2100
scenario_colors = Dict("historical" => :red4, "ssp585" => :red, "ssp245" => :magenta3, "ssp119" => :indigo)
scenario_labels = ["Historical","SSP1-1.9",  "SSP2-4.5", "SSP5-8.5"]

baseline = 0.0
fig = Figure(resolution=(600, 400))
ax = Axis(fig[1,1], xlabel="Year", ylabel="Global Mean Temperature (K)", xticks=1850:50:2100, title="Global Mean Temperature in MPI-ESM1.2-LR")
for (e, scenario) in enumerate(scenarios)
    hfile = h5open(save_directory * field_name * "_" * scenario * "_projection.hdf5", "r")
    global_mean = read(hfile["global mean"])
    close(hfile)
    annuals = zeros(Int(size(global_mean)[1]/12), size(global_mean)[2])
    for i in 1:size(global_mean)[2] #num_ens_members
        gmt = global_mean[:,i]
        annual_avg = mean(reshape(gmt, 12, :), dims=1)[:]
        if any(x -> x < 200, annual_avg)
            println("Warning: $scenario ensemble memeber $i has a value less than 200")
            continue
        end
        annuals[:,i] = annual_avg
        lines!(ax, (scenario=="historical" ? time_history : time_future), annual_avg, color=scenario_colors[scenario], alpha=0.2, linestyle=:dot)
    end
    annuals = annuals[:, 1:end .!= 49 ]
    mean_gmt = mean(annuals, dims=2)
    if scenario == "historical"
        baseline = mean_gmt[1]
    end
    lines!(ax, (scenario=="historical" ? time_history : time_future), mean_gmt[:], color=scenario_colors[scenario], alpha=1, linestyle=:solid, label=scenario_labels[e])
end


axislegend(ax, position=:lt)
save(figure_directory*"figure_1.png", fig)
# display(fig)

