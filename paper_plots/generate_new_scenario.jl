# version used in paper. creates one scenario to be compared to ssp585

time_future = 2015:2100

init_temp = regression_variable("historical")[end]*273
baseline = mean(regression_variable("historical")[1:51])*273
end_temp = baseline + 3

ssp5 = regression_variable("ssp585").*273   
ssp2 = regression_variable("ssp245").*273

new_scenario = zeros(length(time_future))
n= length(time_future) 
for t in 1:length(time_future)
    new_scenario[t] = init_temp * (end_temp / init_temp) ^ ((t - 1) / (n - 1)) 
end

# plot the scenarios
fig = Figure()
ax = Axis(fig[1,1], xlabel="Year", ylabel="Global Mean Temperature (K)", xticks=1850:50:2100)
lines!(ax, time_future, new_scenario, label="Scenario 1", color=:blue)
lines!(ax, time_future, ssp5, label="SSP5-8.5", color=:green)
lines!(ax, time_future, ssp2, label="SSP2-4.5", color=:purple, linestyle=:dash)
axislegend(ax, position=:lt)
display(fig)

hfile = h5open(save_directory * "new_scenarios.hdf5", "w")
write(hfile, "new_scenario", new_scenario)
close(hfile)