time_future = 2015:2100

init_temp = regression_variable("historical")[end]*273
baseline = mean(regression_variable("historical")[1:51])*273
# so limiting it to 2 degrees would be limitign it to 290
end_temp = baseline + 2

ssp5 = regression_variable("ssp585").*273
ssp2 = regression_variable("ssp245").*273

scenario_1 = zeros(length(time_future))
scenario_2 = zeros(length(time_future))

# Define the total temperature change
delta_T = end_temp - init_temp

# Define the speed parameters (k values for logarithmic growth)

k1 = 0.05  # Normal exponential rate for scenario 1
k2 = 0.05  # Faster exponential rate for scenario 2
interim_end_temp = ssp2[16]
n =16
for t in 1:16
    # scenario_1[t] = init_temp * (interim_end_temp / init_temp) ^ ((t - 1) / (n - 1)) ^ k1
    # scenario_2[t] = init_temp * (interim_end_temp / init_temp) ^ ((t - 1) / (n - 1)) ^ k2
    scenario_1[t] = init_temp + (interim_end_temp - init_temp) * (log(1 + k1 * t) / log(1 + k1 * n))
    scenario_2[t] = init_temp + (interim_end_temp - init_temp) * (log(1 + k2 * t) / log(1 + k2 * n))
end

# Generate temperature sequences using a logarithmic function
k1 = 0.02  # Speed of temperature increase for scenario 1
k2 = 0.02  # Slightly faster speed for scenario 2 (to create a difference)

n= length(time_future) - 16
init_temp = scenario_1[16]
for (t, i) in enumerate(17:length(time_future))
    scenario_1[i] = init_temp + (end_temp - init_temp) * (log(1 + k1 * t) / log(1 + k1 * n))
    scenario_2[i] = init_temp + (end_temp+0.3 - init_temp) * (log(1 + k2 * t) / log(1 + k2 * n))
end
# Print or plot the scenarios
fig = Figure()
ax = Axis(fig[1,1], xlabel="Year", ylabel="Global Mean Temperature (K)", xticks=1850:50:2100)
lines!(ax, time_future, scenario_1, label="Scenario 1", color=:blue)
lines!(ax, time_future, scenario_2, label="Scenario 2", color=:red)
# lines!(ax, time_future, ssp5, label="SSP5-8.5", color=:green)
lines!(ax, time_future, ssp2, label="SSP2-4.5", color=:purple, linestyle=:dash)
axislegend(ax, position=:lt)
display(fig)

hfile = h5open("new_scenarios.hdf5", "w")
write(hfile, "scenario_1", scenario_1)
write(hfile, "scenario_2", scenario_2)
close(hfile)