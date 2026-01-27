import Base.Float32
Float32(a::Missing) = Float32(0.0)
ensemble_members = 30

scenario_directories = readdir(data_directory)
scenario_directories = ["historical", "ssp585"] #overwrite

function get_all_fields(scenario; variable = "tas", data_directory = data_directory)
    fields = Array{Float32, 3}[]
    current_path = joinpath(data_directory, scenario)
    local_current_path = joinpath(current_path, variable)
    file_names = readdir(local_current_path)
    if length(file_names) > 0 # sometimes the directory is empty
        file_name = file_names[1] # pick the first file for obtaining varaibles
        file_path = joinpath(local_current_path, file_name)
        for (i, file_name) in ProgressBar(enumerate(file_names[1:ensemble_members]))
            file_path = joinpath(local_current_path, file_name)
            ds = Dataset(file_path)
            field = Float32.(ds[variable][:,:,:])
            push!(fields, field)
        end
    end
    return fields
end

function restructure_fields(fields)
    time1 = size(fields[1])[end]
    time2 = size(fields[end])[end]
    total_time = time1 + time2
    total_time_month = total_time ÷ 12
    all_together = zeros(Float32, size(fields[1])[1], size(fields[1])[2], total_time_month, 12, ensemble_members);

    monthtime1 = time1 ÷ 12
    monthtime2 = time2 ÷ 12
    for ω in ProgressBar(1:ensemble_members)
        for month in 1:12
            all_together[:,:, 1:monthtime1, month, ω] .= fields[ω][:, :, month:12:end]
        end
    end
    for ω in ProgressBar(1:ensemble_members)
        for month in 1:12
            all_together[:,:, monthtime1+1:end, month, ω] .= fields[ω+ensemble_members][:, :, month:12:end]
        end
    end
    return all_together
end

function pattern_scaling_fit(all_together)
    historical = regression_variable("historical")
    ssp585 = regression_variable("ssp585")
    x = cat(historical, ssp585, dims = 1)
    xs = repeat(x, 30)
    X = ones(size(xs, 1), 2)
    X[:, 2] .= xs

    linear_fit = zeros(Float32, size(all_together)[1], size(all_together)[2], 2, 12)
    for month in ProgressBar(1:12)
        for i in 1:size(all_together)[1]
            for j in 1:size(all_together)[2]
                M, N = size(all_together[i,j, :, month, :])
                Y = reshape(all_together[i,j, :, month, :], (M*N)) 
                linear_fit[i, j, :, month] .= X \ Y
            end
        end
    end

    linear_fit_yearly = zeros(Float32, size(all_together)[1], size(all_together)[2], 2)
    for i in 1:size(all_together)[1]
        for j in 1:size(all_together)[2]
            M, N = size(all_together[i,j, :, 1, :])
            Y = reshape(mean(all_together[i,j, :, :, :], dims = 2)[:, 1, :], (M*N)) 
            linear_fit_yearly[i, j, :] .= X \ Y
        end
    end
    return linear_fit, linear_fit_yearly
end

@info "Computing Projection for all variables and all scenarios"
variables = ["tas", "hurs"]
for variable in ProgressBar(variables)
    fields = Array{Float32, 3}[]
    for scenario in ProgressBar(scenario_directories)
        scenario_fields = get_all_fields(scenario; variable)
        fields = vcat(fields, scenario_fields)
    end
    all_together = restructure_fields(fields)
    # pattern scaling
    linear_fit, linear_fit_yearly = pattern_scaling_fit(all_together)
    # save
    filename  = variable * "_pattern_scaling.hdf5"
    hfile = h5open(save_directory * filename, "w")
    hfile["linear fit yearly"] = linear_fit_yearly
    hfile["linear fit monthly"] = linear_fit
    close(hfile)
end