using NCDatasets, HDF5, Statistics, ProgressBars

include("utils.jl")

# --------- CONFIG ---------
# save_directory = "PLEASE/SET/YOUR/SAVE/PATH/HERE/"
# data_directory = "PLEASE/SET/YOUR/DATA/PATH/HERE/"
save_directory = "/net/fs06/d3/mgeo/GaussianEarthData/"
data_directory = "/net/fs06/d3/mgeo/CMIP6/interim/"


output_file = joinpath(save_directory, "ground_truth_bundle.hdf5")

scenarios = ["historical", "ssp119", "ssp245", "ssp585"]
ensemble_members = Dict("tas" => 45, "hurs" => 29)

# --------------------------

function default_chunk(dims)
    if length(dims) == 1
        return (min(dims[1], 4096),)
    else
        return ntuple(i -> min(dims[i], 32), length(dims))
    end
end

function ensure_group(h5, path)
    parts = split(path, "/")
    if length(parts) == 1
        return h5, parts[1]
    end
    grp = h5
    for p in parts[1:end-1]
        if haskey(grp, p)
            grp = grp[p]
        else
            grp = g_create(grp, p)
        end
    end
    return grp, parts[end]
end

function write_compressed(h5, path, data; compress = 3)
    if data isa BitArray
        data = Array{Bool}(data)
    end
    grp, name = ensure_group(h5, path)
    dims = size(data)
    chunk = default_chunk(dims)
    dset = d_create(grp, name, datatype(eltype(data)), dataspace(dims); chunk = chunk, compress = compress)
    write(dset, data)
    close(dset)
end

function list_member_files(scenario, variable; data_directory, n)
    base = joinpath(data_directory, scenario, variable)
    files = sort(readdir(base))
    files = files[1:min(n, length(files))]
    return [joinpath(base, f) for f in files]
end

function compute_mean_std_monthly(scenario, variable; data_directory, n_members)
    files = list_member_files(scenario, variable; data_directory, n=n_members)
    ds0 = Dataset(files[1])
    field0 = ds0[variable]
    nx, ny, nt = size(field0)
    years = nt ÷ 12
    close(ds0)

    mean_monthly = zeros(Float32, nx, ny, years, 12)
    m2_monthly = zeros(Float32, nx, ny, years, 12)

    count = 0
    for file in ProgressBar(files)
        ds = Dataset(file)
        field = Float32.(ds[variable][:,:,:])
        close(ds)

        count += 1
        for month in 1:12
            slice = field[:, :, month:12:end]
            delta = slice .- mean_monthly[:, :, :, month]
            mean_monthly[:, :, :, month] .+= delta ./ count
            m2_monthly[:, :, :, month] .+= delta .* (slice .- mean_monthly[:, :, :, month])
        end
    end

    std_monthly = sqrt.(m2_monthly ./ (count - 1))
    annual_mean = mean(mean_monthly, dims=4)[:, :, :, 1]

    return mean_monthly, std_monthly, annual_mean
end

function load_metric_and_mask(save_directory)
    metric = nothing
    mask = nothing

    basis_file = joinpath(save_directory, "tas_basis.hdf5")
    if ispath(basis_file)
        hfile = h5open(basis_file, "r")
        metric = read(hfile["metric"])
        close(hfile)
    end

    mask_file = joinpath(save_directory, "land_sea_mask.hdf5")
    if ispath(mask_file)
        hfile = h5open(mask_file, "r")
        mask = read(hfile["mask"]) .> 0.5
        close(hfile)
    end

    return metric, mask
end

function hemisphere_means(field_slice, metric; upper=true)
    if upper
        w = metric[:, 49:end]
        f = field_slice[:, 49:end, :]
    else
        w = metric[:, 1:48]
        f = field_slice[:, 1:48, :]
    end
    vals = sum(f .* w, dims = (1, 2)) * 2
    return vec(vals)
end

function land_sea_means(field_slice, metric, mask; land=true)
    if land
        m = .!mask
    else
        m = mask
    end
    w = metric .* m
    denom = sum(w)
    vals = sum(field_slice .* w, dims = (1, 2)) ./ denom
    return vec(vals)
end

function collect_location_samples!(samples, field_slice, lon_idx, lat_idx)
    append!(samples, vec(field_slice[lon_idx, lat_idx, :]))
end

function collect_hist_samples(
    scenario, variable, month, year_inds;
    data_directory, n_members, metric=nothing, mask=nothing,
    locations=Tuple{Int,Int}[], hemispheres=Symbol[], land_ocean=Symbol[]
)
    files = list_member_files(scenario, variable; data_directory, n=n_members)
    samples = Dict{String, Vector{Float32}}()

    for (lon_idx, lat_idx) in locations
        key = "loc_$(lon_idx)_$(lat_idx)"
        samples[key] = Float32[]
    end
    for h in hemispheres
        samples[string(h)] = Float32[]
    end
    for lo in land_ocean
        samples[string(lo)] = Float32[]
    end

    for file in ProgressBar(files)
        ds = Dataset(file)
        field = Float32.(ds[variable][:,:,:])
        close(ds)

        field_month = field[:, :, month:12:end]
        field_slice = field_month[:, :, year_inds]

        for (lon_idx, lat_idx) in locations
            key = "loc_$(lon_idx)_$(lat_idx)"
            collect_location_samples!(samples[key], field_slice, lon_idx, lat_idx)
        end

        for h in hemispheres
            if metric === nothing
                error("metric is required for hemisphere means")
            end
            vals = hemisphere_means(field_slice, metric; upper = (h == :upper))
            append!(samples[string(h)], Float32.(vals))
        end

        for lo in land_ocean
            if metric === nothing || mask === nothing
                error("metric and mask are required for land/ocean means")
            end
            vals = land_sea_means(field_slice, metric, mask; land = (lo == :land))
            append!(samples[string(lo)], Float32.(vals))
        end
    end

    return samples
end

function year_window(center, window, n)
    lo = max(1, center - window)
    hi = min(n, center + window)
    return lo:hi
end

function compute_year_inds_historical(temps)
    year_inds = 1:(argmin(temps) - 1)
    accept_lower = temps .> minimum(temps[year_inds])
    accept_upper = temps .< maximum(temps[year_inds])
    accept = accept_lower .& accept_upper
    return collect(eachindex(temps))[accept]
end

@info "Writing ground-truth bundle to $output_file"
h5 = h5open(output_file, "w")

# Save grid if available
metric, mask = load_metric_and_mask(save_directory)
if metric !== nothing
    write_compressed(h5, "grid/metric", metric)
end
if mask !== nothing
    write_compressed(h5, "grid/mask", mask)
end

# --- Aggregate statistics for tas (mean/std monthly + annual mean) ---
@info "Computing ensemble mean/std for tas"
for scenario in scenarios
    mean_monthly, std_monthly, annual_mean = compute_mean_std_monthly(
        scenario, "tas";
        data_directory, n_members = ensemble_members["tas"]
    )
    write_compressed(h5, "tas/mean_monthly/$scenario", mean_monthly)
    write_compressed(h5, "tas/std_monthly/$scenario", std_monthly)
    write_compressed(h5, "tas/annual_mean/$scenario", annual_mean)
end

# --- Histogram samples for figures 2/7/8 ---
@info "Collecting histogram samples for figures 2/7/8"

# Figure 2 (tas only, ssp119)
temps_ssp119 = regression_variable("ssp119"; directory = save_directory)
n_years_ssp119 = length(temps_ssp119)
window = 2
indmin = 6
indmax = 81
yr_low = year_window(indmin, window, n_years_ssp119)
yr_high = year_window(indmax, window, n_years_ssp119)

for month in (1, 7)
    samples_low = collect_hist_samples(
        "ssp119", "tas", month, yr_low;
        data_directory, n_members = ensemble_members["tas"],
        metric = metric,
        locations = [(157, 46), (19, 87)],
        hemispheres = [:upper, :lower]
    )
    samples_high = collect_hist_samples(
        "ssp119", "tas", month, yr_high;
        data_directory, n_members = ensemble_members["tas"],
        metric = metric,
        locations = [(157, 46), (19, 87)],
        hemispheres = [:upper, :lower]
    )
    for (k, v) in samples_low
        write_compressed(h5, "samples/figure_2/month_$month/$k/low", v)
    end
    for (k, v) in samples_high
        write_compressed(h5, "samples/figure_2/month_$month/$k/high", v)
    end
end

# Figure 7 & 8 (tas + hurs, ssp245/ssp585 vs historical)
temps_hist = regression_variable("historical"; directory = save_directory)
temps_ssp245 = regression_variable("ssp245"; directory = save_directory)
temps_ssp585 = regression_variable("ssp585"; directory = save_directory)
n_years_hist = length(temps_hist)
n_years_ssp245 = length(temps_ssp245)
n_years_ssp585 = length(temps_ssp585)
indmin = 100
indmax = 43
yr_hist = year_window(indmin, window, n_years_hist)
yr_ssp245 = year_window(indmax, window, n_years_ssp245)
yr_ssp585 = year_window(indmax, window, n_years_ssp585)

for month in (1, 7)
    for field in ("tas", "hurs")
        # figure_7: land/ocean averages
        s_hist = collect_hist_samples(
            "historical", field, month, yr_hist;
            data_directory, n_members = ensemble_members[field],
            metric = metric, mask = mask,
            land_ocean = [:land, :ocean]
        )
        s_ssp245 = collect_hist_samples(
            "ssp245", field, month, yr_ssp245;
            data_directory, n_members = ensemble_members[field],
            metric = metric, mask = mask,
            land_ocean = [:land, :ocean]
        )
        s_ssp585 = collect_hist_samples(
            "ssp585", field, month, yr_ssp585;
            data_directory, n_members = ensemble_members[field],
            metric = metric, mask = mask,
            land_ocean = [:land, :ocean]
        )
        for (k, v) in s_hist
            write_compressed(h5, "samples/figure_7/$field/month_$month/$k/historical", v)
        end
        for (k, v) in s_ssp245
            write_compressed(h5, "samples/figure_7/$field/month_$month/$k/ssp245", v)
        end
        for (k, v) in s_ssp585
            write_compressed(h5, "samples/figure_7/$field/month_$month/$k/ssp585", v)
        end

        # figure_8: point locations
        s_hist = collect_hist_samples(
            "historical", field, month, yr_hist;
            data_directory, n_members = ensemble_members[field],
            locations = [(157, 46), (19, 87)]
        )
        s_ssp245 = collect_hist_samples(
            "ssp245", field, month, yr_ssp245;
            data_directory, n_members = ensemble_members[field],
            locations = [(157, 46), (19, 87)]
        )
        s_ssp585 = collect_hist_samples(
            "ssp585", field, month, yr_ssp585;
            data_directory, n_members = ensemble_members[field],
            locations = [(157, 46), (19, 87)]
        )
        for (k, v) in s_hist
            write_compressed(h5, "samples/figure_8/$field/month_$month/$k/historical", v)
        end
        for (k, v) in s_ssp245
            write_compressed(h5, "samples/figure_8/$field/month_$month/$k/ssp245", v)
        end
        for (k, v) in s_ssp585
            write_compressed(h5, "samples/figure_8/$field/month_$month/$k/ssp585", v)
        end
    end
end

# --- Figure 6 and A1: historical tas (month=1) ---
@info "Collecting figure_6 and figure_A1 data"
temps_hist = regression_variable("historical"; directory = save_directory)
year_inds = compute_year_inds_historical(temps_hist)
month = 1

# Figure 6 samples and zonal stats (mean over longitude)
lat_indices = [1, 24, 48, 72, 96]
lat_samples = Dict(i => Float32[] for i in lat_indices)
lat_sum = nothing
lat_sumsq = nothing
lat_count = 0

files = list_member_files("historical", "tas"; data_directory, n=ensemble_members["tas"])
for file in ProgressBar(files)
    global lat_sum, lat_sumsq, lat_count
    ds = Dataset(file)
    field = Float32.(ds["tas"][:,:,:])
    close(ds)
    field_month = field[:, :, month:12:end]
    field_slice = field_month[:, :, year_inds] # (lon, lat, years)

    # mean over longitude -> (lat, years)
    zonal = mean(field_slice, dims = 1)[1, :, :]  # lat x years

    if lat_sum === nothing
        lat_sum = zeros(Float64, size(zonal, 1))
        lat_sumsq = zeros(Float64, size(zonal, 1))
    end
    lat_sum .+= sum(zonal, dims = 2)[:]
    lat_sumsq .+= sum(zonal .^ 2, dims = 2)[:]
    lat_count += size(zonal, 2)

    for i in lat_indices
        append!(lat_samples[i], Float32.(vec(zonal[i, :])))
    end
end

lat_mean = Float32.(lat_sum ./ lat_count)
lat_std = Float32.(sqrt.(lat_sumsq ./ lat_count .- lat_mean .^ 2))
write_compressed(h5, "stats/figure_6/latitude_mean", lat_mean)
write_compressed(h5, "stats/figure_6/latitude_std", lat_std)
write_compressed(h5, "stats/figure_6/year_inds", year_inds)
for i in lat_indices
    write_compressed(h5, "samples/figure_6/lat_$i", lat_samples[i])
end

# Figure A1: kurtosis/skewness and sample points
sum1 = nothing
sum2 = nothing
sum3 = nothing
sum4 = nothing
count = 0
for file in ProgressBar(files)
    global sum1, sum2, sum3, sum4, count
    ds = Dataset(file)
    field = Float32.(ds["tas"][:,:,:])
    close(ds)
    field_month = field[:, :, month:12:end]
    field_slice = field_month[:, :, year_inds] # (lon, lat, years)
    if sum1 === nothing
        sum1 = zeros(Float64, size(field_slice, 1), size(field_slice, 2))
        sum2 = zeros(Float64, size(field_slice, 1), size(field_slice, 2))
        sum3 = zeros(Float64, size(field_slice, 1), size(field_slice, 2))
        sum4 = zeros(Float64, size(field_slice, 1), size(field_slice, 2))
    end
    sum1 .+= sum(field_slice, dims = 3)[:, :, 1]
    sum2 .+= sum(field_slice .^ 2, dims = 3)[:, :, 1]
    sum3 .+= sum(field_slice .^ 3, dims = 3)[:, :, 1]
    sum4 .+= sum(field_slice .^ 4, dims = 3)[:, :, 1]
    count += size(field_slice, 3)
end

μ = sum1 ./ count
σ2 = sum2 ./ count .- μ .^ 2
σ = sqrt.(σ2)
skew = (sum3 ./ count .- 3 .* μ .* σ2 .- μ .^ 3) ./ (σ .^ 3)
kurt = (sum4 ./ count .- 4 .* μ .* (sum3 ./ count) .+ 6 .* μ .^ 2 .* (sum2 ./ count) .- 3 .* μ .^ 4) ./ (σ .^ 4) .- 3

write_compressed(h5, "stats/figure_A1/skewness", Float32.(skew))
write_compressed(h5, "stats/figure_A1/kurtosis", Float32.(kurt))

vec_skew = vec(skew)
vec_kurt = vec(kurt)
gaussian_max = argmin(abs.(vec_kurt) .+ abs.(vec_skew))
skewness_min = argmin(vec_skew)
skewness_max = argmax(vec_skew)
kurtosis_min = argmin(vec_kurt)
kurtosis_max = argmax(vec_kurt)

indices = [gaussian_max, kurtosis_min, skewness_min, kurtosis_max, skewness_max]
write_compressed(h5, "stats/figure_A1/selected_indices", indices)

# store samples for those points
for (j, idx) in enumerate(indices)
    lon_idx = (idx - 1) % size(skew, 1) + 1
    lat_idx = (idx - 1) ÷ size(skew, 1) + 1
    samples = collect_hist_samples(
        "historical", "tas", month, year_inds;
        data_directory, n_members = ensemble_members["tas"],
        locations = [(lon_idx, lat_idx)]
    )
    key = "loc_$(lon_idx)_$(lat_idx)"
    write_compressed(h5, "samples/figure_A1/point_$j", samples[key])
end

close(h5)
@info "Done"
