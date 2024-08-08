using CairoMakie
# Only fit a few modes of the covariance matrix
##
cov1 = cov(modes[:,1,:])
cov13 = cov(modes[:,13,:])
##
regression_indices = [1:139..., 141:251...]
ensemble_members = 1:40
size(modes,2)
month = 1
total_modes = 100
ϵ = var(modes[end, month, :]) * 0.1
tmp = [cholesky(cor(Float64.(modes[:,i, ensemble_members]), dims = 2) + sqrt(eps(1.0)) * I) for i in ProgressBar(month:12:size(modes,2))]
tmp_small = [cholesky(cor(Float64.(modes[1:total_modes,i,ensemble_members]), dims = 2) + sqrt(eps(1.0)) * I ) for i in ProgressBar(month:12:size(modes,2))]
tmp_std = [std(Float64.(modes[:,i, ensemble_members]), dims = 2) for i in ProgressBar(month:12:size(modes,2))]

# tmp[1].U' * tmp[1].U - cov(Float64.(modes[:,1,:]), dims = 2)
tmp_small[1].U' * tmp_small[1].U - cor(Float64.(modes[1:total_modes,1,ensemble_members]), dims = 2)

var(modes[1000, month, ensemble_members])
tmp_small[1].U


##
order = 1
order > 5 ? error("Order must be less than 5") : nothing
rX = ones(Float64, 251, 1 + order)
for i in 1:size(modes[1, 1:12:end, :])[1]
    for j in 1:order
        rX[i, 1+j] = Float64(temperature[i])^j
    end
end

total_modes = 200 # size(tmp[1].U)[1] ÷ 2
mat = zeros(Float64, total_modes, total_modes, order + 1);
for kk in ProgressBar(1:total_modes^2)
    ii = div(kk-1, total_modes) + 1
    jj = mod(kk-1, total_modes) + 1
    if ii ≤ jj
        data = [tmp[k].U[ii, jj] for k in 1:251]
        j = 251
        regression_coefficients = rX \ data
        model = sum([regression_coefficients[k+1]* (temperature[j])^k for k in 0:order])
        mat[ii, jj, :] .= regression_coefficients
    end
end

##
j = 1
ii = 1
jj = 1
truth = tmp[j].U[ii, jj]
model = sum([mat[ii,jj, k+1]* (temperature[j])^k for k in 0:order])
relative_error_percent = abs.((model - truth) / truth) * 100

Umodel = sum([mat[:,:, k+1]* (temperature[j])^k for k in 0:order])
Σ_model = Umodel' * Umodel
Σ_truth = tmp[j].U' * tmp[j].U

Σ_model[ii, jj]
Σ_truth[ii, jj]

data = [tmp[k].U[ii, jj] for k in 1:251]

fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, 1:251, data, color = :blue)
display(fig)

##
tmp_cor = [(cor(Float64.(modes[:,i, ensemble_members]), dims = 2) + sqrt(eps(1.0)) * I) for i in ProgressBar(month:12:size(modes,2))]
tmp_var = [(cov(Float64.(modes[:,i, ensemble_members]), dims = 2) + sqrt(eps(1.0)) * I) for i in ProgressBar(month:12:size(modes,2))]
##
order = 0
order > 5 ? error("Order must be less than 5") : nothing
rX = ones(Float64, 251, 1 + order)
for i in 1:size(modes[1, 1:12:end, :])[1]
    for j in 1:order
        rX[i, 1+j] = Float64(temperature[i])^j
    end
end


total_modes = 200 # size(tmp[1].U)[1] ÷ 2
mat = zeros(Float64, total_modes, total_modes, order + 1);
for kk in ProgressBar(1:total_modes^2)
    ii = div(kk-1, total_modes) + 1
    jj = mod(kk-1, total_modes) + 1
    data = [tmp_cor[k][ii,jj] for k in 1:251]
    j = 251
    regression_coefficients = rX \ data
    model = sum([regression_coefficients[k+1]* (temperature[j])^k for k in 0:order])
    mat[ii, jj, :] .= regression_coefficients
end

##
j = 250
ii = 1
jj = 2
model = sum([mat[:, :, k+1]* (temperature[j])^k for k in 0:order])
truth = tmp_cor[j]
##
min_eigenvalues = [minimum(eigvals(sum([mat[:, :, k+1]* (temperature[j])^k for k in 0:order]))) for j in ProgressBar(1:251)]
max_eigenvalues = [maximum(eigvals(sum([mat[:, :, k+1]* (temperature[j])^k for k in 0:order]))) for j in ProgressBar(1:251)]
##
fig = Figure()
ax = Axis(fig[1, 1])
data = [tmp_cor[k][ii,jj] for k in 1:251]
model = [sum([mat[ii, jj, k+1]* (temperature[j])^k for k in 0:order]) for j in 1:251]
scatter!(ax, temperature, data, color = :blue)
lines!(ax, temperature, model, color = :red)
display(fig)

##
fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, 1:251, temperature, color = :blue)
display(fig)