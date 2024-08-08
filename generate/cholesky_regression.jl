using CairoMakie
##
cov1 = cov(modes[:,1,:])
cov13 = cov(modes[:,13,:])
##
ensemble_members = 1:40
size(modes,2)
month = 1
total_modes = 100
ϵ = var(modes[end, month, :]) * 0.1
tmp = [cholesky(cov(Float64.(modes[:,i, ensemble_members]), dims = 2) + ϵ * I ) for i in ProgressBar(month:12:size(modes,2))]
tmp_small = [cholesky(cov(Float64.(modes[1:total_modes,i,ensemble_members]), dims = 2) + 0.01I ) for i in ProgressBar(month:12:size(modes,2))]

# tmp[1].U' * tmp[1].U - cov(Float64.(modes[:,1,:]), dims = 2)
tmp_small[1].U' * tmp_small[1].U - cov(Float64.(modes[1:total_modes,1,ensemble_members]), dims = 2)

var(modes[1000, month, ensemble_members])
tmp_small[1].U


##
order > 5 ? error("Order must be less than 5") : nothing
rX = ones(Float64, 251, 1 + order)
for i in 1:size(modes[1, 1:12:end, :])[1]
    for j in 1:order
        rX[i, 1+j] = Float64(temperature[i])^j
    end
end

mat = zeros(Float64, total_modes, total_modes, order + 1);
for kk in ProgressBar(1:total_modes^2)
    ii = div(kk-1, total_modes) + 1
    jj = mod(kk-1, total_modes) + 1
    data = [tmp[k].U[ii, jj] for k in 1:251]
    j = 251
    regression_coefficients = rX \ data
    model = sum([regression_coefficients[k+1]* (temperature[j])^k for k in 0:order])
    mat[ii, jj, :] .= regression_coefficients
end

j = 1
ii = 1
jj = 1
truth = tmp[j].U[ii, jj]
model = mat[ii, jj, 1] + mat[ii, jj, 2] * temperature[j]
relative_error_percent = abs.((model - truth) / truth) * 100

Umodel = mat[:,:,1] + mat[:,:,2] * temperature[1] 

Umodel' * Umodel - (tmp[j].U' * tmp[j].U)[1:total_modes, 1:total_modes]