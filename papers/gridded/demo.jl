
using SpectralEstimators

# Here's a neat SDF: you choose the smoothness at the origin, and it is smooth
# at the periodic endpoints. Constrain α∈[1,2], ϕ>0, and ρ>0. Even something
# this simple doesn't have a closed form ACF as far as I know.
function sdf(w, params)
  (phi, rho, alpha) = params
  exp(-(abs(sin(pi*w))^alpha)/rho)
end

# Create a SDFModel, where you encode the points at which the SDF is nonsmooth
# (but is directonally smooth) and what rank you want to use for approximating
# the Whittle correction.
model = SpectralEstimators.SDFModel(sdf, [0.0], rank=72, kernel_tail_method=:nufft)

# Now, for a specific set of parameters, you can create a ParametricSDF object:
parametric_sdf = SpectralEstimators.ParametricSDF(model, (5.0, 0.1, 1.25))

#=
# Finally, you can perform the fast assemble of FFT * CovMatrix(above) * FFT^H:
fftcov = SpectralEstimators.fftcovmat(parametric_sdf, 3_000)

# Now a quick test of accuracy:
dftn = SpectralEstimators.dftmat(3_000)
acf  = SpectralEstimators.evaluate_covariance(parametric_sdf, 3_000, method=:nufft)
naive_matrix = dftn * [acf[abs(j-k)+1] for j in 1:3_000, k in 1:3_000] * dftn'

@show opnorm(Matrix(fftcov) - naive_matrix)/opnorm(naive_matrix)
=#
