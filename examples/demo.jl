
using SpectralEstimators, FFTW

include("model_pairs.jl")

acf(h) = _rough_acf(h)
sdf(w, dummy=nothing) = _rough_sdf(w, dummy)

# some fake data and its (unitary) DFT.
const n   = 2048
const v   = randn(n)
const vft = fftshift(fft(v))./sqrt(n)

# the true negative log-likelihood for vft.
true_nll = let S = [acf(abs(j-k)) for j in 1:n, k in 1:n]
  ft = SpectralEstimators.dftmat(n)
  M  = cholesky!(Hermitian(ft*S*ft'))
  (logdet(M) + sum(abs2, M.U'\vft))/2
end

# the approximated negative log-likelihood for vft.
model    = SpectralEstimators.SDFModel(sdf, [0.0], rank=200)
test_sdf = SpectralEstimators.ParametricSDF(model)
test     = SpectralEstimators.fftcovmat(test_sdf, n)
appx_nll = SpectralEstimators.nll(test, vft)

