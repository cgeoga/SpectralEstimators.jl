# SpectralEstimators.jl

This package provides tools for scalably and accurate evaluating the
log-likelihood of a Gaussian time series using its spectral density. It is the
software companion to [Fast Machine-Precision Spectral Likelihoods for
Stationary Time Series](https://arxiv.org/abs/2404.16583). 

# Demonstration

Here is a very quick demonstration that should get you off the ground and
collecting all the digits quickly:
```julia
using LinearAlgebra, SpectralEstimators

# Write some parametric family of spectral densities you want to fit. Here is a
# neat model that I don't mention in the paper, but is fun to play with because
# it is smooth at the endpoints but potentially rough at the origin, so the
# parameter alpha really directly controls how quickly the ACF decays.
function sdf(w, params)
    (phi, rho, alpha) = params
    exp(-(abs(sinpi(w))^alpha)/rho)
end

# Now for the interface, you create an object called "SpectralModel". You can
# choose how to evaluate the tail sequence of your ACF: :asexp for the
# asymptotic expansion, that will be faster, or :nufft, which will be less
# error-prone to functions coded in such a way that AD doesn't give the right
# answer.
model = SpectralEstimators.SDFModel(sdf,     # your SDF, see above
                                    [0.0],   # rough points of the SDF
                                    rank=72, # fixed rank for the Whittle correction
                                    kernel_tail_method=:nufft)

# Now we're ready to build our fast and machine-exact approximation for
# F*Sigma*F^H, as described in the paper. In this example, we build a 3k x 3k
# matrix.
parametric_sdf = SpectralEstimators.ParametricSDF(model, (5.0, 0.1, 1.25))
fftcov         = SpectralEstimators.fftcovmat(parametric_sdf, 3_000)

# Enjoy! Take a look at the source to see all the methods. Here is a simple case
# of the negative log-likelihood (nll):
@show SpectralEstimators.nll(fftcov, rand(3000))
```

