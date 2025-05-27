module SpectralEstimators

  using LinearAlgebra, FFTW, FINUFFT, FastGaussQuadrature, QuadGK, ForwardDiff, WoodburyMatrices, LowRankApprox

  include("utils.jl")
  include("aquad.jl")
  include("toeplitz.jl")
  include("sdf.jl")
  include("quadkernel.jl")
  include("asymptotic_int.jl")
  include("remainder.jl")
  include("assemble.jl")

end 
