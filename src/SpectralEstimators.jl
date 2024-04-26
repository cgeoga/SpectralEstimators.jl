module SpectralEstimators

  using LinearAlgebra, FFTW, QuadGK, ForwardDiff, WoodburyMatrices

  include("utils.jl")
  include("aquad.jl")
  include("toeplitz.jl")
  include("sdf.jl")
  include("quadkernel.jl")
  include("asymptotic_int.jl")
  include("rfact.jl")
  include("assemble.jl")

end 
