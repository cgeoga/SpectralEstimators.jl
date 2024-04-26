
struct FFTCovMat{W}
  woodbury::W
  method::Symbol
end

Base.Matrix(ftc::FFTCovMat{W}) where{W}  = Matrix(ftc.woodbury)
Base.adjoint(ftc::FFTCovMat{W}) where{W} = ftc

# sometimes A isn't invertible, but I just need to multiply the matrix so who
# cares.
struct WeakWoodbury
  A::Diagonal{Float64, Vector{Float64}}
  U::Matrix{ComplexF64}
  C::Diagonal{Float64, Vector{Float64}}
  V::Adjoint{ComplexF64, Matrix{ComplexF64}}
end

Base.:*(ww::WeakWoodbury, m)       = ww.A*m + ww.U*(ww.C*(ww.V*m))
LinearAlgebra.tr(ww::WeakWoodbury) = tr(ww.A) + tr((ww.V*ww.U)*ww.C)
Base.Matrix(ww::WeakWoodbury)      = ww.A + ww.U*(ww.C*ww.V)
Base.adjoint(ww::WeakWoodbury)     = WeakWoodbury(ww.A, ww.V', ww.C, ww.U')
Base.size(ww::WeakWoodbury, j)     = size(ww.A, j)

function _fftcovmat(sdf::ParametricSDF{S,P}, n, ::Type{T}; 
                    sketchmat=0, method=:reigen) where{S,P,T}
  kv    = evaluate_covariance(sdf, n)
  tm    = PDHermitianToeplitz(kv, nbufcols=sdf.rank+10)
  wgrid = range(-0.5, 0.5, length=n+1)[1:n]
  sdfv  = complex(sdf.(wgrid))
  if sdf.rank == -1 # "debiased" whittle
    D = Diagonal(complex(debiased_whittle_diag(kv)))
    wood = T(D, zeros(n,1), [1.0;;], zeros(n,1)')
  elseif iszero(sdf.rank) # standard Whittle
    D = Diagonal(sdfv)
    wood = T(D, zeros(n,1), [1.0;;], zeros(n,1)')
  else
    lmap  = PgramRemainder(real(sdfv), tm, false)
    if method == :fast
      (C,R) = rfact(lmap, sdf.rank)
      wood = T(Diagonal(sdfv), C, Diagonal(ones(ComplexF64, sdf.rank)), R')
    elseif method == :reigen
      (U,L) = reigen(lmap, sdf.rank)
      wood = T(Diagonal(sdfv), U, Diagonal(L), U')
    else
      throw(error("Method options are :fast (unstructured low-rank approx), or :reigen (eigen)."))
    end
  end
  FFTCovMat(wood, method)
end

# TODO (cg 2024/03/12 16:05): check type stability of this hack.
function fftcovmat(sdf::ParametricSDF{S,P}, n; method=:reigen) where{S,P}
  _fftcovmat(sdf, n, Woodbury; method=method)
end

function weakfftcovmat(sdf::ParametricSDF{S,P}, n; method=:reigen) where{S,P}
  _fftcovmat(sdf, n, WeakWoodbury; method=method)
end

function LinearAlgebra.logdet(ftcov::FFTCovMat{W}) where{W} 
  (ld, sgn_ld) = logabsdet(ftcov.woodbury)
  real(sgn_ld) < 0.0 && @warn "logdet of ::FFTCovMat was < 0, proceed with caution..."
  real(ld)
end

quadform(ftcov::FFTCovMat{W}, v::Vector) where{W} = real(dot(v, ftcov.woodbury\v))

nll(ftcov::FFTCovMat{W}, v::Vector) where{W} = (logdet(ftcov) + quadform(ftcov, v))/2

function nll(ftcov::FFTCovMat{W}, vm::Matrix) where{W}
  ld = logdet(ftcov)
  sv = ftcov.woodbury\vm
  m  = size(vm, 2)
  (m*ld + sum(j->dot(view(vm, :, j), view(sv, :, j)), 1:m))/2
end

function nll(sdfm::SDFModel{S}, p::Vector{Float64}, data) where{S}
  ftcov = fftcovmat(ParametricSDF(sdfm, p), size(data, 1))
  nll(ftcov, data) 
end

# TODO (cg 2024/02/23 12:58): Find a way to switch to the more stable calculation.
function ambikasaranX(U, K)
  L  = cholesky(Hermitian(U'U)).U'
  M  = cholesky(Hermitian(I + L'*K*L)).U'
  iL = inv(L)
  iL'*(M - I)*iL
end

struct FFTCovMatSqrt{W}
  sqrtD::Diagonal{Float64, Vector{Float64}}
  woodbury::W
  adj::Bool
end

# Returns an object that implements * and \.
function symfactor(ftcov::FFTCovMat{W}) where{W}
  @assert ftcov.method == :reigen "This only works if a randomized partial eigenfactorization was used to assemble the approximation."
  wood  = ftcov.woodbury
  sqrtD = sqrt(real(wood.A))
  modU  = sqrtD\wood.U
  X     = ambikasaranX(modU, wood.C)
  FFTCovMatSqrt(sqrtD, Woodbury(Diagonal(ones(size(sqrtD, 1))), modU, X, modU'), false)
end

Base.Matrix(ftcovw::FFTCovMatSqrt{W}) where{W}  = ftcovw.sqrtD*Matrix(ftcovw.woodbury)

# Annoying: if I make this return an Adjoint type, then it is a <:
# AbstractMatrix and all of my custom methods get ignored. That sucks.
Base.adjoint(ftcovw::FFTCovMatSqrt{W}) where{W} = FFTCovMatSqrt(ftcovw.sqrtD, ftcovw.woodbury, !ftcovw.adj)

function Base.:*(ftcovsqrt::FFTCovMatSqrt{W}, m) where{W}
  sqrtD = ftcovsqrt.sqrtD
  wood  = ftcovsqrt.woodbury
  ftcovsqrt.adj ? wood'*(sqrtD*m) : sqrtD*(wood*m)
end

function Base.:\(ftcovsqrt::FFTCovMatSqrt{W}, m) where{W}
  sqrtD = ftcovsqrt.sqrtD
  wood  = ftcovsqrt.woodbury
  ftcovsqrt.adj ? sqrtD\(wood'\m) : wood\(sqrtD\m)
end

struct WoodburyProductLRPart{WT1,WT2}
  adj::Bool
  W1::WT1
  W2::WT2
end

Base.size(wp::WoodburyProductLRPart, j) = size(wp.W1, j)

Base.adjoint(wp::WoodburyProductLRPart) = WoodburyProductLRPart(!wp.adj, wp.W1, wp.W2)

# TODO (cg 2024/03/11 08:31): this could of course be optimized.
function Base.:*(wp::WoodburyProductLRPart, m)
  if wp.adj
    wp.W2'*(wp.W1'*m) - wp.W2.A'*(wp.W1.A'*m)
  else
    wp.W1*(wp.W2*m) - wp.W1.A*(wp.W2.A*m)
  end
end

function woodbury_rmult(w1, w2)
  wprodlr = WoodburyProductLRPart(false, w1, w2)
  (U, V) = rfact(wprodlr, min(size(w1.U, 2)*3, size(w1,1)))
  WeakWoodbury(w1.A*w2.A, U, Diagonal{ComplexF64}(I(size(U,2)))  , V')
end

function nll_gradient(sdfm::SDFModel{S}, p::Vector{Float64}, data; dsdfv=nothing) where{S}
  pt = tuple(p...)
  _nll_gradient(sdfm, pt, data, dsdfv)
end

function _nll_gradient(sdfm::SDFModel{S}, p::NTuple{N,Float64}, data, dsdfv) where{S,N}
  # assemble the primal approximation:
  ftcov  = fftcovmat(ParametricSDF(sdfm, p), size(data, 1))
  iftcov = inv(ftcov.woodbury)
  # loop over the length of the parameters and compute gradient entries:
  if isnothing(dsdfv)
    dsdfv = component_derivatives(sdfm.sdf, Val(N))
  end
  dpsdfv = ntuple(j->ParametricSDF(dsdfv[j], p, sdfm.asexp_intervals, sdfm.rank), N)
  out = map(dpsdfv) do dpsdfj
    # assemble the derivative matrix:
    dftcovj = weakfftcovmat(dpsdfj, size(data, 1)) 
    # compute the trace term:
    tr_prod = tr(woodbury_rmult(iftcov, dftcovj.woodbury))
    # compute the quadratic form term:
    t1 = ftcov.woodbury\data
    t2 = dftcovj.woodbury*t1
    t3 = ftcov.woodbury\t2
    qf = dot(data, t3)
    # return the trace and quadform values:
    real((tr_prod - qf)/2)
  end
  collect(out) # returning a vector is probably best for compat.
end

function nll_fish(sdfm::SDFModel{S}, p::Vector{Float64}, data; 
                  dsdfv=nothing, saa=nothing) where{S}
  pt = tuple(p...)
  _nll_fish(sdfm, pt, data, dsdfv, saa)
end

function _nll_fish(sdfm::SDFModel{S}, p::NTuple{N,Float64}, data, dsdfv, saa) where{S,N}
  # assemble the primal approximation:
  ftcov = fftcovmat(ParametricSDF(sdfm, p), size(data, 1))
  ftcov_fact = symfactor(ftcov)
  # pre-solve the SAA vectors:
  if isnothing(saa)
    saa = rand((-1.0, 1.0), length(data), 72)
  end
  nsaa = size(saa, 2)
  solved_saa = adjoint(ftcov_fact)\saa
  # compute the SDF derivatives if they weren't provided:
  if isnothing(dsdfv)
    dsdfv = component_derivatives(sdfm.sdf, Val(N))
  end
  dpsdfv = ntuple(j->ParametricSDF(dsdfv[j], p, sdfm.asexp_intervals, sdfm.rank), N)
  # compute the half-solved, derivative-applied, half-solved SAA vectors:
  applied_saa = map(dpsdfv) do dpsdfj
    # assemble the derivative matrix:
    dftcovj = weakfftcovmat(dpsdfj, size(data, 1)) 
    # apply to the half-solved SAA vectors:
    ftcov_fact\(dftcovj.woodbury*solved_saa)
  end
  # assemble the matrix:
  h = zeros(ComplexF64, N, N)
  # first pass: computed the symmetrized estimator.
  for j in 1:N
    saaj = applied_saa[j]
    # compute the diagonal value:
    h[j,j] = 0.5*sum(abs2, saaj)/nsaa
    # compute the off-diagonal value:
    for k in 1:(j-1)
      saak   = applied_saa[k]
      h[j,k] = sum(z->abs2(z[1]+z[2]), zip(saaj, saak))/(4*nsaa)
      h[k,j] = conj(h[j,k])
    end
  end
  # second pass: correct the off-diagonal values:
  for j in 1:N
    for k in 1:N
      if j != k
        h[j,k] -= 0.5*(h[j,j] + h[k,k])
      end
    end
  end
  real(h .* size(data, 2))
end

function _nll_fish_exact(sdfm::SDFModel{S}, p::NTuple{N,Float64}, data, dsdfv) where{S,N}
  # assemble the primal approximation:
  ftcov  = fftcovmat(ParametricSDF(sdfm, p), size(data, 1))
  iftcov = inv(ftcov.woodbury)
  # compute the SDF derivatives if they weren't provided:
  dpsdfv = ntuple(j->ParametricSDF(dsdfv[j], p, sdfm.asexp_intervals, sdfm.rank), N)
  # compute the half-solved, derivative-applied, half-solved SAA vectors:
  Sjv = map(dpsdfv) do fj
    woodbury_rmult(iftcov, weakfftcovmat(fj, size(data, 1)).woodbury)
  end
  h = zeros(ComplexF64, N, N)
  for j in 1:N
    for k in 1:j
      h[j,k] = tr(woodbury_rmult(Sjv[j], Sjv[k]))/2
      h[k,j] = conj(h[j,k])
    end
  end

  real(h .* size(data, 2))
end

