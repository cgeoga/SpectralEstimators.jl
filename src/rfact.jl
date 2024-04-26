
struct PgramRemainder{M}
  sdfv::Vector{Float64}
  toep::M
  adj::Bool
end

Base.adjoint(pm::PgramRemainder{M}) where{M} = PgramRemainder(pm.sdfv, pm.toep, !pm.adj)
Base.size(pm::PgramRemainder{M}, j::Int) where{M} = length(pm.sdfv)

function Base.:*(pm::PgramRemainder{M}, m::Matrix{T}) where{M, T<:Complex}
  m1 = ifftshift(m,  1)
  pl = plan_fft!(m1, 1)
  if pm.adj
    pl*m1
  else
    pl\m1
  end
  m2 = pm.toep*m1
  if pm.adj
    pl\m2
  else
    pl*m2
  end
  fftshift!(m1, m2, 1)
  @inbounds begin
    for k in 1:size(m, 2), j in 1:size(m,1)
      m1[j,k] -= pm.sdfv[j]*m[j,k]
    end
  end
  m1
end

# Simple fixed-rank range finder approximation:
function rangefinder(A, cnd::Int64, p::Int64)
  cnd > size(A, 2) && error("Invalid fixed rank specification: A of size $(size(A,1)), requested rank of $cnd.") 
  T = Matrix(A*complex(randn(size(A, 1), cnd + p)))
  any(isnan, T) && throw(error("NaN detected in A*S."))
  Matrix(qr(T, ColumnNorm()).Q)[:,1:cnd] 
end

function reigen(A, cnd::Int64, p::Int64=8; trim=1e-13, warn_discard=true)
  Q    = rangefinder(A, cnd, p) 
  tmp  = adjoint(Q)*(A*Q)
  tmpe = eigen(Hermitian(tmp))
  U    = Q*tmpe.vectors
  if trim > 0.0 && minimum(abs, tmpe.values) < trim
    ix = findall(x->abs(x) > trim, tmpe.values)
    if warn_discard
      @info "Discarding eigenvalues/vectors corresponding to eigenvalues \
      with magnitude below $trim ($(cnd - length(ix))/$(cnd)). If this keeps \
      happening, consider reducing the fixed rank of your approximation \
      to avoid doing pointless work." maxlog=5
    end
    (U[:,ix], tmpe.values[ix])
  else
    (U, tmpe.values)
  end
end

# Slightly cheaper: gives a randomized low-rank approximation A ≈ C*R', where C
# is an ortogonal matrix but R is not. For simplicity, this function absorbs the
# "rangefinder" function, and ultimately I may move this function and the struct
# into assembly.jl and delete this source file.
function rfact(A, cnd::Int64, p::Int64=8)
  cnd > size(A, 2) && error("Invalid fixed rank specification: A of size $(size(A,1)), requested rank of $cnd.") 
  T = Matrix(A*complex(randn(size(A, 1), cnd + p)))
  any(isnan, T) && throw(error("NaN detected in A*S."))
  C = Matrix(qr(T, ColumnNorm()).Q)[:,1:cnd] 
  AC = adjoint(A)*C
  (C, AC)
end

