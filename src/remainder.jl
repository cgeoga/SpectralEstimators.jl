
struct PgramRemainder{M}
  sdfv::Vector{Float64}
  toep::M
end

LinearAlgebra.ishermitian(pm::PgramRemainder{M}) where{M} = true
LinearAlgebra.issymmetric(pm::PgramRemainder{M}) where{M} = false
Base.eltype(pm::PgramRemainder{M}) where{M}       = ComplexF64
Base.size(pm::PgramRemainder{M}) where{M}         = (length(pm.sdfv), length(pm.sdfv))
Base.size(pm::PgramRemainder{M}, j::Int) where{M} = size(pm)[j]

function LinearAlgebra.mul!(buf, pm::PgramRemainder{M}, m) where{M}
  m1 = ifftshift(m,  1)
  pl = plan_fft!(m1, 1)
  pl\m1
  m2 = pm.toep*m1
  pl*m2
  fftshift!(m1, m2, 1)
  @inbounds begin
    for k in 1:size(m, 2), j in 1:size(m,1)
      m1[j,k] -= pm.sdfv[j]*m[j,k]
    end
  end
  buf .= m1
end

