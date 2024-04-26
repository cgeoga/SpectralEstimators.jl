
struct PDHermitianToeplitz{F,T<:Complex}
  c::Vector{T}
  c_ext_ft::Vector{T}
  buf::Matrix{T}
  plan::F
end

function PDHermitianToeplitz(v; nbufcols=1)
  cv    = copy(complex(v))
  c_ext = vcat(cv, conj.(reverse(cv)[1:(end-1)]))
  buf   = zeros(complex(eltype(v)), 2*length(v)-1, nbufcols)
  plan  = plan_fft!(buf, 1)
  fft!(c_ext)
  PDHermitianToeplitz(cv, c_ext, buf, plan)
end

Base.size(toep::PDHermitianToeplitz{F,T}, j) where{F,T} = length(toep.c)

function LinearAlgebra.mul!(buf::Matrix{X}, 
                            toep::PDHermitianToeplitz{F,T}, 
                            x::Matrix{X}) where{F,T,X}
  if size(toep.buf, 2) < size(x, 2)
    throw(error("This routine assumes that you have pre-allocated enough columns to do all matvecs at once. Please rebuild your PDHermitianToeplitz matrix with the nbufcols kwarg >= $(size(x,2)) for this call."))
  end
  if size(x, 1) != size(toep, 1)
    throw(error("Input argument and matrix size do not agree."))
  end
  (sx1, sx2) = size(x)
  fill!(toep.buf, zero(T))
  copyto!(view(toep.buf, 1:sx1, 1:sx2), x)
  toep.plan*toep.buf
  @inbounds begin
    for j in 1:size(toep.buf,1)
      cj = toep.c_ext_ft[j]
      for k in 1:size(x,2)
        toep.buf[j,k] *= cj
      end
    end
  end
  toep.plan\toep.buf
  @inbounds begin
    for j in 1:size(buf, 1)
      for k in 1:size(buf, 2)
        buf[j,k] = X(toep.buf[j,k])
      end
    end
  end
  buf
end

Base.:*(toep::PDHermitianToeplitz{F,T}, x) where{F,T} = mul!(copy(x), toep, x)

