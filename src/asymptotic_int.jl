
fwd_ho(x, ::Val{0}) = x
fwd_ho(x, ::Val{1}) = ForwardDiff.Dual{typeof(fwd_ho)}(x, (1.0,))
fwd_ho(x, ::Val{N}) where{N} = ForwardDiff.Dual{typeof(fwd_ho)}(fwd_ho(x, Val(N-1)), (1.0,))

get_partial(x, ::Val{0}) = x
get_partial(x, ::Val{1}) = x.partials[1]
get_partial(x, ::Val{N}) where{N} = get_partial(x.partials[1], Val(N-1))

get_value(x, ::Val{0})  = x
get_value(x, ::Val{1}) = x.value
get_value(x, ::Val{N}) where{N} = get_value(x.value, Val(N-1))

@generated function hoderivatives(fn::F, x, ::Val{O}) where{F,O}
  quote
    dx  = fwd_ho(x, Val($O));
    fdx = fn(dx)
    d0  = get_partial(fdx, Val($O))
    tmp = fdx
    dxs = Base.Cartesian.@ntuple $O j-> begin
      layer = get_partial(tmp, Val($O - j + 1))
      tmp   = tmp.value
      layer
    end
    out = reverse((dxs..., tmp))
    any(isnan, out) && throw(error("Derivatives up to $O at x=$x have at least one NaN."))
    out
  end
end

# This computes an integral of \int_{a}^b e^{i*w*x} f(x) dx with an asymptotic
# expansion of order O.
@generated function asymptotic_ft(fn::F, w::Float64, a, b, ::Val{O}) where{F,O}
  quote
    (out, k) = (0.0, 0)
    dfdo_a = hoderivatives(fn, a, Val($O))
    dfdo_b = hoderivatives(fn, b, Val($O))
    Base.Cartesian.@nexprs $O j-> begin
      out += (dfdo_b[j]*cis(b*w) - dfdo_a[j]*cis(a*w))/((-im*w)^(k+1))
      k   += 1
    end
    -out
  end
end

# specific function for computing \int_{-1/2}^{1/2} fn(w) e^{2 pi i k w} dw,
# which has a lot of nice simplifications and cancellations.
@generated function asymptotic_ft_sdf_ff!(buf::AbstractVector{Float64}, fn::F, 
                                          kv::AbstractVector{Int64}, ::Val{O}) where{F,O}
  quote
    fill!(buf, 0.0)
    dfdo_a = hoderivatives(fn, -1/2, Val($O))
    dfdo_b = Base.Cartesian.@ntuple $O j -> isodd(j) ? dfdo_a[j] : -dfdo_a[j]
    for l in eachindex(buf, kv)
      k = kv[l]
      w = 2*pi*k
      nimw = -im*w
      nimw_pow = nimw
      tmp = zero(ComplexF64)
      Base.Cartesian.@nexprs $O j -> begin
        tmp += (dfdo_b[j] - dfdo_a[j])/nimw_pow
        nimw_pow *= nimw
      end
      tmp = ifelse(isodd(k), tmp, -tmp)
      buf[l] = real(tmp)
    end
  end
end

@generated function asymptotic_ft!(buf::AbstractVector{T}, fn::F, wv::AbstractVector, 
                                   abv::AbstractVector, ::Val{O}) where{T,F,O}
  valfun = (T <: Complex) ? identity : real
  quote
    fill!(buf, zero(T))
    for (a, b) in abv
      dfdo_a = hoderivatives(fn, a, Val($O))
      dfdo_b = hoderivatives(fn, b, Val($O))
      for l in eachindex(buf, wv)
        w = wv[l]
        (tmp, k) = (zero(ComplexF64), 0)
        Base.Cartesian.@nexprs $O j -> begin
          tmp += (dfdo_b[j]*cis(b*w) - dfdo_a[j]*cis(a*w))/((-im*w)^(k+1))
          k   += 1
        end
        buf[l] -= $valfun(tmp)
      end
    end
  end
end

