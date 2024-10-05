
function evaluate_covariance(sdf::ParametricSDF{S}, n; cutoff=2000, method=:nufft) where{S}
  T        = eltype(sdf(0.0))
  out      = zeros(T, n)
  krv      = [Kronrod(T) for _ in 1:Threads.nthreads()]
  quadlags = 0:min((n-1), cutoff)
  chunks   = collect(Iterators.partition(eachindex(quadlags), cld(length(quadlags), length(krv))))
  @sync for (kj, cj) in zip(krv, chunks)
    Threads.@spawn begin
      for k in cj
        out[k] = 2*integrate_interval(0.0, 0.5, kj, w->sdf(w)*cos(2*pi*w*quadlags[k]))
      end
    end
  end
  quadlags[end] == (n-1) && return out
  ixs     = (length(quadlags)+1):n
  lags    = (quadlags[end]+1):(n-1)
  lags2pi = lags.*(2*pi)
  if method == :asexp
    if length(sdf.asexp_intervals) == 1
      asymptotic_ft_sdf_ff!(view(out, ixs), sdf, lags, Val(5))
    else
      asymptotic_ft!(view(out, ixs), sdf, lags2pi, sdf.asexp_intervals, Val(5))
    end
  elseif method == :nufft
    no_wt_v = map(sdf.asexp_intervals) do (a,b)
      nodes_required = Int(ceil(2*(b-a)*(2*pi*n)))
      (bmad2, bpad2) = ((b-a)/2, (b+a)/2)
      (no, wt) = gausslegendre(nodes_required)
      buf = complex(wt)
      for k in eachindex(no, wt)
        no[k] = bmad2*no[k] + bpad2
        buf[k] = buf[k]*bmad2*sdf(no[k])
      end
      (no, buf)
    end
    no   = reduce(vcat, getindex.(no_wt_v, 1))
    buf  = reduce(vcat, getindex.(no_wt_v, 2))
    itg = real(vec(nufft1d3(no, buf, +1, 1e-15, collect(lags2pi))))
    view(out, ixs) .= itg
  else
    throw(error("The two available methods right now are :asexp and :nufft."))
  end
  out
end

