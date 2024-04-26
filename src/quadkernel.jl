
function evaluate_covariance(sdf::ParametricSDF{S}, n; asexp_cutoff=2000) where{S}
  T        = eltype(sdf(0.0))
  out      = zeros(T, n)
  krv      = [Kronrod(T) for _ in 1:Threads.nthreads()]
  quadlags = 0:min((n-1), asexp_cutoff)
  chunks   = collect(Iterators.partition(eachindex(quadlags), cld(length(quadlags), length(krv))))
  @sync for (kj, cj) in zip(krv, chunks)
    Threads.@spawn begin
      for k in cj
        out[k] = 2*integrate_interval(0.0, 0.5, kj, w->sdf(w)*cos(2*pi*w*quadlags[k]))
      end
    end
  end
  quadlags[end] == (n-1) && return out
  asexp_ixs  = (length(quadlags)+1):n
  asexp_lags = (quadlags[end]+1):(n-1)
  if length(sdf.asexp_intervals) == 1
    asymptotic_ft_sdf_ff!(view(out, asexp_ixs), sdf, asexp_lags, Val(5))
  else
    lags2pi = asexp_lags.*(2*pi)
    asymptotic_ft!(view(out, asexp_ixs), sdf, lags2pi, sdf.asexp_intervals, Val(5))
  end
  out
end

