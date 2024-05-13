
using SpectralEstimators, FFTW, BenchmarkTools, Serialization, GPMaxlik, Printf, StableRNGs, ForwardDiff, JuMP, StandaloneKNITRO

FFTW.set_num_threads(12)

if !(@isdefined _rng)
  const _rng    = StableRNG(12345)
  const n       = 10_000
  const ntrials = 50
  const rnk     = 128
  const true_params = [5.0, 35.0]
end

# Zero derivatives at the origin:
function acf(k, p) 
  (a, b) = p
  num = 2*a*exp(-b/2)*(-b*cospi(k) + exp(b/2)*b + 2*pi*k*sinpi(k))
  den = b^2 + (2*pi*k)^2
  num/den
end
kfn(j,k,p) = acf(abs(j-k), p)

sdf(w, p)     = p[1]*exp(-p[2]*abs(w))
sdf_dp1(w, p) =  sdf(w, (1.0, p[2]))
sdf_dp2(w, p) = -sdf(w, p)*abs(w)
const dsdfv   = (sdf_dp1, sdf_dp2)

# Wrap this in a function so that the GC knows it can throw away the big matrix buffer:
function simulate(acf, true_params)
  buf = GPMaxlik.build_covmat(1:n, acf, true_params)
  buf.M.U'*randn(_rng, n, ntrials)
end

function fit_model(sdfmodel, dsdfv, data; init=[5.0, 15.0], acf=nothing)
  if !isnothing(acf)
    # if an ACF is provided, use the exact log-likelihood.
    nll = p -> begin
      try
        GPMaxlik.gnll_forwarddiff(p, 1:n, data, kfn)
      catch
        NaN
      end
    end
  else
    # otherwise use the approximated one:
    nll = p -> begin
      try
        SpectralEstimators.nll(sdfmodel, p, data; dsdfv=dsdfv)
      catch
        NaN
      end
    end
  end
  knitro_optimize(nll, init; box_lower=fill(1e-8, 2), box_upper=fill(100.0,2), 
                  param_file="sqp_loose.opt")
end


if !isinteractive()

  dummy = let scopetrick = 0

    sims = simulate(kfn, true_params)

    debiased_model = SpectralEstimators.SDFModel(sdf, [0.0], rank=-1)
    whittle_model  = SpectralEstimators.SDFModel(sdf, [0.0], rank=0)
    our_model      = SpectralEstimators.SDFModel(sdf, [0.0], rank=rnk)

    results = map(1:ntrials) do j

      println("\n\n\n*******************************")
      println("*********TRIAL $j/$ntrials************")
      println("*******************************")
      sj  = sims[:,j]
      fsj = fftshift(fft(sj))./sqrt(n)

      # compute estimators:
      println("\n\nTrue nll:")
      tru_time = @elapsed true_mle = fit_model(our_model, dsdfv, sj, acf=acf)
      println("\n\nWhittle nll:")
      whi_time = @elapsed whittle_mle = fit_model(whittle_model, dsdfv, fsj)
      println("\n\nDebiased Whittle nll:")
      deb_time = @elapsed debiased_mle = fit_model(debiased_model, dsdfv, fsj)
      println("\n\nOur nll:")
      our_time = @elapsed our_mle = fit_model(our_model, dsdfv, fsj)

      # compute nll of each estimator:
      tru_nll  =  GPMaxlik.gnll_forwarddiff(true_mle.minimizer,     1:n, sj, kfn)
      whi_nll  =  GPMaxlik.gnll_forwarddiff(whittle_mle.minimizer,  1:n, sj, kfn)
      deb_nll  =  GPMaxlik.gnll_forwarddiff(debiased_mle.minimizer, 1:n, sj, kfn)
      our_nll  =  GPMaxlik.gnll_forwarddiff(our_mle.minimizer,      1:n, sj, kfn)

      (;true_mle, whittle_mle, debiased_mle, our_mle, 
       tru_time,  whi_time,    deb_time,     our_time,
       tru_nll,   whi_nll,     deb_nll,      our_nll)

    end

    serialize("fit_study.jls", results)

  end

end

