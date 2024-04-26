
using SpectralEstimators, FFTW, BenchmarkTools, Serialization, GPMaxlik, ForwardDiff, Printf
using ForwardDiff.DiffResults
import GPMaxlik: gnll_forwarddiff

GPMaxlik.SETTINGS.STOCH_GRAD = false

FFTW.set_num_threads(12)

# Analytic:
const _PHI = 0.9
_smooth_sdf(w, p) = p[1]*inv(1 - 2*p[2]*cos(2*pi*w) + p[2]^2)
_smooth_acf(k, p)  = p[1]*p[2]^(abs(k))/(1 - p[2]^2)
_smooth_sdf_dp1(w, p) = _smooth_sdf(w, (1.0, p[2]))
_smooth_sdf_dp2(w, p) = p[1]*2*(cos(2*pi*w) - p[2])/abs2(2*cos(2*pi*w)*p[2] - p[2]^2 - 1)

# Zero derivatives at the origin:
_rough_sdf(w, p) = p[1]*exp(-p[2]*abs(w))
function _rough_acf(k, p) 
  (a, b) = p
  num = 2*a*exp(-b/2)*(-b*cospi(k) + exp(b/2)*b + 2*pi*k*sinpi(k))
  den = b^2 + (2*pi*k)^2
  num/den
end
_rough_sdf_dp1(w, p) = _rough_sdf(w, (1.0, p[2]))
_rough_sdf_dp2(w, p) = -_rough_sdf(w, p)*abs(w)

const nv    = (500, 1_000, 2_000, 4_000, 8_000, 16_000, 24_000, 48_000) # will throw out 500 since it is just precompile burn-in
const ranks = (2, 4, 8, 16, 32, 64, 128)
const cases = ((_smooth_sdf, _smooth_acf, :smooth, (_smooth_sdf_dp1, _smooth_sdf_dp2), [1.0, 0.9]), 
               (_rough_sdf,  _rough_acf,  :rough,  (_rough_sdf_dp1,  _rough_sdf_dp2),  [5.0, 10.0]))


results = let scopetrick = 0

  out = Dict{Tuple{Int64, Symbol, Int64}, 
             @NamedTuple{grad_time::Float64, grad_err::Float64, 
                         stoch_fish_time::Float64, stoch_fish_err::Float64,
                         exact_fish_time::Float64, exact_fish_err::Float64}}()

  for n in nv

    saa = rand((-1.0, 1.0), n, 72)

    # not just straight white noise: some vaguely dependent process that _sort
    # of_ looks like an OU process. I'm lazily throwing away the "burn-in" in
    # the beginning.
    v   = (cumsum(randn(2*n))./sqrt.(1:(2*n)))[(n+1):end]
    vft = fftshift(fft(v))./sqrt(n)

    for (_sdf, _acf, _case, _dsdf_dpv, testp) in cases

      println("\n\n**Model: $_sdf**")

      # Get the exact gradient and efish:
      _ker      = (j, k, p) -> _acf(abs(j-k), p)
      _ker_dp1  = (x, y, p) -> ForwardDiff.derivative(p1 -> _ker(x, y, (p1, p[2])), p[1])
      _ker_dp2  = (x, y, p) -> ForwardDiff.derivative(p2 -> _ker(x, y, (p[1], p2)), p[2])
      refres    = gnll(1:n, v, _ker, (_ker_dp1, _ker_dp2), testp; 
                       nll=false, grad=true, fish=true, saa=nothing)
      refgrad   = refres.grad
      reffish   = refres.fish

      for rank in ranks

        model = if _case == :smooth || _case == :naive
          SpectralEstimators.SDFModel(_sdf, rank=rank)
        elseif _case == :rough || _case == :hard
          SpectralEstimators.SDFModel(_sdf, [0.0], rank=rank)
        else
          throw(error("what?"))
        end

        testp_t = tuple(testp...)

        grad_time       = @elapsed       grad = SpectralEstimators.nll_gradient(model, testp, vft, dsdfv=_dsdf_dpv)
        stoch_fish_time = @elapsed stoch_fish = SpectralEstimators.nll_fish(model, testp, vft, dsdfv=_dsdf_dpv, saa=saa)
        exact_fish_time = @elapsed exact_fish = SpectralEstimators._nll_fish_exact(model, testp_t, vft, _dsdf_dpv)

        grad_err        = maximum(abs.(grad - refgrad)./abs.(refgrad))
        stoch_fish_err  = opnorm(reffish - stoch_fish)/opnorm(reffish)
        exact_fish_err  = opnorm(reffish - exact_fish)/opnorm(reffish)

        @printf "\nGradient     (%i, %s, %i): %2.5e error, %3.5e seconds\n" n _case rank grad_err grad_time
        @printf "Stoch Fisher (%i, %s, %i): %2.5e error, %3.5e seconds\n" n _case rank stoch_fish_err stoch_fish_time
        @printf "Exact Fisher (%i, %s, %i): %2.5e error, %3.5e seconds\n" n _case rank exact_fish_err exact_fish_time


        out[(n, _case, rank)] = (;grad_time, grad_err, stoch_fish_time, stoch_fish_err, exact_fish_time, exact_fish_err)

      end

    end

  end

  out

end

serialize("gradfish_accuracy_scaling_tests.jls", results)

