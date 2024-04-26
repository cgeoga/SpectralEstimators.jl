
using SpectralEstimators, FFTW, BenchmarkTools, Serialization, GPMaxlik

FFTW.set_num_threads(12)

# Analytic:
const _PHI = 0.9
_smooth_sdf(w, dummy) = inv(1 - 2*_PHI*cos(2*pi*w) + _PHI^2)
_smooth_acf(h)  = _PHI^(abs(h))/(1 - _PHI^2)

# Zero derivatives at the origin:
_rough_sdf(w, dummy) = 10*exp(-10*abs(w))
function _rough_acf(k) 
  econ = exp(5)
  num  = pi*k*sinpi(k) - 5*cospi(k) + 5*econ
  den  = econ*((pi*k)^2 + 25)
  10*num/den
end

# Hard because of a huge dynamic rannge:
function _hard_sdf(w, dummy) 
  (s, v, a)  = (5.0, 1.0, 0.0112)
  (s2, a2)   = (s*s, a*a)
  inv(exp((v+1/2)*log(a2 + w^2)))/10000
end
struct HardACF
  v::Vector{Float64}
end
(ha::HardACF)(j::Int) = ha.v[j+1]
const hard_model = SpectralEstimators.SDFModel(_hard_sdf, rank=100) 
const hard_sdf   = SpectralEstimators.ParametricSDF(hard_model)
const _hard_acf  = HardACF(SpectralEstimators.evaluate_covariance(hard_sdf, 100_001))



const nv    = (1_000, 2_000, 4_000, 8_000, 16_000, 24_000, 48_000, 96_000)
const ranks = (-1, 0, 2, 4, 8, 16, 32, 64, 128)
const cases = ((_smooth_sdf, _smooth_acf, :smooth), (_rough_sdf,  _rough_acf,  :naive),
               (_rough_sdf,  _rough_acf,  :rough),  (_hard_sdf,   _hard_acf,   :hard))

results = let scopetrick = 0

  out = Dict{Tuple{Int64, Symbol, Int64}, 
             @NamedTuple{assemble_time::Float64, nll_time::Float64, err::Float64}}()

  for n in nv

    v   = randn(n)
    vft = fftshift(fft(v))./sqrt(n)

    for (_sdf, _acf, _case) in cases

      println("\n\n**Model: $_sdf**")

      refnll = GPMaxlik.gnll_forwarddiff(nothing, 1:n, v, (j,k,p)->_acf(abs(j-k)))

      println("Precompile burn-in...\n\n")
      _model   = SpectralEstimators.SDFModel(_sdf, [0.0], rank=10)
      test_sdf = SpectralEstimators.ParametricSDF(_model)
      test     = SpectralEstimators.fftcovmat(test_sdf, n)
      SpectralEstimators.nll(test, vft);

      for rank in ranks

        model = if _case == :smooth || _case == :naive
          SpectralEstimators.SDFModel(_sdf, rank=rank)
        elseif _case == :rough || _case == :hard
          SpectralEstimators.SDFModel(_sdf, [0.0], rank=rank)
        else
          throw(error("what?"))
        end

        test_sdf = SpectralEstimators.ParametricSDF(model)
        assemble_time = @elapsed test = SpectralEstimators.fftcovmat(test_sdf, n)

        nll_time = @elapsed SpectralEstimators.nll(test, vft);

        cand = SpectralEstimators.nll(test, vft)
        err  = abs(cand - refnll)/abs(refnll)

        println("\nCASE: (n, model, rank) = ($n, $_case, $rank): $(abs(log10(err))) digits")

        out[(n, _case, rank)] = (;assemble_time, nll_time, err)

      end

    end

  end

  out

end

serialize("accuracy_scaling_tests.jls", results)

