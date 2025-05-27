
abstract type QuadRule end

struct Kronrod{S} <: QuadRule
  segbuf::S
end

function Kronrod(::Type{T}=Float64, n::Int64=5000) where{T}
  segbuf = QuadGK.alloc_segbuf(Float64, T, Float64; size=n)
  Kronrod(segbuf)
end

function integrate_interval(a, b, rule::Kronrod, f::F; tol=1e-14) where{F}
  # TODO (cg 2023/11/05 15:59): order not picked precisely here. Could optimize.
  quadgk(f, a, b; order=32, atol=tol, segbuf=rule.segbuf)[1]
end

