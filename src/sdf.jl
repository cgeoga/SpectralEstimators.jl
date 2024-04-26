
struct SDFModel{S}
  sdf::S
  asexp_intervals::Vector{Tuple{Float64, Float64}}
  rank::Int64 # rank of the remainder term
end

function SDFModel(sdf::S, roughpoints=Float64[]; rank=50) where{S} 
  if in(-1/2, roughpoints) || in(1/2, roughpoints)
    throw(error("This code currently does not handle roughness at the endpoints correctly, sorry. Please use a different model."))
  end
  rpts = vcat(-1/2, roughpoints, 1/2)
  asexp_intervals = map(z->(nextfloat(z[1]), prevfloat(z[2])), zip(rpts, rpts[2:end]))
  SDFModel(sdf, asexp_intervals, rank)
end

struct ParametricSDF{S,P}
  sdf::S
  params::P
  asexp_intervals::Vector{Tuple{Float64, Float64}}
  rank::Int64
end
(psdf::ParametricSDF{S,P})(w) where{S,P} = psdf.sdf(w, psdf.params)

function ParametricSDF(sdfm::SDFModel{S}, params::P=nothing) where{S,P} 
  ParametricSDF(sdfm.sdf, params, sdfm.asexp_intervals, sdfm.rank)
end

struct ComponentFunction{J,F,P,X}
  f::F
  p::P
  x::X
end

struct ComponentDerivative{J,F}
  f::F
end

@generated function splice(x::NTuple{N,T}, newval::G, idx::Val{J}) where{J,N,T,G}
  quote
    ($([:(x[$j]) for j in 1:(J-1)]...), newval, $([:(x[$j]) for j in (J+1):N]...))
  end
end

function (k::ComponentFunction{J,F,P,X})(pj) where{J,F,P,X}
  k.f(k.x, splice(k.p, pj, Val(J)))
end

function (dk::ComponentDerivative{J,F})(x::X, p::P) where{J,F,X,P}
  cf = ComponentFunction{J,F,P,X}(dk.f, p, x)
  ForwardDiff.derivative(cf, p[J])
end

@generated function component_derivatives(f::F, ::Val{N}) where{F,N}
  quote
    Base.Cartesian.@ntuple $N j->ComponentDerivative{j,F}(f)
  end
end

