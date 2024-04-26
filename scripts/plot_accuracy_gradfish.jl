
using LiterallyGnuplot, Serialization, DelimitedFiles

if !(@isdefined ranks)
  const RAW   = deserialize("./output/gradfish_accuracy_scaling_tests.jls")
  const sizes = sort(unique(getindex.(collect(keys(RAW)), 1)))[2:end]
  const cases = sort(unique(getindex.(collect(keys(RAW)), 2)))
  const ranks = sort(unique(getindex.(collect(keys(RAW)), 3)))
end

function slicefield(slfield, sloutput, dv)
  if slfield == :size
    [getfield(RAW[(sj, dv[:case], dv[:rank])], sloutput) for sj in sizes]
  elseif slfield == :case
    [getfield(RAW[(dv[:size], cj, dv[:rank])], sloutput) for cj in cases]
  else
    [getfield(RAW[(dv[:size], dv[:case], rj)], sloutput) for rj in ranks]
  end
end

if !isinteractive()
  for case in (:grad_, :stoch_fish_, :exact_fish_)

    key = Symbol(case, :err)

    smooth_2_errs    = slicefield(:size, key, Dict((:case=>:smooth, :rank=>2)))
    smooth_4_errs    = slicefield(:size, key, Dict((:case=>:smooth, :rank=>4)))
    smooth_32_errs   = slicefield(:size, key, Dict((:case=>:smooth, :rank=>32)))
    smooth_64_errs   = slicefield(:size, key, Dict((:case=>:smooth, :rank=>64)))
    smooth_128_errs  = slicefield(:size, key, Dict((:case=>:smooth, :rank=>128)))

    rough_2_errs     = slicefield(:size, key, Dict((:case=>:rough, :rank=>2)))
    rough_4_errs     = slicefield(:size, key, Dict((:case=>:rough, :rank=>4)))
    rough_32_errs    = slicefield(:size, key, Dict((:case=>:rough, :rank=>32)))
    rough_64_errs    = slicefield(:size, key, Dict((:case=>:rough, :rank=>64)))
    rough_128_errs   = slicefield(:size, key, Dict((:case=>:rough, :rank=>128)))

    writedlm(string("./output/analyticsdf_", case, "errors.csv"), 
             abs.(hcat(sizes, smooth_2_errs, smooth_4_errs, smooth_32_errs, smooth_64_errs, smooth_128_errs)), ',')

    writedlm(string("./output/rough_split_sdf_",  case, "errors.csv"), 
             abs.(hcat(sizes, rough_2_errs, rough_4_errs, rough_32_errs, rough_64_errs, rough_128_errs)), ',')

    key = Symbol(case, :time)

    rough_2_times   = slicefield(:size, key, Dict((:case=>:rough, :rank=>2)))
    rough_4_times   = slicefield(:size, key, Dict((:case=>:rough, :rank=>4)))
    rough_32_times  = slicefield(:size, key, Dict((:case=>:rough, :rank=>32)))
    rough_64_times  = slicefield(:size, key, Dict((:case=>:rough, :rank=>64)))
    rough_128_times = slicefield(:size, key, Dict((:case=>:rough, :rank=>128)))

    # a quick regression line:
    X     = hcat(ones(length(sizes)), [s*log(s) for s in sizes])
    rline = X*(X\rough_64_times)

    writedlm(string("./output/rough_sdf_", case, "times.csv"), 
             hcat(sizes, rline, rough_2_times, rough_4_times, rough_32_times, rough_64_times, rough_128_times), ',')

  end
end
