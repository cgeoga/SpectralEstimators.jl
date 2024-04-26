
using LiterallyGnuplot, Serialization, DelimitedFiles

if !(@isdefined ranks)
  const RAW   = deserialize("./output/accuracy_scaling_tests.jls")
  const sizes = sort(unique(getindex.(collect(keys(RAW)), 1)))
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

smooth_d_errs    = slicefield(:size, :err, Dict((:case=>:smooth, :rank=>-1)))
smooth_0_errs    = slicefield(:size, :err, Dict((:case=>:smooth, :rank=>0)))
smooth_2_errs    = slicefield(:size, :err, Dict((:case=>:smooth, :rank=>2)))
smooth_32_errs   = slicefield(:size, :err, Dict((:case=>:smooth, :rank=>32)))
smooth_64_errs   = slicefield(:size, :err, Dict((:case=>:smooth, :rank=>64)))
smooth_128_errs  = slicefield(:size, :err, Dict((:case=>:smooth, :rank=>128)))

naive_d_errs     = slicefield(:size, :err, Dict((:case=>:naive, :rank=>-1)))
naive_0_errs     = slicefield(:size, :err, Dict((:case=>:naive, :rank=>0)))
naive_2_errs     = slicefield(:size, :err, Dict((:case=>:naive, :rank=>2)))
naive_32_errs    = slicefield(:size, :err, Dict((:case=>:naive, :rank=>32)))
naive_64_errs    = slicefield(:size, :err, Dict((:case=>:naive, :rank=>64)))
naive_128_errs   = slicefield(:size, :err, Dict((:case=>:naive, :rank=>128)))

rough_d_errs     = slicefield(:size, :err, Dict((:case=>:rough, :rank=>-1)))
rough_0_errs     = slicefield(:size, :err, Dict((:case=>:rough, :rank=>0)))
rough_2_errs     = slicefield(:size, :err, Dict((:case=>:rough, :rank=>2)))
rough_32_errs    = slicefield(:size, :err, Dict((:case=>:rough, :rank=>32)))
rough_64_errs    = slicefield(:size, :err, Dict((:case=>:rough, :rank=>64)))
rough_128_errs   = slicefield(:size, :err, Dict((:case=>:rough, :rank=>128)))

writedlm("./output/analyticsdf_errors.csv", 
         hcat(sizes, smooth_d_errs, smooth_0_errs, smooth_2_errs, smooth_32_errs, smooth_64_errs, smooth_128_errs), ',')

writedlm("./output/rough_naive_sdf_errors.csv", 
         hcat(sizes, naive_d_errs, naive_0_errs, naive_2_errs, naive_32_errs, naive_64_errs, naive_128_errs), ',')

writedlm("./output/rough_split_sdf_errors.csv", 
         hcat(sizes, rough_d_errs, rough_0_errs, rough_2_errs, rough_32_errs, rough_64_errs, rough_128_errs), ',')




naive_d_times   = slicefield(:size, :assemble_time, Dict((:case=>:naive, :rank=>-1)))
naive_0_times   = slicefield(:size, :assemble_time, Dict((:case=>:naive, :rank=>0)))
naive_2_times   = slicefield(:size, :assemble_time, Dict((:case=>:naive, :rank=>2)))
naive_32_times  = slicefield(:size, :assemble_time, Dict((:case=>:naive, :rank=>32)))
naive_64_times  = slicefield(:size, :assemble_time, Dict((:case=>:naive, :rank=>64)))
naive_128_times = slicefield(:size, :assemble_time, Dict((:case=>:naive, :rank=>128)))

# a quick regression line:
X     = hcat(ones(length(sizes)), [s*log(s) for s in sizes])
rline = X*(X\naive_64_times)

writedlm("./output/rough_naive_sdf_times.csv", 
         hcat(sizes, rline, naive_d_times, naive_0_times, naive_2_times, naive_32_times, naive_64_times, naive_128_times), ',')

