
function dftmat(n)
  ir  = vcat(div(n,2):(n-1), 0:(div(n,2)-1))
  [exp(-2*pi*im*j*k/n)/sqrt(n) for j in ir, k in 0:(n-1)]
end

function pgram_covmat(kfn::F, n) where{F}
  dft = dftmat(n)
  S   = [kfn(abs(j-k)) for j in 0:(n-1), k in 0:(n-1)]
  Hermitian(dft*S*dft')
end

function nll(M::AbstractMatrix, v::Vector)
  Mf = cholesky(Hermitian(M))
  (logdet(Mf) + sum(abs2, Mf.U'\v))/2
end

# This is piracy, should just open a PR.
LinearAlgebra.tr(wm::Woodbury) = tr(wm.A) + tr((wm.V*wm.U)*wm.C)

function debiased_whittle_diag(kv)
  n   = length(kv)
  tri = [1 - t/n for t in 0:(length(kv)-1)]
  fftshift(2*real(fft(tri.*kv)) .- kv[1])
end

