

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

