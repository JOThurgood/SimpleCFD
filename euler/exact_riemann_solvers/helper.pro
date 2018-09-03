gamma = 1.4d0 
p0 = 0.55d0 

fl = 2d0 / (gamma-1d0) * SQRT(gamma) * $ 
  (p0^((gamma-1d0)/2d0/gamma) - 1d0)
print,fl

ar = 2d0 / (gamma+1d0) / 0.125d0 
br = (gamma-1d0)/(gamma+1d0) * 0.1
fr = (p0 - 0.1d0 ) * sqrt( ar / (p0 + br))
print,fr 

print,fl+fr

end

