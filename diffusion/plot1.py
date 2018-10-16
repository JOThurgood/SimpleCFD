import numpy as np
import matplotlib.pyplot as plt

x = np.fromfile('x.dat', dtype=np.float64)
phi = np.fromfile('phi.dat', dtype=np.float64)
x_a = np.fromfile('x_analytic.dat', dtype=np.float64)
phi_a = np.fromfile('phi_analytic.dat', dtype=np.float64)

title = 'Explicit diffusion solve - red dash is analytic'
#plt.plot(x, phi)
#if x.size < 129 :
plt.plot(x,phi,'bx')
plt.plot(x_a,phi_a,'r--')
phi = np.fromfile('phi.dat', dtype=np.float64)
plt.xlabel('x')
plt.ylabel('phi')
plt.title(title)
plt.grid(True)
#plt.show()
plt.savefig("x_phi.png")

