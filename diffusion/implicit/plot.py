import numpy as np
import matplotlib.pyplot as plt

xc = np.fromfile('x.dat', dtype = np.float64)
yc  = np.fromfile('y.dat', dtype = np.float64)
phi = np.fromfile('phi.dat', dtype = np.float64)
phi_an = np.fromfile('phi_an.dat', dtype = np.float64)

phi = phi.reshape(yc.size, xc.size)
phi_an = phi_an.reshape(yc.size, xc.size)


plt.clf()
plt.contourf(xc, yc, phi, 255) 
#plt.grid(True)
plt.xlabel('x')
plt.ylabel('y') 
plt.title('phi (simulated)')
plt.colorbar()
#plt.show()
plt.savefig('phi_simulated.png')

plt.clf()
plt.contourf(xc, yc, phi_an, 255) 
#plt.grid(True)
plt.xlabel('x')
plt.ylabel('y') 
plt.title('phi (analytic)')
plt.colorbar()
#plt.show()
plt.savefig('phi_analytic.png')



