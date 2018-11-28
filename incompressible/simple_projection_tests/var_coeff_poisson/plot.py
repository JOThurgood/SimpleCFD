import numpy as np
import matplotlib.pyplot as plt

xc = np.fromfile('x.dat', dtype = np.float64)
yc  = np.fromfile('y.dat', dtype = np.float64)
phi = np.fromfile('phi.dat', dtype = np.float64)
eta = np.fromfile('eta.dat', dtype = np.float64)
f = np.fromfile('f.dat', dtype = np.float64)



phi = phi.reshape(yc.size, xc.size)
eta = eta.reshape(yc.size, xc.size)
f = f.reshape(yc.size, xc.size)


plt.clf()
plt.contourf(xc, yc, phi, 255) 
plt.xlabel('x')
plt.ylabel('y') 
plt.title('phi')
plt.colorbar()
plt.savefig('phi.png')


plt.clf()
plt.contourf(xc, yc, eta, 255) 
plt.xlabel('x')
plt.ylabel('y') 
plt.title('eta')
plt.colorbar()
plt.savefig('eta.png')

plt.clf()
plt.contourf(xc, yc, f, 255) 
plt.xlabel('x')
plt.ylabel('y') 
plt.title('f')
plt.colorbar()
plt.savefig('f.png')

