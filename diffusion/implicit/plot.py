import numpy as np
import matplotlib.pyplot as plt

xc = np.fromfile('x.dat', dtype = np.float64)
yc  = np.fromfile('y.dat', dtype = np.float64)
phi = np.fromfile('phi.dat', dtype = np.float64)

phi = phi.reshape(yc.size, xc.size)


plt.clf()
plt.contourf(xc, yc, phi, 255) 
#plt.grid(True)
plt.xlabel('x')
plt.ylabel('y') 
plt.title('phi')
plt.colorbar()
#plt.show()
plt.savefig('phi.png')




