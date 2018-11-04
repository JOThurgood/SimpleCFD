import numpy as np
import matplotlib.pyplot as plt

xc = np.fromfile('xc.dat', dtype = np.float64)
yc  = np.fromfile('yc.dat', dtype = np.float64)
u = np.fromfile('u.dat', dtype = np.float64)
v = np.fromfile('v.dat', dtype = np.float64)
divu = np.fromfile('divu.dat', dtype = np.float64)
curlu = np.fromfile('curlu.dat', dtype = np.float64)
#phi = np.fromfile('phi.dat', dtype = np.float64)

u = u.reshape(yc.size, xc.size)
v = v.reshape(yc.size, xc.size)
divu = divu.reshape(yc.size, xc.size)
curlu = curlu.reshape(yc.size, xc.size)
#phi = phi.reshape(yc.size, xc.size)


plt.clf()
plt.contourf(xc, yc, u, 255, cmap=plt.cm.RdBu_r)
#plt.grid(True)
plt.xlabel('x')
plt.ylabel('y') 
plt.title('u')
plt.colorbar()
#plt.show()
plt.savefig('u.png')




plt.clf()
plt.contourf(xc, yc, v, 255, cmap=plt.cm.RdBu_r)
#plt.grid(True)
plt.xlabel('x')
plt.ylabel('y') 
plt.title('v')
plt.colorbar()
#plt.show()
plt.savefig('v.png')

plt.clf()
plt.contourf(xc, yc, divu, 255, cmap=plt.cm.RdBu_r)
#plt.grid(Trve)
plt.xlabel('x')
plt.ylabel('y') 
plt.title('divu')
plt.colorbar()
#plt.show()
plt.savefig('divu.png')

plt.clf()
plt.contourf(xc, yc, curlu, 255, cmap=plt.cm.RdBu_r)
#plt.grid(Trve)
plt.xlabel('x')
plt.ylabel('y') 
plt.title('curl u')
plt.colorbar()
#plt.show()
plt.savefig('curlu.png')

#plt.clf()
#plt.contourf(xc, yc, phi, 255) 
##plt.grid(Trve)
#plt.xlabel('x')
#plt.ylabel('y') 
#plt.title('phi')
#plt.colorbar()
##plt.show()
#plt.savefig('phi.png')



