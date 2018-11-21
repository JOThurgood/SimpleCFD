import numpy as np
import matplotlib.pyplot as plt

xc = np.fromfile('xc.dat', dtype = np.float64)
yc  = np.fromfile('yc.dat', dtype = np.float64)
gradp_x = np.fromfile('gradp_x.dat', dtype = np.float64)
gradp_y = np.fromfile('gradp_y.dat', dtype = np.float64)
phi = np.fromfile('phi.dat', dtype = np.float64)

gradp_x = gradp_x.reshape(yc.size, xc.size)
gradp_y = gradp_y.reshape(yc.size, xc.size)
phi = phi.reshape(yc.size, xc.size)

plt.clf()
plt.contourf(xc, yc, gradp_x, 255,cmap=plt.cm.RdBu_r) 
#plt.imshow(divu, cmap=plt.cm.RdBu_r)
plt.grid(True)
plt.xlabel('x')
plt.ylabel('y') 
plt.title('gradp_x')
plt.colorbar()
#plt.show()
plt.savefig('gradp_x.png')

plt.clf()
plt.contourf(xc, yc, gradp_y, 255,cmap=plt.cm.RdBu_r) 
#plt.imshow(divu, cmap=plt.cm.RdBu_r)
plt.grid(True)
plt.xlabel('x')
plt.ylabel('y') 
plt.title('gradp_y')
plt.colorbar()
#plt.show()
plt.savefig('gradp_y.png')

plt.clf()
plt.contourf(xc, yc, phi, 255,cmap=plt.cm.RdBu_r) 
#plt.imshow(divu, cmap=plt.cm.RdBu_r)
plt.grid(True)
plt.xlabel('x')
plt.ylabel('y') 
plt.title('phi')
plt.colorbar()
#plt.show()
plt.savefig('phi.png')

plt.clf()
plt.plot(yc,phi[:,xc.size//2])
plt.scatter(yc,phi[:,xc.size//2],color='red')
plt.axvline(x=0.0, color='r', linestyle='-')
plt.axvline(x=1.0, color='r', linestyle='-')
plt.xlabel('y')
plt.ylabel('phi(x=middle)')
plt.savefig('phi_cut.png')

np.set_printoptions(threshold=np.nan)

print("gradpx")
print(gradp_x)
print("")
print("")
print("")
print("gradpy")
print(gradp_y)
print("")
print("")
print("")
print("phi")
print(phi)

