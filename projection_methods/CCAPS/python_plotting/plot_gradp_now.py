import numpy as np
import matplotlib.pyplot as plt

xc = np.fromfile('xc.dat', dtype = np.float64)
yc  = np.fromfile('yc.dat', dtype = np.float64)
gradp_x = np.fromfile('gradp_x.dat', dtype = np.float64)
gradp_y = np.fromfile('gradp_y.dat', dtype = np.float64)

gradp_x = gradp_x.reshape(yc.size, xc.size)
gradp_y = gradp_y.reshape(yc.size, xc.size)

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



np.set_printoptions(threshold=np.nan)

print("gradpx")
print(gradp_x)
print("")
print("")
print("")
print("gradpy")
print(gradp_y)
