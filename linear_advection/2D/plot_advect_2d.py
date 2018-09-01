import numpy as np
import matplotlib.pyplot as plt

xc = np.fromfile('xc.dat', dtype = np.float64)
yc  = np.fromfile('yc.dat', dtype = np.float64)
a = np.fromfile('a.dat', dtype = np.float64)
a = a.reshape(yc.size, xc.size)


plt.clf()
plt.contour(xc, yc, a, 10) 
plt.grid(True)
plt.xlabel('x')
plt.ylabel('y') 
plt.title('a(x,y,t_end)')
#plt.show()
plt.savefig('output_contour.png')

plt.clf()
plt.contourf(xc, yc, a, 255, 
  cmap =  plt.cm.bone) 
plt.colorbar()
plt.grid(False)
plt.xlabel('x')
plt.ylabel('y') 
plt.title('a(x,y,t_end)')
#plt.show()
plt.savefig('output_contourf.png')

plt.clf()
plt.imshow(a)
plt.xlabel('x')
plt.ylabel('y') 
plt.title('a(x,y,t_end)')
plt.colorbar()
plt.savefig('output_pixels.png')

