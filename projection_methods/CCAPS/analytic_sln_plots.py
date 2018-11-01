import numpy as np
import matplotlib.pyplot as plt

xc = np.fromfile('xc.dat', dtype = np.float64)
yc  = np.fromfile('yc.dat', dtype = np.float64)
u = np.fromfile('u.dat', dtype = np.float64)
v = np.fromfile('v.dat', dtype = np.float64)
utest = np.fromfile('utest.dat', dtype = np.float64)
vtest = np.fromfile('vtest.dat', dtype = np.float64)

u = u.reshape(yc.size, xc.size)
v = v.reshape(yc.size, xc.size)
utest = utest.reshape(yc.size, xc.size)
vtest = vtest.reshape(yc.size, xc.size)


plt.clf()
plt.contourf(xc, yc, u, 255) 
#plt.grid(True)
plt.xlabel('x')
plt.ylabel('y') 
plt.title('u (simulated)')
#plt.show()
plt.savefig('u_simulated.png')

plt.clf()
plt.contourf(xc, yc, v, 255) 
#plt.grid(True)
plt.xlabel('x')
plt.ylabel('y') 
plt.title('v (simulated)')
#plt.show()
plt.savefig('v_simulated.png')

plt.clf()
plt.contourf(xc, yc, utest, 255) 
#plt.grid(True)
plt.xlabel('x')
plt.ylabel('y') 
plt.title('u (analytic)')
#plt.show()
plt.savefig('u_analytic.png')

plt.clf()
plt.contourf(xc, yc, vtest, 255) 
#plt.grid(True)
plt.xlabel('x')
plt.ylabel('y') 
plt.title('v (analytic)')
#plt.show()
plt.savefig('v_analytic.png')


