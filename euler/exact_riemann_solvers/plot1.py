import numpy as np
import matplotlib.pyplot as plt

x = np.fromfile('x.dat', dtype=np.float64)
rho = np.fromfile('rho.dat', dtype=np.float64)
p = np.fromfile('p.dat', dtype=np.float64)

plt.plot(x, rho)
if x.size < 129 :
  plt.plot(x,rho,'bx')
plt.xlabel('x')
plt.ylabel('rho')
plt.title('title')
plt.grid(True)
#plt.show()
plt.savefig("rho.png")

plt.clf()
plt.plot(x, p)
if x.size < 129 :
  plt.plot(x,p,'bx')
plt.xlabel('x')
plt.ylabel('p')
plt.title('title')
plt.grid(True)
#plt.show()
plt.savefig("p.png")

