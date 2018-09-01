import numpy as np
import matplotlib.pyplot as plt

xc = np.fromfile('xc.dat', dtype=np.float64)
a = np.fromfile('u.dat', dtype=np.float64)

plt.plot(xc, a)
if xc.size < 129 :
  plt.plot(xc,a,'bx')
plt.xlabel('x')
plt.ylabel('u')
plt.title('u(x,t) at t_end')
plt.grid(True)
#plt.show()
plt.savefig("output.png")
