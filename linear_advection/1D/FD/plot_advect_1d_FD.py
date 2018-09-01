import numpy as np
import matplotlib.pyplot as plt


xb = np.fromfile('xb.dat', dtype=np.float64)
a = np.fromfile('a.dat', dtype=np.float64)

plt.plot(xb, a)
if xb.size < 130 :
  plt.plot(xb,a,'bx')
plt.xlabel('x')
plt.ylabel('a(x)')
plt.title('a(x) at t_end')
plt.grid(True)
#plt.show()
plt.savefig("output.png")
