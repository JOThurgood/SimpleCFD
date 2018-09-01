import numpy as np
import matplotlib.pyplot as plt

n = 64+1+2 #ncells + 1 for at boundary + 2 for ghosts

data = np.fromfile('xb.dat', dtype=np.float64, count=n)
xb = data.reshape(n)

data = np.fromfile('a.dat', dtype=np.float64, count=n)
a = data.reshape(n)

plt.plot(xb, a)
if n < 129 :
  plt.plot(xb,a,'bx')
plt.xlabel('x')
plt.ylabel('a(x)')
plt.title('a(x) at t_end')
plt.grid(True)
#plt.show()
plt.savefig("test.png")
