import numpy as np
import matplotlib.pyplot as plt

n = 64+3
data = np.fromfile('xb.dat', dtype=np.float64, count=n)
xb = data.reshape(n)
data = np.fromfile('a.dat', dtype=np.float64, count=n)
a = data.reshape(n)

plt.plot(xb, a)
#plt.show()
plt.savefig("test.png")
