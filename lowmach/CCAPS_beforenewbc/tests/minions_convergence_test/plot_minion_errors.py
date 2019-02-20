import numpy as np
import matplotlib.pyplot as plt

nx = np.array([32., 64., 128.])
L2 = np.array([9.8070375464557186E-003, 1.7188907566484436E-003, 3.5289939673991170E-004])


ax = plt.subplot(111)
ax.set_xscale('log')
ax.set_yscale('log')

plt.scatter(nx, L2, marker="x")
plt.plot(nx,L2[0]*(nx[0]/nx)**2,"--",color="r")

plt.xlabel("n cells")
plt.ylabel("L2 norm - simulation vs analytical")

plt.savefig("convergence.png")
