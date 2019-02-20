import numpy as np
import matplotlib.pyplot as plt
import os 


def getfilenameindex(filename):
  filenum = 0
  while (os.path.exists(filename+"{:04}".format(filenum)+".png")):
    filenum+=1
  filename = "{}{:04}.png".format(filename,filenum)
  return filename

# Load the data

time = np.fromfile('time.dat', dtype = np.float64)
xc = np.fromfile('xc.dat', dtype = np.float64)
yc  = np.fromfile('yc.dat', dtype = np.float64)
u = np.fromfile('u.dat', dtype = np.float64)
v = np.fromfile('v.dat', dtype = np.float64)

# Reshape arrays
u = u.reshape(yc.size, xc.size)
v = v.reshape(yc.size, xc.size)

cen = yc.size // 2
plt.plot(yc,(u[:,cen]+u[:,cen+1])/2.)


y_ghia = np.array([
  0.0000, 
  0.0547,
  0.0625,
  0.0703,
  0.1016,
  0.1719, 
  0.2813,
  0.4531,
  0.5,
  0.6172,
  0.7344,
  0.8516,
  0.9531,
  0.9609,
  0.9688,
  0.9766,
  1.])

u_center_re1000 = np.array([
1.0, 
0.65928,
0.57492,
0.51117,
0.46604,
0.33304,
0.18719,
0.05702,
-0.06080,
-0.10648,
-0.27805,
-0.38289,
-0.29730,
-0.22220,
-0.20196,
-0.18109,
0.00000 ])

u_center_re1000 = np.flip(u_center_re1000,0)
plt.scatter(y_ghia,u_center_re1000)
plt.scatter(yc,(u[:,cen]+u[:,cen+1])/2.,color='red')
plt.title('u on central vertical axis at t = {:.3f}'.format(float(time)))
plt.savefig(getfilenameindex("output/ghia_u_cent_"))
