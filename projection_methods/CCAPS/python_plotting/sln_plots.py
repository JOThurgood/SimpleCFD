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
divu = np.fromfile('divu.dat', dtype = np.float64)
curlu = np.fromfile('curlu.dat', dtype = np.float64)
#phi = np.fromfile('phi.dat', dtype = np.float64)

# Reshape arrays
u = u.reshape(yc.size, xc.size)
v = v.reshape(yc.size, xc.size)
divu = divu.reshape(yc.size, xc.size)
curlu = curlu.reshape(yc.size, xc.size)
#phi = phi.reshape(yc.size, xc.size)

# Check the output directory exists and mkdir if nescessary
if not os.path.exists("output"):
        os.makedirs("output")

# Stream plot with abs(vel)

plt.streamplot(xc,yc,u,v)
plt.contourf(xc,yc,abs(u**2+v**2),255, cmap=plt.cm.Greys)
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y') 
plt.title('abs vel + streamlines at t = {:.3f}'.format(float(time)))
plt.axis([xc[0],xc[xc.size-1],yc[0],yc[yc.size-1]])
plt.savefig(getfilenameindex("output/stream_"))


# u data

plt.clf()
plt.contourf(xc, yc, u, 255, cmap=plt.cm.RdBu_r)
plt.xlabel('x')
plt.ylabel('y') 
plt.title('u at t = {:.3f}'.format(float(time)))
plt.colorbar()
plt.savefig(getfilenameindex("output/u_"))

# v data

plt.clf()
plt.contourf(xc, yc, v, 255, cmap=plt.cm.RdBu_r)
plt.xlabel('x')
plt.ylabel('y') 
plt.title('v at t = {:.3f}'.format(float(time)))
plt.colorbar()
plt.savefig(getfilenameindex("output/v_"))


plt.clf()
plt.contourf(xc, yc, divu, 255, cmap=plt.cm.RdBu_r)
plt.xlabel('x')
plt.ylabel('y') 
plt.title('div(U) at t = {:.3f}'.format(float(time)))
plt.colorbar()
plt.savefig(getfilenameindex("output/div_"))


plt.clf()
plt.contourf(xc, yc, curlu, 255, cmap=plt.cm.RdBu_r)
plt.xlabel('x')
plt.ylabel('y') 
plt.title('curl(U) at t = {:.3f}'.format(float(time)))
plt.colorbar()
plt.savefig(getfilenameindex("output/curl_"))

#plt.clf()
#plt.contourf(xc, yc, phi, 255) 
##plt.grid(Trve)
#plt.xlabel('x')
#plt.ylabel('y') 
#plt.title('phi')
#plt.colorbar()
##plt.show()
#plt.savefig('phi.png')



