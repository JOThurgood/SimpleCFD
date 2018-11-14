import numpy as np
import matplotlib.pyplot as plt
import os 

xc = np.fromfile('xc.dat', dtype = np.float64)
yc  = np.fromfile('yc.dat', dtype = np.float64)
u = np.fromfile('u.dat', dtype = np.float64)
v = np.fromfile('v.dat', dtype = np.float64)
divu = np.fromfile('divu.dat', dtype = np.float64)
curlu = np.fromfile('curlu.dat', dtype = np.float64)
#phi = np.fromfile('phi.dat', dtype = np.float64)

u = u.reshape(yc.size, xc.size)
v = v.reshape(yc.size, xc.size)
divu = divu.reshape(yc.size, xc.size)
curlu = curlu.reshape(yc.size, xc.size)
#phi = phi.reshape(yc.size, xc.size)


if not os.path.exists("output"):
        os.makedirs("output")


filename = "output/stream_"
filenum = 0 
while (os.path.exists(filename+str(filenum)+".png")):
    filenum+=1

filename = "{}{}.png".format(filename,filenum)

plt.streamplot(xc,yc,u,v)
plt.contourf(xc,yc,abs(u**2+v**2))
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y') 
plt.title('abs(U) + stream at dump #{}'.format(filenum))
plt.axis([xc[0],xc[xc.size-1],yc[0],yc[yc.size-1]])
plt.savefig(filename)



filename = "output/u_"
filenum = 0 
while (os.path.exists(filename+str(filenum)+".png")):
    filenum+=1

filename = "{}{}.png".format(filename,filenum)


plt.clf()
plt.contourf(xc, yc, u, 255, cmap=plt.cm.RdBu_r)
#plt.grid(True)
plt.xlabel('x')
plt.ylabel('y') 
plt.title('u')
plt.colorbar()
#plt.show()
plt.savefig(filename)


filename = "output/v_"
filenum = 0 
while (os.path.exists(filename+str(filenum)+".png")):
    filenum+=1

filename = "{}{}.png".format(filename,filenum)


plt.clf()
plt.contourf(xc, yc, v, 255, cmap=plt.cm.RdBu_r)
#plt.grid(True)
plt.xlabel('x')
plt.ylabel('y') 
plt.title('v')
plt.colorbar()
#plt.show()
plt.savefig(filename)

filename = "output/divu_"
filenum = 0 
while (os.path.exists(filename+str(filenum)+".png")):
    filenum+=1

filename = "{}{}.png".format(filename,filenum)
plt.clf()
plt.contourf(xc, yc, divu, 255, cmap=plt.cm.RdBu_r)
#plt.grid(Trve)
plt.xlabel('x')
plt.ylabel('y') 
plt.title('divu')
plt.colorbar()
#plt.show()
plt.savefig(filename)

filename = "output/curlu_"
filenum = 0 
while (os.path.exists(filename+str(filenum)+".png")):
    filenum+=1
filename = "{}{}.png".format(filename,filenum)

plt.clf()
plt.contourf(xc, yc, curlu, 255, cmap=plt.cm.RdBu_r)
#plt.grid(Trve)
plt.xlabel('x')
plt.ylabel('y') 
plt.title('curl u')
plt.colorbar()
#plt.show()
plt.savefig(filename)

#plt.clf()
#plt.contourf(xc, yc, phi, 255) 
##plt.grid(Trve)
#plt.xlabel('x')
#plt.ylabel('y') 
#plt.title('phi')
#plt.colorbar()
##plt.show()
#plt.savefig('phi.png')



