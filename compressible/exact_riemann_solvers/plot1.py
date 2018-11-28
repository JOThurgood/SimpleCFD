import numpy as np
import matplotlib.pyplot as plt

x = np.fromfile('x.dat', dtype=np.float64)
rho = np.fromfile('rho.dat', dtype=np.float64)
p = np.fromfile('p.dat', dtype=np.float64)
u = np.fromfile('u.dat', dtype=np.float64)
en = np.fromfile('en.dat', dtype=np.float64)

#title = 'Test 1 - Sod at t=0.25'
#title = 'Test 2 - 1-2-3 problem at t=0.15'
#title = 'Test 3 - left half W+C blast at t=0.012'
#title = 'Test 4 - right half W+C blast at t=0.035'
#title = 'Test 5 - full W+C blast at t=0.035'
title = 'your title here!'
plt.plot(x, rho)
if x.size < 129 :
  plt.plot(x,rho,'bx')
plt.xlabel('x')
plt.ylabel('rho')
plt.title(title)
plt.grid(True)
#plt.show()
plt.savefig("rho.png")

plt.clf()
plt.plot(x, p)
if x.size < 129 :
  plt.plot(x,p,'bx')
plt.xlabel('x')
plt.ylabel('p')
plt.title(title)
plt.grid(True)
#plt.show()
plt.savefig("p.png")

plt.clf()
plt.plot(x, u)
if x.size < 129 :
  plt.plot(x,u,'bx')
plt.xlabel('x')
plt.ylabel('u')
plt.title(title)
plt.grid(True)
#plt.show()
plt.savefig("u.png")

plt.clf()
plt.plot(x, en)
if x.size < 129 :
  plt.plot(x,en,'bx')
plt.xlabel('x')
plt.ylabel('en')
plt.title(title)
plt.grid(True)
#plt.show()
plt.savefig("en.png")


