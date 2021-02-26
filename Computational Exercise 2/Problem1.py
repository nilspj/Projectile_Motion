# -*- coding: utf-8 -*-
"""
Created on Sat Oct 20 16:24:07 2018

@author: Nils og Therese
"""

import common as co
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

P =320
N =int(300/0.1)
h =P/N

phi0 = 49.605
theta0 = 3.514722
x0 = 6371*10**3*np.sin(co.DtR(phi0))*np.cos(co.DtR(theta0))
y0 = 6371*10**3*np.sin(co.DtR(phi0))*np.sin(co.DtR(theta0))
z0 = 6371*10**3*np.cos(co.DtR(phi0))
r0 = np.sqrt(x0**2+y0**2+z0**2)

v0= 1640
angle = 20
v0x = v0*np.sin(co.DtR(phi0+angle))*np.cos(co.DtR(theta0)) #m/s^2   start velocity in x-direction
v0y = v0*np.sin(co.DtR(phi0+angle)*np.sin(co.DtR(theta0))) #m/ss
v0z = v0*np.cos(co.DtR(phi0+angle))
vrot = 7.29*10**(-5)


w = np.array([x0,y0,z0,v0x,v0y,v0z])


fig = plt.figure()
ax = fig.gca(projection='3d')

Projectile1 = co.rungekutta(co.A,h,w,N)
print(Projectile1[7])

ax.plot(Projectile1[0], Projectile1[1], Projectile1[2], label='parametric curve')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.legend()
plt.show()

#sphere
v, u = np.mgrid[co.DtR(49):co.DtR(54):20j, co.DtR(1):co.DtR(5):10j]
x = r0*np.cos(u)*np.sin(v)
y = r0*np.sin(u)*np.sin(v)
z = r0*np.cos(v)
ax.plot_wireframe(x, y, z, color="r")
plt.show()
