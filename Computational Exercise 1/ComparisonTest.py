# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import common as co

N = 971
P = 101
h = P/N

y0 = 0 #start position in y-direction
x0 = 0  #start postition in x-direction
t0 = 0 #yr
v0= 700
v0x = v0*np.cos(co.DtR(45)) #m/s^2   start velocity in x-direction
v0y = v0*np.sin(co.DtR(45)) #m/ss

B_m=4*10**(-5)
g = 9.81 #m/s^2
Y0= 10**(4)
a=6.5*10**(-3)
alpha = 2.5
T0= 288 #kelvin

w = np.array([x0,y0,v0x,v0y])

pNum1 = co.rungekutta(co.F,h,w,N)
pAna = co.analytical()
plt.plot(pNum1[0],pNum1[1],linewidth = 0.7, label='Numerical')
plt.plot(pAna[0],pAna[1],linewidth = 0.7, label='Analytical')
plt.ylabel(r'$y\,[\mathrm{m}]$', fontsize=14)
plt.xlabel(r'$x\,[\mathrm{m}]$', fontsize=14)
plt.ylim(ymin=0) #Stops plot from plotting negative y-values
plt.legend(loc=2)
plt.savefig("NumAna.pdf", bbox_inches='tight')