# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import common as co

P =300
N =int(300/0.1)
h =P/N

y0 = 0 #start position in y-direction
x0 = 0  #start postition in x-direction
t0 = 0 #yr
v0= 1640
v0x = v0*np.cos(co.DtR(45)) #m/s^2   start velocity in x-direction
v0y = v0*np.sin(co.DtR(45)) #m/ss

w = np.array([x0,y0,v0x,v0y])

B_m=4*10**(-5)*(50/106)


def BigBerthaAngle():
    rangeList = np.zeros(91)
    heightList = np.zeros(91)
    timeList = np.zeros(91)
    degList = np.linspace(0,90,91)
    for i in range(0,91):
        y0 = 0 
        x0 = 0 
        v0x =v0*np.cos(co.DtR(i)) 
        v0y = v0*np.sin(co.DtR(i))
        w = np.array([x0,y0,v0x,v0y]) 
        rangeList[i] = co.rungekutta(co.A,h,w,N)[5] # 5 xlanding, 6 ymax, 7 tim
        heightList[i] = co.rungekutta(co.A,h,w,N)[6]
        timeList[i] = co.rungekutta(co.A,h,w,N)[7]
    print("Max Range:",rangeList[np.argmax(rangeList)],",Angle:",np.argmax(rangeList))
    print("Max height:",heightList[np.argmax(heightList)],",Angle:",np.argmax(heightList))
    print("total time:",timeList[np.argmax(timeList)],",Angle:",np.argmax(timeList))
    return degList,rangeList,heightList,timeList

BigB = co.rungekutta(co.A,h,w,N)
BB = BigBerthaAngle()
'''
#projectile path
plt.plot(BigB[0],BigB[1],linewidth = 0.7, label='BigBertha')
plt.ylabel(r'$y\,[\mathrm{m}]$', fontsize=14)
plt.xlabel(r'$x\,[\mathrm{m}]$', fontsize=14)
plt.ylim(ymin=0) #Stops plot from plotting negative y-values
#plt.legend(loc=2)
plt.savefig("BB.pdf", bbox_inches='tight')

#Maximum landing distance
plt.plot(BB[0],BB[1],linewidth = 0.7, label='BigBerthaX')
plt.ylabel(r'$x_{L}\,[\mathrm{m}]$', fontsize=14)
plt.xlabel(r'$\theta\,[\mathrm{^{\circ}}]$', fontsize=14)
#plt.legend(loc=2)
plt.savefig("BBx.pdf", bbox_inches='tight')

#Maximum heigh 
plt.plot(BB[0],BB[2],linewidth = 0.7, label='BigBerthaY')
plt.ylabel(r'$y_{max}\,[\mathrm{m}]$', fontsize=14)
plt.xlabel(r'$\theta\,[\mathrm{^{\circ}}]$', fontsize=14)
#plt.legend(loc=2)
plt.savefig("BBy.pdf", bbox_inches='tight')
'''
#Total time 
plt.plot(BB[0],BB[3],linewidth = 0.7, label='BigBerthaT')
plt.ylabel(r'$t_{tot}\,[\mathrm{m}]$', fontsize=14)
plt.xlabel(r'$\theta\,[\mathrm{^{\circ}}]$', fontsize=14)
#plt.legend(loc=2)
plt.savefig("BBt.pdf", bbox_inches='tight')
