# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import common as co

P = 250
N =int(P/0.1)
h =P/N

def OptimalAngle(f):
    rangeList = np.zeros(91)
    #heightList = np.zeros(91)
    #timeList = np.zeros(91)
    degList = np.linspace(0,90,91)
    for i in range(0,91):
        y0 = 0 
        x0 = 0 
        v0x =700*np.cos(co.DtR(i)) 
        v0y = 700*np.sin(co.DtR(i))
        w = np.array([x0,y0,v0x,v0y]) 
        rangeList[i] = co.rungekutta(f,h,w,N)[5] # 5 xlanding, 6 ymax, 7 tim
        #heightList[i] = co.rungekutta(f,h,w,N)[6]
        #timeList[i] = co.rungekutta(f,h,w,N)[7]
    print("Max Range:",rangeList[np.argmax(rangeList)],",Angle:",np.argmax(rangeList))
    #print("Max Range:",heightList[np.argmax(rangeList)],",Angle:",np.argmax(rangeList))
    #print("Max Range:",timeList[np.argmax(rangeList)],",Angle:",np.argmax(rangeList))
    return rangeList,degList

Angle1 = OptimalAngle(co.F)
Angle2 =OptimalAngle(co.D)
Angle3 =OptimalAngle(co.I)
Angle4 =OptimalAngle(co.A)

plt.plot(Angle1[1],Angle1[0],linewidth = 0.7, label='No Drag')
plt.plot(Angle2[1],Angle2[0],linewidth = 0.7, label='Drag')
plt.plot(Angle3[1],Angle3[0],linewidth = 0.7, label='Isothermal')
plt.plot(Angle4[1],Angle4[0],linewidth = 0.7, label='Adiabatic')
plt.ylabel(r'$x_{L}\,[\mathrm{m}]$', fontsize=14)
plt.xlabel(r'$\theta\,[\mathrm{^{\circ}}]$', fontsize=14)
plt.legend(loc=2)
plt.savefig("OptimalAngle.pdf", bbox_inches='tight')