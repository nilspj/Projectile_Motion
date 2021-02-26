# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import common as co

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

def ConvTest(): #Convergence test 
    N = 1000#number of different time steps
    P = 101
    stepList = np.zeros(N)
    diffList = np.zeros(N)
    for i in range(10,N+10):
        h = P/i
        stepList[i-10]=i
        diffList[i-10]= np.abs(co.rungekutta(co.F,h,w,i)[2][i-1]-co.analytical()[2][999])
        
        if diffList[i-10] <= 1:
            print("Enough points:", i-10,", difference:",diffList[i-10],",Time step:",h)
            
    print("Final number of points:", N-1,", difference:",diffList[N-1],",Time step:",P/(N-1))
    return stepList,diffList 

Test1=ConvTest()
plt.plot(Test1[0],Test1[1],linewidth = 1.5)
plt.ylabel(r'$\Delta r\,[\mathrm{m}]$', fontsize=14)
plt.xlabel(r'$N\,[\mathrm{1}]$', fontsize=14)
plt.savefig("ConvTest.pdf", bbox_inches='tight')