# -*- coding: utf-8 -*-
import numpy as np

def DtR(n):  #converts from degree to rad
    return n*np.pi/180

T = 101# time period #101 normal #205 big bertha
N = 350#250#h*N=1 number of points
h = T/N #time step

y0 = 0 #start position in y-direction
x0 = 0  #start postition in x-direction
t0 = 0 #yr
v0= 700
v0x = v0*np.cos(DtR(45)) #m/s^2   start velocity in x-direction
v0y = v0*np.sin(DtR(45)) #m/ss

B_m=4*10**(-5)*(50/106)
g = 9.81 #m/s^2
Y0= 10**(4)
a=6.5*10**(-3)
alpha = 2.5
T0= 288 #kelvin

w = np.array([x0,y0,v0x,v0y])

def DtR(n):  #converts from degree to rad
    return n*np.pi/180

def F(w):#without drag
    xdot = w[2]
    ydot = w[3]
    ax = 0
    ay = -1*g 
    return np.array([xdot,ydot,ax,ay])

def D(w):#with drag
    xdot = w[2]
    ydot = w[3]
    ax = -B_m*(w[2]**2)
    ay = -1*g -B_m*(w[3]**2) 
    return np.array([xdot,ydot,ax,ay])


def I(w):#Isothermal
    xdot = w[2]
    ydot = w[3]
    ax = -B_m*(w[2]**2)*np.exp(-w[1]/Y0)
    ay = -1*g -B_m*(w[3]**2)*np.exp(-w[1]/Y0)
    return np.array([xdot,ydot,ax,ay])

def A(w):#Adiabatic
    xdot = w[2]
    ydot = w[3]
    if(w[1]<(T0/a)):
        ax = -B_m*(w[2]**2)*(np.abs(1-(a*w[1])/T0))**alpha
        ay = -1*g -B_m*(w[3]**2)*(np.abs(1-(a*w[1])/T0))**alpha
    else:
        ax= 0
        ay= -1*g
    return np.array([xdot,ydot,ax,ay])

def RK4(f,h,w):
    
    s1 = f(w)
    s2 = f(w + h*s1/2)
    s3 = f(w + h*s2/2)
    s4 = f(w + h*s3)
    return (w + h*(s1 + 2*s2 + 2*s3 + s4)/6)

def rungekutta(f,h,w,N):
    RKxList = np.zeros(N)
    RKyList = np.zeros(N)
    RKrList = np.zeros(N)
    RKvxList = np.zeros(N)
    RKvyList = np.zeros(N)
    for i in range(N):
        #t = t0 + h
        #t0 = t
        #r = dist(x0,y0)
        RKxList[i]=w[0]
        RKyList[i]=w[1]
        RKrList[i]=np.sqrt((w[0])**2 + (w[1])**2)
        RKvxList[i] = w[2]
        RKvyList[i] = w[3]
        w1 = RK4(f,h,w)
        w=w1
      
    x_landing=0
    flight_time=0
    for i in range(1,N):
        if RKyList[i]<0:
            r=-RKyList[i-1]/RKyList[i]
            x_landing=(RKxList[i-1]+r*RKxList[i])/(r+1)
            flight_time=i*h
            break
    #print(x_landing)   
    
    y_highest=np.max(RKyList)
    #print(y_highest)
    return RKxList,RKyList,RKrList,RKvxList,RKvyList,x_landing,y_highest,flight_time

def analytical():
    N=1000#number of different time steps
    P = 101
    h=P/N
    vy = np.zeros(N)
    vx=np.zeros(N)
    sx=np.zeros(N)
    sy=np.zeros(N)
    sr=np.zeros(N)
    
    for i in range(0,N):
        
        vx[i]=v0x
        vy[i]=-g*h*i+v0y
        sx[i]=v0x*h*i
        sy[i]=-0.5*g*((h*i)**2)+v0y*h*i
        sr[i]=np.sqrt(sx[i]**2+sy[i]**2)
    
    #print(v0x*(2*v0y/g))  #landing point
    
    return sx,sy,sr


def OptimalAngle(f):
    rangeList = np.zeros(91)
    degList = np.linspace(0,90,91)
    for i in range(0,91):
        y0 = 0 
        x0 = 0 
        v0x =700*np.cos(DtR(i)) 
        v0y = 700*np.sin(DtR(i))
        w = np.array([x0,y0,v0x,v0y]) 
        rangeList[i] = rungekutta(f,h,w,N)[5] # 5 xlanding, 6 ymax, 7 time
    return rangeList,degList

def BigBertha():
    v0= 1640
    y0 = 0 
    x0 = 0 
    v0x =v0 *np.cos(DtR(45)) 
    v0y = v0* np.sin(DtR(45))
    w = np.array([x0,y0,v0x,v0y])
    RKbb = rungekutta(A,h,w,N)
    
    x_landing=RKbb[5]
    print(x_landing)
    y_highest=RKbb[6]
    
    time=RKbb[7]
    
    return RKbb[0],RKbb[1]