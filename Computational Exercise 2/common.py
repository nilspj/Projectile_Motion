# -*- coding: utf-8 -*-
import numpy as np

def DtR(n):  #converts from degree to rad
    return n*np.pi/180

T = 101# time period #101 normal #205 big bertha
N = 350#250#h*N=1 number of points
h = T/N #time step

phi0 = 49.605
theta0 = 3.514722
x0 = 6371*10**3*np.sin(DtR(phi0))*np.cos(DtR(theta0))
y0 = 6371*10**3*np.sin(DtR(phi0))*np.sin(DtR(theta0))
z0 = 6371*10**3*np.cos(DtR(phi0))
r0 = np.sqrt(x0**2+y0**2+z0**2)


t0 = 0 #yr
v0= 1640
v0x = v0*np.sin(DtR(45))*np.cos(DtR(0)) #m/s^2   start velocity in x-direction
v0y = v0*np.sin(DtR(45))*np.sin(DtR(0)) #m/ss
v0z = v0*np.cos(DtR(45))
vrot = 7.29*10**(-5)

B_m=4*10**(-5)*(50/106)
g = 9.81 #m/s^2
Y0= 10**(4)
a=6.5*10**(-3)
alpha = 2.5
T0= 288 #kelvin

w = np.array([x0,y0,z0,v0x,v0y,v0z])

def DtR(n):  #converts from degree to rad
    return n*np.pi/180

def F(w):#without drag
    xdot = w[3]
    ydot = w[4]
    zdot = w[5]
    ax = 0
    ay = 0
    az = -1*g
    return np.array([xdot,ydot,zdot,ax,ay,az])

def D(w):#with drag
    xdot = w[3]
    ydot = w[4]
    zdot = w[5]
    ax = -B_m*(w[2]**2)
    ay = -B_m*(w[4]**2)
    az = -1*g -B_m*(w[5]**2)
    return np.array([xdot,ydot,zdot,ax,ay,az])


def I(w):#Isothermal
    xdot = w[3]
    ydot = w[4]
    zdot = w[5]
    ax = -B_m*(w[3]**2)*np.exp(-w[2]/Y0)
    ay = -B_m*(w[4]**2)*np.exp(-w[2]/Y0)
    az = -1*g -B_m*(w[5]**2)*np.exp(-w[2]/Y0)
    return np.array([xdot,ydot,zdot,ax,ay,az])

def A(w):#Adiabatic
    xdot = w[3]
    ydot = w[4]
    zdot = w[5]
    if(w[1]<(T0/a)):
        ax =-1*g*w[0]/r0-B_m*(w[3]**2)*(np.abs(1-(a*w[2])/T0))**alpha
        ay = -1*g*w[1]/r0 -B_m*(w[4]**2)*(np.abs(1-(a*w[2])/T0))**alpha
        az = -1*g*w[2]/r0 -B_m*(w[5]**2)*(np.abs(1-(a*w[2])/T0))**alpha
    else:
        ax= -1*g*w[0]/r0
        ay= -1*g*w[1]/r0
        az = -1*g*w[2]/r0
    return np.array([xdot,ydot,zdot,ax,ay,az])

def RK4(f,h,w):
    
    s1 = f(w)
    s2 = f(w + h*s1/2)
    s3 = f(w + h*s2/2)
    s4 = f(w + h*s3)
    return (w + h*(s1 + 2*s2 + 2*s3 + s4)/6)

def rungekutta(f,h,w,N):
    RKxList = np.zeros(N)
    RKyList = np.zeros(N)
    RKzList = np.zeros(N)
    RKrList = np.zeros(N)
    RKvxList = np.zeros(N)
    RKvyList = np.zeros(N)
    RKvzList = np.zeros(N)
    for i in range(N):
        #t = t0 + h
        #t0 = t
        #r = dist(x0,y0)
        RKxList[i]=w[0]
        RKyList[i]=w[1]
        RKzList[i]=w[2]
        RKrList[i]=np.sqrt(w[0]**2 + w[1]**2 +w[2]**2)
        RKvxList[i] = w[3]
        RKvyList[i] = w[4]
        RKvzList[i] = w[5]
        w1 = RK4(f,h,w)
        w=w1
      
    d_landing=0
    x_landing=0
    y_landing=0
    z_landing=0
    phi_landing=0
    flight_time=0
    for i in range(1,N):
        if RKrList[i]<r0:
            r= RKrList[i-1]/RKrList[i]
            x_landing=(RKxList[i-1]+r*RKxList[i])/(r+1)
            y_landing=(RKyList[i-1]+r*RKyList[i])/(r+1)
            z_landing=(RKzList[i-1]+r*RKzList[i])/(r+1)
            phi_landing = np.arccos(z_landing/r0)-DtR(phi0) 
            print(phi_landing)
            d_landing = r0*phi_landing
            flight_time=i*h
            break
    
    
    r_highest=np.max(RKrList)
    #print(y_highest)
    return RKxList,RKyList,RKzList,RKrList,RKvxList,RKvyList,RKvzList,d_landing,r_highest,flight_time

def analytical():
    N=1000#number of different time steps
    P = 101
    h=P/N
    vy = np.zeros(N)
    vx=np.zeros(N)
    vz=np.zeros(N)
    sx=np.zeros(N)
    sy=np.zeros(N)
    sz=np.zeros(N)
    sr=np.zeros(N)
    
    for i in range(0,N):
        
        vx[i]=v0x
        vy[i]=v0y
        vz[i]=-g*h*i+v0z
        sx[i]=v0x*h*i
        sy[i]=v0y*h*i
        sz[i]=-0.5*g*((h*i)**2)+v0z*h*i
        sr[i]=np.sqrt(sx[i]**2+sy[i]**2+sz[i]**2)
    
    #print(v0x*(2*v0y/g))  #landing point
    
    return sx,sy,sz,sr


def OptimalAngle(f):
    rangeList = np.zeros(91)
    degList = np.linspace(0,90,91)
    for i in range(0,91):
        y0 = 0 
        x0 = 0
        z0 = 0
        v0x =700*np.sin(DtR(i))*np.cos(DtR(0)) 
        v0y = 700*np.sin(DtR(i))*np.sin(DtR(0)) 
        v0z = 700*np.cos(DtR(i))
        w = np.array([x0,y0,z0,v0x,v0y,v0z]) 
        rangeList[i] = rungekutta(f,h,w,N)[5] # 5 xlanding, 6 ymax, 7 time
    return rangeList,degList

def BigBertha():
    v0 = 1640
    z0 = 0
    y0 = 0 
    x0 = 0 
    v0x =v0 *np.sin(DtR(45))*np.cos(DtR(0)) 
    v0y = v0* np.sin(DtR(45))*np.sin(DtR(0))
    v0z = v0* np.cos(DtR(45))
    w = np.array([x0,y0,z0,v0x,v0y,v0z])
    RKbb = rungekutta(A,h,w,N)
    
    z_landing=RKbb[7]
    print(z_landing)
    z_highest=RKbb[8]
    
    time=RKbb[9]
    
    return RKbb[0],RKbb[1],RKbb[2]