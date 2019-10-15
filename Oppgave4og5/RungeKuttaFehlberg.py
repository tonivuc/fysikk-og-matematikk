#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 20:26:05 2018

@author: rivertz
"""
from SaturnV import SaturnV
import numpy as np
import math as m
import sys
import time
import matplotlib.pyplot as plot

#Mass of different stages (kilograms)
m1 = 2970000;
m2 = 680000;
m3 = 183000;
#Burn time of different stages (seconds)
b1 = 168;
b2 = 360;
b3 = 165;
#Fuel consumption of different stages (kilograms per second)
c1 = 12857.1;
c2 = 1266.9;
c3 = 219;
#Thrust of different stages (Newtons)
t1 = 35100000;
t2 = 5141000;
t3 = 1033100;

class RungeKuttaFehlberg54:
    localError = 0;
    A=np.array(
        [[   0     ,    0     ,    0      ,    0     ,  0   ,0],
         [   1/4   ,    0     ,    0      ,    0     ,  0   ,0],
         [   3/32  ,    9/32  ,    0      ,    0     ,  0   ,0],
         [1932/2197,-7200/2197, 7296/2197 ,    0     ,  0   ,0],
         [ 439/216 ,   -8     , 3680/513  , -845/4104,  0   ,0],
         [  -8/27  ,    2     ,-3544/2565 , 1859/4104,-11/40,0]])

    B=np.array(
        [[  25/216 ,    0     , 1408/2565 , 2197/4104 ,-1/5 ,0],
         [  16/135 ,    0     , 6656/12825,28561/56430,-9/50,2/55]]);

    def __init__(self,
                 function,
                 dimension,
                 stepsize,
                 tolerance):
        self.F=function;
        self.dim=dimension;
        self.h=stepsize;
        self.tol=tolerance;

    def step(self,
             Win):
        s=np.zeros((6,self.dim))

        for i in range(0,6):
            s[i,:]=self.F(Win+self.h*self.A[i,0:i].dot(s[0:i,:]))

        Zout=Win+self.h*(self.B[0,:].dot(s));
        Wout=Win+self.h*(self.B[1,:].dot(s));

        E=np.linalg.norm(Wout-Zout,2)/np.linalg.norm(Wout,2);
        return Wout, E

    def safeStep(self,
                 Win):
        Wout,E = self.step(Win);
        # Check if the error is tolerable
        if(not self.isErrorTolerated(E)):
            #Try to adjust the optimal step length
            self.adjustStep(E);
            Wout,E = self.step(Win);
        # If the error is still not tolerable
        counter=0;
        while(not self.isErrorTolerated(E)):
            #Try if dividing the steplength with 2 helps.
            self.divideStepByTwo();
            Wout,E = self.step(Win);
            counter = counter + 1;
            if(counter>10):
                sys.exit(-1);

        self.adjustStep(E);
        self.localError += E;
        return Wout, E

    def isErrorTolerated(self,E):
        return E<self.tol;

    def adjustStep(self,E):
        if(E==0):
            s=2;
        else:
            s=m.pow(self.tol*self.h/(2*E),0.25);
        self.h = s * self.h;

    def divideStepByTwo(self):
        self.h=self.h/2;

    def setStepLength(self,stepLength):
        self.h=stepLength;

#The actual equations to be "solved" numerically and their initial values are to be set here
#F = ydot = funksjon for y derviert og y dobbelderivert
"""
def F(Y):
    G = self.GravConst
    Gm = G * self.m
    t = x[0]
    px2 = 0
    py2 = 0
    px1 = x[1]
    py1 = x[3]
    vx1 = x[2]
    vy1 = x[4]
    dist = np.sqrt((px2 - px1) ** 2 + (py2 - py1) ** 2)
    z = np.zeros(5)

    # Force from gravity on rocket divided by rocket mass
    Fg_x = (Gm * (px2 - px1)) / (dist ** 3)
    Fg_y = (Gm * (py2 - py1)) / (dist ** 3)

    # Force from air drag on rocket divided by rocket mass
    absolute_velocity = np.sqrt(vx1*vx1 + vy1*vy1)
    saturnV = SaturnV(m1,c1,b1,t1,m2,c2,b2,t2,m3,c3,b3,t3)
    Fd = self.get_air_drag(self.moh(), 1/2, 50, absolute_velocity) / saturnV.calculateMass(t)

    # Force from thrusters on rocket divided by rocket mass
    F = saturnV.calculateThrust(t)

    self.acceleration = F/saturnV.calculateMass(t) - (Fg_y + Fd)

    z[0] = 1
    z[1] = vx1
    z[2] = Fg_x
    z[3] = vy1
    z[4] = F/saturnV.calculateMass(t) - (Fg_y + Fd)

    self.xy[0].append(self.get_position()[0])
    self.xy[1].append(self.get_position()[1])

    return z
"""
def main():
    # W = initial values
    #          [t,x,vx,y,vy]
    W  =np.array([0.0,0.0, 1082, 362570e3, 0.0]);
    h=0.25; #Step size (0.1 t)
    tol=05e-14; #RelativeS error, or just error?
    tEnd=10.0; #Value for t where we stop the approximation
    rkf54 = RungeKuttaFehlberg54(F,W.size,h,tol) #function, dimension, stepsize, tolerance. Why is dimension = 4?

    print('kjører main?')
    while(W[0]<tEnd): #Why
        W , E = rkf54.safeStep(W); #Returns a new vector of approximated y-values W. And an error vector?

    print("Now we are fucking the step length")
    rkf54.setStepLength(tEnd-W[0]); #Makes no sense. This would make the step length 1.9 after the first iteration.
    W,E = rkf54.step(W); #Make one step forward.

    print(W,E);

if __name__ == "__main__":
    # execute only if run as a script
    main()
