#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 20:26:05 2018

@author: rivertz
"""

import numpy as np
import math as m
import sys
import time
import matplotlib.pyplot as plot

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
def F(Y):
    G=6.67408e-11
    m1=7.3477e22 #Mass Moon
    m2=5.9736e24 #Mass Earth
    Gm2=G*m2;

    px2=0;py2=0; #Position of the earth. Needs to be changed in two-planet problem!
    px1=Y[1]; #x0
    py1=Y[3]; #y0
    vx1=Y[2]; #vx0
    vy1=Y[4]; #vy0
    dist=np.sqrt((px2-px1)**2+(py2-py1)**2);
    z=np.zeros(Y.size);

    #Equation from the book for one-body-problem
    #Gives the change in velocity
    #t0,x0,vx0,y0,vx0
    z[0]=1 #Vi vil ikke endre, fordi det skal ganges med h for å få neste t.
    z[1]=vx1 #x' = vx
    z[2]=(Gm2*(px2-px1))/(dist**3) # (vx)' = (Change in velocity in x-direction per time)
    z[3]=vy1 #y' = vy
    z[4]=(Gm2*(py2-py1))/(dist**3) #(vy)' = (Change in velocity in y-direction per time)
    return z #Returns s1 or s2

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