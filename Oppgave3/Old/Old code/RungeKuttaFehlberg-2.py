#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 20:26:05 2018

@author: rivertz
"""

import numpy as np
import math as m
import sys

class RungeKuttaFehlberg54:

    #Constants used together with w_i
    A=np.array(
        [[   0     ,    0     ,    0      ,    0     ,  0   ,0],
         [   1/4   ,    0     ,    0      ,    0     ,  0   ,0],
         [   3/32  ,    9/32  ,    0      ,    0     ,  0   ,0],
         [1932/2197,-7200/2197, 7296/2197 ,    0     ,  0   ,0],
         [ 439/216 ,   -8     , 3680/513  , -845/4104,  0   ,0],
         [  -8/27  ,    2     ,-3544/2565 , 1859/4104,-11/40,0]])

    #Constants for Zout and Wout
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

    #This function is used all the time
    def step(self,
             Win):
        s=np.zeros((6,self.dim))
        v=np.zeros((6,3))

        #Calculate all 6 s-values/vectors
        for i in range(0,6):
            s[i,:]=self.F(Win+self.h*self.A[i,0:i].dot(s[0:i,:]),v)
            print(s[i,:])

        Zout=Win+self.h*(self.B[0,:].dot(s)); #In the book, Zout is better than Wout, given that Zout is the locally extrapolated version.
        Wout=Win+self.h*(self.B[1,:].dot(s));

        E=np.linalg.norm(Wout-Zout,2)/np.linalg.norm(Wout,2);
        return Wout, E

    def safeStep(self,
                 Win): #Win are the current approximations for y
        Wout,E = self.step(Win); #Makes a normal step
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
                sys.exit(-1); #Some sort of way to get out of an inifinite step size reducing loop

        self.adjustStep(E);

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
#F = ydot = y derivert
def F(Y, velocitiesIn):

    #Y er en vektor med 3 verdier w; x y og z


    """
    xDot = v_x \\
    v_x' &= -G\frac{m_2 x}{(x^2+y^2)^{3/2}} \\
    y' &= v_y \\
    v_y' &= -G\frac{m_2 y}{(x^2+y^2)^{3/2}}
    """
    M=np.array([[0.49119653, 0.32513304, 0.98057799],
                [0.20768544, 0.97699416, 0.18220559],
                [0.96407071, 0.18373237, 0.95307793]]);
    res=np.ones(4);
    res[1:4]=M.dot(Y[1:4]); #Mutiply each x, y, z with
    return res;


def main():
    W  =np.array([0,1,1,1]); #Initial values for [t,x,y,z]
    velocities = np.array([0,-1,0]) #Initial velocities in x, y and z direction
    h=0.1; #Step size
    tol=05e-14; #RelativeS error, or just error?
    tEnd=2.0; #Value for t where we stop the approximation
    rkf54 = RungeKuttaFehlberg54(F,4,h,tol) #function, dimension, stepsize, tolerance. Why is dimension = 4?

    while(W[0]<tEnd): #Why use W[0] here? Is it the current time?
        W , E = rkf54.safeStep(W); #Returns a new vector of approximated y-values W. And an error vector?

    rkf54.setStepLength(tEnd-W[0]); #Makes no sense. This would make the step length 1.9 after the first iteration.
    W,E = rkf54.step(W); #Make one step forward.

    print(W,E);

if __name__ == "__main__":
    # execute only if run as a script
    main()
