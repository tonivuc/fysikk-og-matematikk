#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Oppgave 2

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

def F(Y):
    # The system of equations that we want to solve
    res = np.zeros(3)
    res[0] = 1 # time
    res[1] = Y[1] + Y[2] #Y1'
    res[2] = -Y[1] + Y[2] #Y2'
    return res

def main():
    W=np.array([0, 1, 0]) #Initial values
    h=0.25
    tol=02e-14
    tEnd=1.0

    rkf54 = RungeKuttaFehlberg54(F,3,h,tol)

    while(W[0]<tEnd):
        W , E = rkf54.safeStep(W);

    rkf54.setStepLength(tEnd-W[0]);
    W,E = rkf54.step(W);

    y = correctSolution(tEnd);
    print(W,E);
    print('Global feil i y_1: ', globalError(W[1], y[1]))
    print('Global feil i y_2: ', globalError(W[2], y[2]))
    print('Akkumulert lokal error: ', rkf54.localError)

    antTester = 10;
    beregningstid = np.zeros(antTester);
    toleranser = np.zeros(antTester);
    for i in range(0, antTester):
        W=np.array([0, 1, 0]);
        rkf54 = RungeKuttaFehlberg54(F,3,h,tol)
        start_time = time.time()
        while(W[0]<tEnd):
            W , E = rkf54.safeStep(W);

        rkf54.setStepLength(tEnd-W[0]);
        W,E = rkf54.step(W);
        stop_time = time.time()
        beregningstid[i] = stop_time-start_time;
        toleranser[i] = tol;
        print("tid: ", stop_time-start_time);

        tol *= 5;

    plot.loglog(toleranser, beregningstid)
    plot.xlabel("Toleranse")
    plot.ylabel("Modeltid/beregningstid")
    plot.title("Funksjon av toleranse")
    plot.grid(True)
    plot.show()

def globalError(w, y) :
    return abs(w-y);

def correctSolution(tEnd) :
    y = np.zeros(3);
    y[1] = m.pow(m.e, tEnd)*m.cos(tEnd);
    y[2] = -m.pow(m.e, tEnd)*m.sin(tEnd);
    return y;

if __name__ == "__main__":
    # execute only if run as a script
    main()
