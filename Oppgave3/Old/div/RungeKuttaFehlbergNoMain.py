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
             Win, velocitiesIn):
        s=np.zeros((6,self.dim))
        velocitiesOut=np.zeros((6,3)) #vx, vy, vz
        print("v:")
        print(velocitiesOut)

        #Calculate all 6 s-values/vectors
        for i in range(0,6):
            s[i,:]=self.F(Win+self.h*self.A[i,0:i].dot(s[0:i,:]),velocitiesIn)
            print(s[i,:])

        Zout=Win+self.h*(self.B[0,:].dot(s)); #In the book, Zout is better than Wout, given that Zout is the locally extrapolated version.
        Wout=Win+self.h*(self.B[1,:].dot(s));

        E=np.linalg.norm(Wout-Zout,2)/np.linalg.norm(Wout,2);

        velocitiesOut = velocitiesIn
        return Wout, E, velocitiesOut

    def safeStep(self,
                 Win, velocityIn): #Win are the current approximations for y
        Wout,E, velocitiesOut = self.step(Win, velocityIn); #Makes a normal step
        # Check if the error is tolerable
        if(not self.isErrorTolerated(E)):
            #Try to adjust the optimal step length
            self.adjustStep(E);
            Wout,E, velocitiesOut = self.step(Win, velocityIn);
        # If the error is still not tolerable
        counter=0;
        while(not self.isErrorTolerated(E)):
            #Try if dividing the steplength with 2 helps.
            self.divideStepByTwo();

            Wout,E, velocitiesOut = self.step(Win, velocityIn);
            counter = counter + 1;
            if(counter>10):
                sys.exit(-1); #Some sort of way to get out of an inifinite step size reducing loop

        self.adjustStep(E);

        return Wout, E, velocitiesOut

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
