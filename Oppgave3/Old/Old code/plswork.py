
from numpy import sqrt
from RungeKuttaFehlberg_2 import RungeKuttaFehlberg54
from RungeKuttaFehlberg_2 import F
import time

import numpy as np
import scipy.integrate as integrate #Not in use?

import matplotlib.pyplot as plot
import matplotlib.animation as animation

class NumericalSolver:

    @staticmethod
    def step(planets, h):
        """Trapezoid method"""

        i = 0
        while i < 1: #Once upon a time this was a two-body problem
            x=planets[i].state
            s1=NumericalSolver.ydot(planets[i],x)
            s2=NumericalSolver.ydot(planets[i],x+h*s1)
            w = x+h*(s1+s2)/2 #Vectors being multiplied and added together. Gives vector with position and velocity
            planets[i].state = w

            i = i + 1
            if(i==len(planets)):
                break

    @staticmethod
    def ydot(self,x): #y derivert
        G=self.GravConst
        m2=self.mSol
        Gm2=G*m2;

        px2=0;py2=0; #Position of the sun. Needs to be changed in two-planet problem!
        px1=x[1]; #x0
        py1=x[3]; #y0
        vx1=x[2]; #vx0
        vy1=x[4]; #vy0
        dist=sqrt((px2-px1)**2+(py2-py1)**2);
        z=np.zeros(5);

        #Equation from the book for one-body-problem
        #Gives the change in velocity
        #t0,x0,vx0,y0,vx0
        z[0]=1
        z[1]=vx1 #x' = vx
        z[2]=(Gm2*(px2-px1))/(dist**3) # (vx)' = (Change in velocity in x-direction per time)
        z[3]=vy1 #y' = vy
        z[4]=(Gm2*(py2-py1))/(dist**3) #(vy)' = (Change in velocity in y-direction per time)
        return z #Returns s1 or s2



class Orbit:
    """

    Orbit Class

    Original init-state: init_state = [0, 0, 1, 2, 0],

    init_state is [t0,x0,vx0,y0,vx0],
    where (x0,y0) is the initial position
    , (vx0,vy0) is the initial velocity
    and t0 is the initial time
    """
    def __init__(self,       #t0,x0,vx0,y0,vx0
                 init_state = [0, 0, 1, 2, 0],
                 G=1,
                 m1=1,
                 m2=3):
        self.GravConst = G
        self.mSol = m2 #Mass of sun
        self.mPlanet = m1 #Mass of planet
        self.state = np.asarray(init_state, dtype='float')

    def position(self):
        """compute the current x,y positions of the pendulum arms"""
        x = self.state[1] #x0
        y = self.state[3] #vx0, starting velocity in x direction
        return (x, y)

    def energy(self):
        x = self.state[1] #x0
        y = self.state[3] #y0
        vx = self.state[2] #vx0
        vy = self.state[4] #xy0
        m1 = self.mPlanet
        m2 = self.mSol
        G = self.GravConst
        U=-G*m1*m2/sqrt(x**2+y**2) #Potential energy?
        K= m1*(vx**2+vy**2)/2 #Kinetic energy?
        return K+U #Total energy

    def time_elapsed(self):
        return self.state[0] #t



# make an Orbit instance

#List of all planets in the system
# [t0,x0,vx0,y0,vx0]
list = [Orbit([0.0,0.0, 1.0, 2.0, 0.0],1,1,3),Orbit([0.0,0.0, 0.0, 0.0, 0.0],1,3,1)]
planetz = np.array(list)

#h=0.1; #Step size
tol=05e-14; #RelativeS error, or just error?
tEnd=10.0; #Value for t where we stop the approximation
#dt = 1./30 # 30 frames per second
dt = 1./90 # 30 frames per second
rkf54 = RungeKuttaFehlberg54(F,5,dt,tol) #function, dimension, stepsize, tolerance. Why is dimension = 4?


# The figure is set
fig = plot.figure() # matplotlib.pyplot = plot
axes = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                     xlim=(-3, 3), ylim=(-3, 3))

line1, = axes.plot([], [], 'o-g', lw=2) # A green planet
line2, = axes.plot([], [], 'o-y', lw=2) # A yellow sun
time_text = axes.text(0.02, 0.95, '', transform=axes.transAxes)
energy_text = axes.text(0.02, 0.90, '', transform=axes.transAxes)

def init():
    """initialize animation""" #Just empty values
    line1.set_data([], [])
    line2.set_data([], [])
    time_text.set_text('')
    energy_text.set_text('')
    return line1, line2, time_text, energy_text


#Basically main
def animate(i):
    """perform animation step"""
    global orbit, orbit2, dt, planetz #Allows to modify a variable outside of the current scope

    W , E = rkf54.safeStep(planetz[0].state);
    planetz[0].state = W

    #NumericalSolver.step(planetz,dt)

    line1.set_data(*planetz[0].position()) #Green planet * operator means to take the array data x and y, and use them as parameters in the set_data function.
    line2.set_data([0,0]) #Sun. Position is constant.
    time_text.set_text('time = %.1f' % planetz[0].time_elapsed())
    energy_text.set_text('energy = %.3f J' % planetz[0].energy())
    return line1,line2, time_text, energy_text

# choose the interval based on dt and the time to animate one step
# Take the time for one call of the animate.
t0 = time.time()
animate(0)
t1 = time.time()

delay = 1000 * dt - (t1 - t0)

anim=animation.FuncAnimation(fig,        # figure to plot in
                        animate,    # function that is called on each frame
                        frames=300, # total number of frames
                        interval=delay, # time to wait between each frame.
                        repeat=False,
                        blit=True,
                        init_func=init # initialization
                        )

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
#anim.save('orbit.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plot.show()
