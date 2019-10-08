
from RungeKuttaFehlberg_3 import RungeKuttaFehlberg54
from RungeKuttaFehlberg_3 import F
import time

import numpy as np

import matplotlib.pyplot as plot
import matplotlib.animation as animation

""" class NumericalSolver:

    @staticmethod
    def step(planets, h):
        "Trapezoid method

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


 """
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
                 G=6.67e-11,
                 m1=7.3477e22,
                 m2=5.97e24):
        self.GravConst = G
        self.mSol = m2 #Mass of sun
        self.mPlanet = m1 #Mass of planet
        self.state = np.asarray(init_state, dtype='float')
        self.xy = [[0], [0.4055e9]]

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
        U=-G*m1*m2/ np.sqrt(x**2+y**2) #Potential energy?
        K= m1*(vx**2+vy**2)/2 #Kinetic energy?
        return K+U #Total energy

    def time_elapsed(self):
        return self.state[0] #t

    def ydot(self, x):
        G=self.GravConst
        m2=self.mSol
        Gm2=G*m2;

        px2=0;py2=0; #Position of the sun. Needs to be changed in two-planet problem!
        px1=x[1]; #x0
        py1=x[3]; #y0
        vx1=x[2]; #vx0
        vy1=x[4]; #vy0
        dist=np.sqrt((px2-px1)**2+(py2-py1)**2);
        z=np.zeros(5);

        #Equation from the book for one-body-problem
        #Gives the change in velocity
        #t0,x0,vx0,y0,vx0
        z[0]=1
        z[1]=vx1 #x' = vx
        z[2]=(Gm2*(px2-px1))/(dist**3) # (vx)' = (Change in velocity in x-direction per time)
        z[3]=vy1 #y' = vy
        z[4]=(Gm2*(py2-py1))/(dist**3) #(vy)' = (Change in velocity in y-direction per time)

        self.xy[0].append(self.position()[0])
        self.xy[1].append(self.position()[1])
        return z #Returns s1 or s2


# make an Orbit instance

#List of all planets in the system
# [t0,x0,vx0,y0,vx0]
list = [Orbit([0.0,0.0, 1082, 362570e3, 0.0]),Orbit([0.0,0.0, 0.0, 0.0, 0.0])]
planetz = np.array(list)

#h=0.1; #Step size
tol=05e-14; #RelativeS error, or just error?
#tEnd=10.0; #Value for t where we stop the approximation
#dt = 1./30 # 30 frames per second
dt = 1.0/30.0 # 30 frames per second
rkf54 = RungeKuttaFehlberg54(planetz[0].ydot,planetz[0].state.size,dt,tol) #function, dimension, stepsize, tolerance.


# The figure is set
fig = plot.figure() # matplotlib.pyplot = plot
axes = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                     xlim=(-4.5e8, 4.5e8), ylim=(-4.5e8, 4.5e8))

line1, = axes.plot([], [], 'o-g', lw=2) # A green planet
line1_2, = axes.plot([], [], 'r--', linewidth=0.5)
line2, = axes.plot([], [], 'o-y', lw=2) # A yellow sun
time_text = axes.text(0.02, 0.95, '', transform=axes.transAxes)
position_text = axes.text(0.02, 0.85, '', transform=axes.transAxes)

def init():
    """initialize animation""" #Just empty values
    line1.set_data([], [])
    line1_2.set_data([], [])
    line2.set_data([], [])
    time_text.set_text('')
    position_text.set_text('')
    return line1, line1_2, line2, time_text, position_text


#Basically main
def animate(i):
    """perform animation step"""
    global orbit, orbit2, dt, planetz #Allows to modify a variable outside of the current scope

    W , E = rkf54.safeStep(planetz[0].state);
    #print(W)
    planetz[0].state = W

    #NumericalSolver.step(planetz,dt)

    line1.set_data(*planetz[0].position()) #Green planet * operator means to take the array data x and y, and use them as parameters in the set_data function.
    line1_2.set_data(planetz[0].xy)
    line2.set_data([0,0]) #Sun. Position is constant.
    time_text.set_text('time = %.1f' % planetz[0].time_elapsed())
    position_text.set_text('avstand = %.3f* km' % ((np.sqrt(planetz[0].position()[0]**2 + planetz[0].position()[1]**2))/1e3))
    return line1, line1_2, line2, time_text, position_text

# choose the interval based on dt and the time to animate one step
# Take the time for one call of the animate.
t0 = time.time()
animate(0)
t1 = time.time()

print((t1 - t0))

delay = 1000 * dt - (t1 - t0)
print(delay)

anim=animation.FuncAnimation(fig,        # figure to plot in
                        animate,    # function that is called on each frame
                        frames=1000, # total number of frames
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
