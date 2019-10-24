
from RungeKuttaFehlberg import RungeKuttaFehlberg54
import time

import numpy as np

import matplotlib.pyplot as plot
import matplotlib.animation as animation

class Orbit:
    """
    Oppgave 3

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
        self.mEarth = m2 #Mass of earth
        self.mMoon = m1 #Mass of moon
        self.state = np.asarray(init_state, dtype='float')
        self.xy = [[self.state[1]], [self.state[3]]]

    def position(self):
        """compute the current x,y positions of the pendulum arms"""
        x = self.state[1] #position in x-direction
        y = self.state[3] #position in y-direction
        return (x, y)

    def energy(self):
        x = self.state[1]
        y = self.state[3]
        vx = self.state[2]
        vy = self.state[4]
        m1 = self.mMoon
        m2 = self.mEarth
        G = self.GravConst
        U=-G*m1*m2/ np.sqrt(x**2+y**2) #Calculates potential energy
        K= m1*(vx**2+vy**2)/2 #Calculates kinetic energy
        return K+U #Total energy

    def time_elapsed(self):
        return self.state[0] #time

    def ydot(self, x):
        G=self.GravConst
        m2=self.mEarth
        Gm2=G*m2

        px2=0
        py2=0 #Position of the centre of earth.
        px1=x[1]
        py1=x[3]
        vx1=x[2]
        vy1=x[4]
        dist=np.sqrt((px2-px1)**2+(py2-py1)**2) #Calculates distance from the moon to the earth
        z=np.zeros(5);

        z[0]=1
        z[1]=vx1 #x' = vx
        z[2]=(Gm2*(px2-px1))/(dist**3) # (vx)' = (Change in velocity in x-direction per time)
        z[3]=vy1 #y' = vy
        z[4]=(Gm2*(py2-py1))/(dist**3) #(vy)' = (Change in velocity in y-direction per time)

        self.xy[0].append(self.position()[0])
        self.xy[1].append(self.position()[1])
        return z


#List of all planets in the system
list = [Orbit([0.0,0.0, 1.082e3, 362570e3, 0.0]),Orbit([0.0,0.0, 0.0, 0.0, 0.0])]
planetz = np.array(list)

tol=02e-14; #The error tolerance
dt = 1.0/30.0 # 30 frames per second
rkf54 = RungeKuttaFehlberg54(planetz[0].ydot,planetz[0].state.size,dt,tol) #function, dimension, stepsize, tolerance.


# The figure is set
fig = plot.figure()
axes = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                     xlim=(-4.5e8, 4.5e8), ylim=(-4.5e8, 4.5e8))

line1, = axes.plot([], [], 'o-g', lw=2) # A green moon
line1_2, = axes.plot([], [], 'r--', linewidth=0.5)
line2, = axes.plot([], [], 'o-y', lw=2) # A yellow earth
time_text = axes.text(0.02, 0.95, '', transform=axes.transAxes)
position_text = axes.text(0.02, 0.85, '', transform=axes.transAxes)

def init():
    """initialize animation"""
    line1.set_data([], [])
    line1_2.set_data([], [])
    line2.set_data([], [])
    time_text.set_text('')
    position_text.set_text('')
    return line1, line1_2, line2, time_text, position_text


#Basically main
def animate(i):
    """perform animation step"""
    global dt, planetz #Allows to modify a variable outside of the current scope

    W , E = rkf54.safeStep(planetz[0].state)

    planetz[0].state = W

    if ((W[1] > -1000000) and (W[1] < 1000000)):
        print(W[1])
        print(W[3])

    line1.set_data(*planetz[0].position()) #The moons position, * operator means to take the array data x and y, and use them as parameters in the set_data function.
    line1_2.set_data(planetz[0].xy) #Finds the line where the moon has been.
    line2.set_data([0,0]) #Earth. Position is constant.
    time_text.set_text('time = %.1f' % planetz[0].time_elapsed())
    position_text.set_text('avstand = %.3f* km' % ((np.sqrt(planetz[0].position()[0]**2 + planetz[0].position()[1]**2))/1e3))
    return line1, line1_2, line2, time_text, position_text

# Take the time for one call of the animate.
t0 = time.time()
animate(0)
t1 = time.time()

delay = 1000 * dt - (t1 - t0)

anim=animation.FuncAnimation(fig,
                        animate,
                        frames=1500,
                        interval=delay,
                        repeat=False,
                        blit=True,
                        init_func=init
                        )

plot.show()
