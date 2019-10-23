
from RungeKuttaFehlberg import RungeKuttaFehlberg54
from SaturnV import SaturnV
import time

import numpy as np
import math

import matplotlib.pyplot as plot
import matplotlib.animation as animation

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

class Orbit:

    def __init__(self,
                 init_state = [0, 0, 1, 2, 0],
                 G=6.67e-11,
                 m=5.97e24):
        self.GravConst = G
        self.m = m #Mass of sun
        self.state = np.asarray(init_state, dtype='float')
        self.xy = [[0], [(12756.28/2) * 1000]]
        self.acceleration = 0
        self.angle1 = math.pi/2
        self.angle2 = 0
        self.angle = 1.30

    def position(self):
        x = self.state[1]
        y = self.state[2]
        return (x, y)

    def moh(self):
        x = self.state[1]
        y = self.state[2]
        earth_radius = (12756.28/2) * 1000
        distance_from_earth_center = np.sqrt(x*x + y*y)

        return distance_from_earth_center - earth_radius

    def time_elapsed(self):
        return self.state[0]

    def get_velocity(self):
        return [self.state[3], self.state[4]]

    def get_acceleration(self):
        return self.acceleration

    def get_position(self):
        return [self.state[1], self.state[2]]

    def ydot(self, x):
        G = self.GravConst
        Gm = G * self.m
        t = x[0]
        px2 = 0
        py2 = 0
        px1 = x[1]
        py1 = x[2]
        vx1 = x[3]
        vy1 = x[4]
        dist = np.sqrt((px2 - px1) ** 2 + (py2 - py1) ** 2)

        # Force from gravity on rocket divided by rocket mass
        Fg_x = (Gm * (px2 - px1)) / (dist ** 3)
        Fg_y = (Gm * (py2 - py1)) / (dist ** 3)

        # Force from air drag on rocket divided by rocket mass
        absolute_velocity = np.sqrt(vx1*vx1 + vy1*vy1)
        saturnV = SaturnV(m1,c1,b1,t1,m2,c2,b2,t2,m3,c3,b3,t3)
        Fd = self.get_air_drag(self.moh(), 0.2, self.get_area(t), absolute_velocity) / saturnV.calculateMass(t)

        F = saturnV.calculateThrust(t)

        print(math.cos(self.angle), F*math.cos(self.angle) + Fg_x - Fd*math.cos(self.angle))

        self.acceleration = (F + Fg_y - Fd)/saturnV.calculateMass(t) #Merk fortegnene inne i ligningen
        z = np.zeros(5)
        z[0] = 1
        z[1] = vx1
        z[2] = vy1
        z[3] = (F*math.cos(self.angle) + Fg_x - Fd*math.cos(self.angle))/saturnV.calculateMass(t)
        z[4] = (F*math.sin(self.angle) + Fg_y - Fd*math.sin(self.angle))/saturnV.calculateMass(t) #Merk fortegnene inne i ligningen

        self.xy[0].append(self.get_position()[0])
        self.xy[1].append(self.get_position()[1])

        return z

    def get_air_drag(self, h, Cd, A, v):
        return 1/2 * Cd * self.get_air_resistance(h) * A * v**2


    def get_air_resistance(self, h):
        return self.get_air_pressure(h) / self.get_temperature(h) * 3.4855

    def get_area(self, t):
        time_stage_one = b1
        time_stage_two = b1 + b2

        if t < time_stage_one:
            return 113
        elif t >= time_stage_one and t < time_stage_two:
            return 80.1
        else:
            return 34.2

    # Returns the air pressure in Pascal based on height above sea level
    def get_air_pressure(self, h):
        if 0 <= h <= 11000:  # Troposphere
            return 101.29e3 * (self.get_temperature(h)/288.08)**5.256

        elif 11000 < h <= 25000:  # Lower stratosphere
            return 127.76e3 * np.e**(-0.000157 * h)

        elif 25000 < h:  # Higher stratosphere
            return (2.488e3 * (self.get_temperature(h)/216.6)**(-11.388))


    # Returns the air temperature in Kelvin based on height above sea level
    def get_temperature(self, h):

        if 0 <= h <= 11000:  # Troposphere
            return 288.19 - 0.00649 * h

        elif 11000 < h <= 25000:  # Lower stratosphere
            return 216.69

        elif 25000 < h:  # Higher stratosphere
            return (141.94 + 0.00299 * h)


#List of all planets in the system
list = [Orbit([0.0, 0.0,  (12756.28/2) * 1000, 0, 0.0]),Orbit([0.0,0.0, 0.0, 0.0, 0.0])]
planetz = np.array(list)


tol=02e-14
animation_time = 0
time_0 = 0
time_difference = 0


dt = 1.0/30.0 # 30 frames per second
rkf54 = RungeKuttaFehlberg54(planetz[0].ydot, planetz[0].state.size, dt, tol) #function, dimension, stepsize, tolerance.


plotScale = (12756.28/2) * 1000 # meters
# The figure is set
fig = plot.figure() # matplotlib.pyplot = plot

axes = fig.add_subplot(111, aspect='equal', autoscale_on=False, xlim=(-plotScale, plotScale), ylim=(-plotScale, plotScale * 1.5))


earth = plot.Circle((0, 0), (12756.28/2) * 1000, color='blue', alpha=0.2)
axes.add_artist(earth)


line1, = axes.plot([], [], 'o-g', linewidth=1, markersize=1)
line1_2, = axes.plot([], [], 'r--', linewidth=0.5)
line2, = axes.plot([], [], 'o-y', lw=2)
time_text = axes.text(0.02, 0.95, '', transform=axes.transAxes)
position_text = axes.text(0.02, 0.90, '', transform=axes.transAxes)
acceleration_text = axes.text(0.02, 0.80, '', transform=axes.transAxes)

def init():
    """initialize animation"""
    line1.set_data([], [])
    line1_2.set_data([], [])
    line2.set_data([], [])
    time_text.set_text('')
    position_text.set_text('')
    acceleration_text.set_text('')
    return line1, line1_2, line2, time_text, position_text, acceleration_text



def animate(i):
    """perform animation step"""
    global dt, planetz, time_0, time_difference
    t0 = time_0
    time_1 = planetz[0].state[0]
    time_0 = time_1
    time_difference = time_1 - t0
    time_to_sleep = time_difference / dt - 1
    """
    if time_to_sleep > 0:
        time.sleep(time_to_sleep * dt)
    """


    W , E = rkf54.safeStep(planetz[0].state)
    planetz[0].angle -= 0.0093
    print('t')

    planetz[0].state = W

    line1.set_data(*planetz[0].position())
    line1_2.set_data(planetz[0].xy)
    line2.set_data([0, 406356640]) # moon position
    time_text.set_text('time = %.1f' % planetz[0].time_elapsed())
    acceleration_text.set_text('Acceleration = %.3f m/s^2' % planetz[0].get_acceleration())
    position_text.set_text('x = %.3f km' % (planetz[0].get_position()[0]/1e3) + ', y = %.3f km' % ((planetz[0].get_position()[1] - (12756.28/2) * 1000)/1e3))
    return line1, line1_2, line2, time_text, position_text, acceleration_text


t0 = time.time()
animate(0)
time_1 = time.time()

delay = 100 * dt - (time_1 - t0)

anim=animation.FuncAnimation(fig,        # figure to plot in
                        animate,    # function that is called on each frame
                        frames=15000, # total number of frames
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
