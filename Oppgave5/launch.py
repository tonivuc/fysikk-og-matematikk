import numpy as np
import math
from SaturnV import SaturnV
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from PIL import Image
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

earth_mass = 5.9736 * (10**24) # Jordens masse i kilogram
earth_radius = (12756.28/2) * 1000 # Jordens radius ved ekvator i meter
gravity_constant = 6.67428 * (10**(-11)) # Gravitasjonskonstanten
area_stage_one = 113 # Overflateområde steg en for Saturn V
area_stage_two = 80.1 # Overflateområde steg to for Saturn V
area_stage_three = 34.2 # Overflateområde steg en tre Saturn V

#Mass of different stages (kilograms)
m1 = 2970000
m2 = 680000
m3 = 183000
#Burn time of different stages (seconds)
b1 = 168
b2 = 360
b3 = 165
#Fuel consumption of different stages (kilograms per second)
c1 = 12857.1
c2 = 1266.9
c3 = 219
#Thrust of different stages (Newtons)
t1 = 35100000
t2 = 5141000
t3 = 1033100

#Create an object to use:
saturnV = SaturnV(m1,c1,b1,t1,m2,c2,b2,t2,m3,c3,b3,t3)

time_stage_one = b1
time_stage_two = b1 + b2
time_stage_three = b1 + b2 + b3

def height(t, v):
    return (1/2)*(v[t] + v[t - 1])*t

def distance(t, v):
    return (1/2)*(v[t] + v[t - 1])*t + earth_radius

def mass(t):
    return saturnV.calculateMass(t)

def Fs(t):
    return saturnV.calculateThrust(t)

def Fg(t, v):
    h = height(t, v)
    return gravity_constant * ((earth_mass * mass(t))/(distance(t, v)**2))

def pressure(h):
    if h < 11000:
        return (101.29) * ((temperature(h)/288.08)**5.256)
    elif h >= 11000 and h < 25000:
        return 127.76*(math.exp(-0.000157*h))
    else:
        return 2.488*((temperature(h)/216.6)**-11.388)

def temperature(h):
    if h < 11000:
        return 288.19 - (0.00649)*h
    elif h >= 11000 and h < 25000:
        return 216.69
    else:
        return 141.94 + (0.00299)*h

def area(t):
    if t < time_stage_one:
        return area_stage_one
    elif t >= time_stage_one and t < time_stage_two:
        return area_stage_two
    else:
        return area_stage_three

def speed(t, a, v):
    return a[t - 1] + v[t - 1]

def density(t, v):
    h = height(t, v)
    return (pressure(h)/temperature(h)) * 3.4855

def Fd(t, a, v):
    return (1/2)*(1/2)*density(t, v)*area(t)*(speed(t, a, v)**2)


def imscatter(x, y, image, ax=None, zoom=1):
    if ax is None:
        ax = plt.gca()
    try:
        image = plt.imread(image)
    except TypeError:
        # Likely already an array...
        pass
    im = OffsetImage(image, zoom=zoom)
    x, y = np.atleast_1d(x, y)
    artists = []
    for x0, y0 in zip(x, y):
        ab = AnnotationBbox(im, (x0, y0), xycoords='data', frameon=False)
        artists.append(ax.add_artist(ab))
    ax.update_datalim(np.column_stack([x, y]))
    ax.autoscale()
    return artists


def main():
    # Antall tidstrinn
    steps = b1 + b2 + b3

    # Initialiserer startverdier
    h = 0.001 # Initialiseres med en verdi større enn null for å komme inn i løkken. Verdien påvirker ingen av kalkulasjonene.
    F = []
    a = []
    m = []
    v = []
    m.append(mass(0))
    v.append(0)
    F.append(Fs(0) - gravity_constant * ((earth_mass * m[0])/(earth_radius**2)))
    a.append(F[0] / m[0])
    t = 1
    img_path_rocket_down = '/home/vebovs/Desktop/fysikk-og-matematikk/Oppgave5/rocket_down.png'
    img_path_rocket_up_no_fuel = '/home/vebovs/Desktop/fysikk-og-matematikk/Oppgave5/rocket_up_no_fuel.png'
    img_path_rocket_up = '/home/vebovs/Desktop/fysikk-og-matematikk/Oppgave5/rocket_up.png'
    img_path = img_path_rocket_up
    down = False
    while h >= 0:
        last_h = h
        v.append(speed(t, a, v))
        F.append(Fs(t) - (Fg(t, v) + Fd(t, a, v)))
        m.append(mass(t))
        a.append(F[t] / m[t])
        h = height(t, v)
        plt.clf()
        plt.xlim(-500, 500)
        plt.ylim(0, 3000000)
        plt.ylabel('Høyde (m)')
        if t > time_stage_three and not down:
            img_path = img_path_rocket_up_no_fuel
        if last_h > h and not down:
            img_path = img_path_rocket_down
            down = True
        img = OffsetImage(plt.imread(img_path))
        img = mpimg.imread(img_path)
        plt.imshow(img, aspect = 'auto', extent = [-40, 40, h - 100000, h + 200000])
        plt.pause(0.000000001)
        t = t + 1

    plt.show()

if __name__ == "__main__":
    main()