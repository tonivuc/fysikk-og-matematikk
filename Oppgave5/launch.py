import numpy as np
import math
from SaturnV import SaturnV
import matplotlib.pyplot as plt

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

def height(t, v):
    if(t == 169):
        print("speed: ", v[t])
        print("speed0: ", v[t - 1])
        print("res: ", (1/2)*(v[t] + v[t - 1])*t)
    return (1/2)*(v[t] + v[t - 1])*t

def distance(t, v):
    return (1/2)*(v[t] + v[t - 1])*t + earth_radius

def mass(t):
    return saturnV.calculateMass(t)

def Fs(t):
    return saturnV.calculateThrust(t)

def Fg(t, v):
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
    time_stage_one = b1
    time_stage_two = b1 + b2

    if t < time_stage_one:
        return area_stage_one
    elif t >= time_stage_one and t < time_stage_two:
        return area_stage_two
    else:
        return area_stage_three

def speed(t, a, v):#live on the edge monjaro is for fucks with no xXxskillzxXx mlg 1337
    return a[t - 1] + v[t - 1]

def density(t, v):
    h = height(t, v)
    if(t==169):
        print("height: ", h)
        print("temp: ", temperature(h))
        print("pressure: ", pressure(h))
    return (pressure(h)/temperature(h)) * 3.4855

def Fd(t, a, v):
    return (1/2)*(1/2)*density(t, v)*area(t)*(speed(t, a, v)**2)


def main():
    # Lav jordbane går cirka til 2000 kilometer
    target_height = 2000000 # Målhøyde i meter
    # Antall tidstrinn
    steps = b1 + b2 + b3

    # Lager tabellene med plassholdere
    F = np.zeros(steps + 1)
    a = np.zeros(steps + 1)
    m = np.zeros(steps + 1)
    v = np.zeros(steps + 1)
    arr = np.zeros(steps + 1)
    arr[0] = 0
    # Initialiserer startverdier
    m[0] = mass(0)
    v[0] = 0
    # Tyngdekraften er betydelig større enn skyvekraften
    F[0] = Fs(0) - gravity_constant * ((earth_mass * m[0])/(earth_radius**2))
    a[0] = F[0] / m[0]

    for t in range(1, steps + 1):
        arr[t] = t
        v[t] = speed(t, a, v)
        F[t] = Fs(t) - (Fg(t, v) + Fd(t, a, v))
        m[t] = mass(t)
        print("Second: ", t)
        print("Air: ", Fd(t, a, v))
        print("Gravity: ", Fg(t, v))
        print("Forces: ", F[t])
        print("Mass: ", m[t])
        a[t] = F[t] / m[t]
        print("accel: ", a[t])

    print("ac: ", a[99])
    print("speed: ", v[70])
    print(height(steps, v))
    print("Height in km: ", height(50, v)/1000)
    print("DENSITY: ", (pressure(1500)/temperature(1500)) * 3.4855)
    plt.scatter(arr, a)
    plt.show()

if __name__ == "__main__":
    main()