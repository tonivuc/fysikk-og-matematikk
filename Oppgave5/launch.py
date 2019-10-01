import numpy as np
from SaturnV import SaturnV

earth_mass = 5.9736 * (10**24) # Jordens masse i kilogram
earth_radius = (12756.28/2) * 1000 # Jordens radius ved ekvator
gravity_constant = 6.67428 * (10**(-11)) # Gravitasjonskonstanten
area_stage_one = 113 # Overflateområde steg en for Saturn V
area_stage_two = 79.5 # Overflateområde steg to for Saturn V
area_stage_three = 34.4 # Overflateområde steg en tre Saturn V

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

#Create an object to use:
saturnV = SaturnV(m1,c1,b1,t1,m2,c2,b2,t2,m3,c3,b3,t3);

def height(t, a):
    return (1/2)*a*(t*t)

def distance(t, a):
    return (1/2)*a*(t*t) + earth_radius

def mass(t):
    return saturnV.calculateMass(t)

def Fs(t):
    return saturnV.calculateThrust(t)

def Fg(t, a):
    return gravity_constant * ((earth_mass * mass(t))/distance(t, a))

def speed(t, a, v):
    return t * a[t] * v[t - 1]

def density(t, a):
    h = height(t, a)


def Fd(t, a, v):
    return (1/2)*(1/2)*density(t, a)*area*(speed(t, a, v)*speed(t, a, v))


def main():
    # Lav jordbane går cirka til 2000 kilometer
    target_height = 2000000 # Målhøyde i meter
    # Antall tidstrinn
    steps = 10000000

    # Lager tabellene med plassholdere
    F = np.zeros(steps)
    a = np.zeros(steps)
    m = np.zeros(steps)
    v = np.zeros(steps)

    # Initialiserer startverdier
    v[0] = 0
    F[0] = Fs(0) - Fg(0, 0)
    a[0] = F[0] / m[0]

if __name__ == "__main__":
    main()