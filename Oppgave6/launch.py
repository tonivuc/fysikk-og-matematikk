import math
from SaturnV import SaturnV
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os

earth_mass = 5.9736 * (10**24) # Jordens masse i kilogram
earth_radius = (12756.28/2) * 1000 # Jordens radius ved ekvator i meter
gravity_constant = 6.67428 * (10**(-11)) # Gravitasjonskonstanten
area_stage_one = 113 # Overflateområde steg en for Saturn V
area_stage_two = 80.1 # Overflateområde steg to for Saturn V
area_stage_three = 34.2 # Overflateområde steg en tre Saturn V
num_height = 0 # Starthøyde
num_distance = earth_radius # Avstand fra sentrum av jordkloden

#Burn time of different stages (seconds)
b1 = 168
b2 = 360
b3 = 165
#Mass of different stages (kilograms)
m1 = 130000
m2 = 40100
m3 = 13500
#Fuel in different stages (kilograms)
d1 = 2160000
d2 = 456100
d3 = 109500
#Fuel consumption of different stages (kilograms per second)
c1 = 12857.1
c2 = 1266.9
c3 = 219
#Thrust of different stages (ts1 = thrust-sealevel-stage1, ts2 = thrust-vacuum-stage1) in newton
ts1 = 33850000
tv1 = 38850000
ts2 = 2431000
tv2 = 5165500
ts3 = 486200
tv3 = 1033100

#Create an object to use:
saturnV = SaturnV(m1,c1,d1,ts1,tv1,m2,c2,d2,ts2,tv2,m3,c3,d3,ts3,tv3)

# Tiden hver motoravkopling oppstår (sekund)
time_stage_one = b1
time_stage_two = b1 + b2
time_stage_three = b1 + b2 + b3

# Regner ut høyde traversert etter ett tidstrinn utfra endring i fart og legger det til summen av av alle tidstrinnene
# Brukes til å posisjonere raketten og utregning av luftmotstand
# Målenhet: Meter (m)
def height(t, v):
    if(v == 0):
        return 0
    global num_height
    num_height = num_height + ((1/2)*(v[t] + v[t - 1]))
    return num_height

# Samme som utregning av høyden, men legger til radiusen for jordkloden
# Brukes til å regne ut avstanden mellom jordkloden og raketten
# Målenhet: Meter (m)
def distance(t, v):
    global num_distance
    num_distance = num_distance + ((1/2)*(v[t] + v[t - 1]))
    return num_distance

# Returnerer massen til raketten gitt et tidstrinn
# Målenhet: Kilogram (kg)
def mass(t):
    return saturnV.calculateMass(t)

# Returnerer skyvekraften til raketten gitt et tidstrinn
# Målenhet: Newton (N)
def Fs(t, v):
    return saturnV.calculateThrust(t, pressure(height(t, v))/100)

# Analytisk utregning av tyngdekraften mellom to legem
# Brukes til å regne ut tyngdekraften mellom jordkloden og raketten
# Målenhet: Newton (N)
def Fg(t, v):
    return gravity_constant * ((earth_mass * mass(t))/(distance(t, v)**2))

# Returnerer det atmosfæriske trykket gitt en høyde
# Målenhet: Pascal (Pa)
def pressure(h):
    if h < 11000:
        return (101.29)*((temperature(h)/288.08)**5.256)
    elif h >= 11000 and h < 25000:
        return (127.76)*(math.exp(-0.000157*h))
    else:
        return (2.488)*((temperature(h)/216.6)**-11.388)

# Returnerer temperaturen i atmosfæren gitt en høyde
# Målenhet: Kelvin (K)
def temperature(h):
    if h < 11000:
        return 288.19 - (0.00649)*h
    elif h >= 11000 and h < 25000:
        return 216.69
    else:
        return 141.94 + (0.00299)*h

# Returnerer overflatesområdet til raketten under de forskjellige trinnene
# Målenhet: Kvadratmeter (m²)
def area(t):
    if t < time_stage_one:
        return area_stage_one
    elif t >= time_stage_one and t < time_stage_two:
        return area_stage_two
    else:
        return area_stage_three

# Regner ut fart i ett tidstrinn utfra farten og akselerasjonen under forrige tidstrinn
# Målenhet: Meter per sekund (m/s)
def speed(t, a, v):
    return a[t - 1] + v[t - 1]

# Regner ut tettheten til luften i atmosfæren gitt en høyde
# Målenhet: Kilogram per kubikkmeter (kg/m³)
def density(t, v):
    h = height(t, v)
    return (pressure(h)/temperature(h)) * 3.4855

# Returnerer luftmostanden på raketten
def Fd(t, a, v):
    return (1/2)*(1/2)*density(t, v)*area(t)*(speed(t, a, v)**2)

def main():
    script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in

    # Henter bildene til figuren
    img_path_rocket_down = os.path.join(script_dir, "../Oppgave4og5/rocket_down.png")
    img_path_rocket_up_no_fuel = os.path.join(script_dir, "../Oppgave4og5/rocket_up_no_fuel.png")#'../Oppgave4og5/rocket_up_no_fuel.png'
    img_path_rocket_up = os.path.join(script_dir, "../Oppgave4og5/rocket_up.png")#'../Oppgave4og5/rocket_up.png'
    img_path_explosion = os.path.join(script_dir, "../Oppgave4og5/explosion.png")#'../Oppgave4og5/explosion.png'
    img_path = img_path_rocket_up # Setter startbilde

    # Lister som skal brukes i kalkulasjonene
    F, a, m, v = [], [], [], []

    # Initialiserer startverdier
    h = 0.001 # Initialiseres med en verdi større enn null for å komme inn i løkken. Verdien påvirker ingen av kalkulasjonene.
    m.append(mass(0)) # Massen til raketten umiddelbart under oppskytning
    v.append(0) # Raketten har enda ikke noe fart ved t = 0
    # Regner ut summen av kreftene på raketten umiddelbart under oppskytning. Ingen luftmostanden siden v = 0
    F.append(Fs(0, 0) - gravity_constant * ((earth_mass * m[0])/(earth_radius**2)))
    a.append(F[0] / m[0]) # Utregning av rakettens akselerasjon (F = ma)

    t = 1 # Tidstrinn i sekund
    down = False # Bolsk variabel for å sjekke om raketten har snudd

    while h >= 0: # Løkken kjører så lenge vi er over havnivå
        last_h = h # last_h blir bruk til å sjekke om raketten enda går oppover eller om den har snudd
        v.append(speed(t, a, v)) # Utregning av nåværende tidstrinnets fart

        # Sjekker om farten er positiv eller negativ
        # Er den positiv er raketten på vei opp. Er den negativ er raketten på vei ned.
        # Hvis raketten er på vei ned vil luftmotstanden virke på raketten med retning oppover, som er valgt til å være positiv retning
        if v[t] > 0 and not down:
            F.append(Fs(t, v) - (Fg(t, v) + Fd(t, a, v)))
        else:
            F.append(-Fg(t, v) + Fd(t, a, v))

        m.append(saturnV.massAddition()) # Henter massen til raketten i nåværende tidstrinn
        a.append(F[t] / m[t]) # Regner ut akselerasjonen i nåværende tidstrinn
        h = height(t, v) # Regner ut den nye høyden fra overflaten

        plt.clf() # Fjerner tidligere plot
        str = '\n'.join((
            r'Tid: %.2f$s$' % (t, ),
            r'Høyde: %.2f$m$' % (h, ),
            r'Fart: %.2f$\frac{m}{s}$' % (v[t], ),
            r'Akselerasjon: %.2f$\frac{m}{s²}$' % (a[t], )))
        props = dict(boxstyle='round', facecolor='blue', alpha=0.5)
        plt.text(-475, 23800000, str, fontsize=10, verticalalignment='top', bbox=props)
        plt.xlim(-500, 500)
        plt.xlabel('Jordens overflate')
        plt.ylim(0, 24000000)
        plt.ylabel('Høyde (m)')
        # Endrer på bildet til figuren basert på om raketten er tom for drivstoff eller om den er på vei nedover
        if t > time_stage_three and not down:
            img_path = img_path_rocket_up_no_fuel
        if last_h > h and not down:
            img_path = img_path_rocket_down
            down = True
        if h < 0:
            img_path = img_path_explosion
        img = mpimg.imread(img_path)
        plt.imshow(img, aspect = 'auto', extent = [-40, 40, h - 100000, h + 1000000])
        plt.pause(0.01)
        t = t + 1 # Videre til neste tidstrinn

    plt.show()

if __name__ == "__main__":
    main()
