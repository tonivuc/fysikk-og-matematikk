class SaturnV:
	#Variables used with some example values
	step_one = 2970000;
	consumption_one = 12857.1;
	burn_one = 168;
	thrust_one = 35100000;
	step_two = 680000;
	consumption_two = 1266.9;
	burn_two = 360;
	thrust_two = 5141000;
	step_three = 183000;
	consumption_three = 219;
	burn_three = 500;
	thrust_three = 1033100;
	
	# Init with all the necessary values
	def __init__(self,
				step_one,
				consumption_one,
				burn_one,
				thrust_one,
				step_two,
				consumption_two,
				burn_two,
				thrust_two,
				step_three,
				consumption_three,
				burn_three,
				thrust_three):
		self.step_one = step_one;
		self.consumption_one=consumption_one;
		self.burn_one=burn_one;
		self.thrust_one=thrust_one;
		self.step_two=step_two;
		self.consumption_two=consumption_two;
		self.burn_two=burn_two;
		self.thrust_two=thrust_two;
		self.step_three=step_three;
		self.consumption_three=consumption_three;
		self.burn_three=burn_three;
		self.thrust_three=thrust_three;

	#Calculates the mass of the rocket at a given time t
	def calculateMass(self,t):
		mass = 0;
		if t <= self.burn_one:
			mass = self.step_one - (self.consumption_one*t);
		elif t <= (self.burn_one + self.burn_two):
			mass = self.step_two - ((t-self.burn_one)*(self.consumption_two));
		elif t <= (self.burn_one + self.burn_two + self.burn_three):
			mass = self.step_three - ((t-(self.burn_one + self.burn_two))*(self.consumption_three));
		else:
			mass = self.step_three - ((self.burn_three)*(self.consumption_three));
		return mass

	#Calculates the thrust of the rocket using the mass found in calculateMass at a certain time t
	def calculateThrust(self, t):
		mass = self.calculateMass(t);
		thrust = 0;
		if t <= self.burn_one:
			thrust = self.thrust_one / mass;
		elif t <= (self.burn_one + self.burn_two):
			thrust = self.thrust_two / mass;
		elif t <= (self.burn_one + self.burn_two + self.burn_three):
			thrust = self.thrust_three / mass;
		else:
			thrust = 0;
		return thrust;

			
def main():
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
	
	#Time variable
	t = 15;

	#Using the variables above with printouts:
	mass_t = saturnV.calculateMass(t);
	thrust_t = saturnV.calculateThrust(t)
	print("The mass is: " + str(mass_t) + "kg")
	print("The thrust is: " + str(thrust_t) + "m/s")


if __name__ == "__main__":
    # execute only if run as a script
    main();