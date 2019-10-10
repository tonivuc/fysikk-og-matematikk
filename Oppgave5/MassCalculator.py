class MassCalculator:
	step_one = 2970000;
	one_consumption = 12857.1;
	burn_one = 168;
	step_two = 680000;
	two_consumption = 1266.9;
	burn_two = 360;
	step_three = 183000;
	three_consumption = 219;
	burn_three = 500;
	
	def __init__(self,
				step_one,
				one_consumption,
				burn_one,
				step_two,
				two_consumption,
				burn_two,
				step_three,
				three_consumption,
				burn_three):
		self.step_one = step_one;
		self.one_consumption=one_consumption;
		self.burn_one=burn_one;
		self.step_two=step_two;
		self.two_consumption=two_consumption;
		self.burn_two=burn_two;
		self.step_three=step_three;
		self.three_consumption=three_consumption;
		self.burn_three=burn_three;

	def calculateMass(self,t):
		mass = 0
		if t <= self.burn_one:
			mass = self.step_one - (self.one_consumption*t);
		elif t <= (self.burn_one + self.burn_two):
			mass = self.step_two - ((t-self.burn_one)*(self.two_consumption));
		elif t <= (self.burn_one + self.burn_two + self.burn_three):
			mass = self.step_three - ((t-(self.burn_one + self.burn_two))*(self.three_consumption));
		else:
			mass = self.step_three - ((self.burn_three)*(self.three_consumption));
		return mass

def main():
	m1 = 2970000;
	m2 = 680000;
	m3 = 183000;
	b1 = 168;
	b2 = 360;
	b3 = 500;
	c1 = 12857.1;
	c2 = 1266.9;
	c3 = 219;

	saturnV = MassCalculator(m1,c1,b1,m2,c2,b2,m3,c3,b3);
	mass_t = saturnV.calculateMass(15);
	print("The mass is: " + str(mass_t) + "kg")

if __name__ == "__main__":
    # execute only if run as a script
    main();