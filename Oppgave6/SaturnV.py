class SaturnV:
	# Init with all the necessary values
	def __init__(self,
				step_one,
				consumption_one,
				fuel_one,
				thrust_sea_one,
				thrust_vac_one,
				step_two,
				consumption_two,
				fuel_two,
				thrust_sea_two,
				thrust_vac_two,
				step_three,
				consumption_three,
				fuel_three,
				thrust_sea_three,
				thrust_vac_three):
		self.step_one = step_one
		self.consumption_one=consumption_one
		self.fuel_one=fuel_one
		self.thrust_sea_one=thrust_sea_one
		self.thrust_vac_one=thrust_vac_one
		self.step_two=step_two
		self.consumption_two=consumption_two
		self.fuel_two=fuel_two
		self.thrust_sea_two=thrust_sea_two
		self.thrust_vac_two=thrust_vac_two
		self.step_three=step_three
		self.consumption_three=consumption_three
		self.fuel_three=fuel_three
		self.thrust_sea_three=thrust_sea_three
		self.thrust_vac_three=thrust_vac_three

	def massAddition(self):
		if self.fuel_one > 0:
			return self.step_one + self.step_two + self.step_three + self.fuel_one + self.fuel_two + self.fuel_three
		elif self.fuel_two > 0:
			return self.step_two + self.step_three + self.fuel_two + self.fuel_three
		elif self.fuel_three > 0:
			return self.step_three + self.fuel_three
		else:
			return self.step_three

	#Calculates the mass of the rocket at a given time t
	def calculateMass(self,t):
		if 0 < self.fuel_one:
			self.fuel_one -= self.consumption_one
		elif 0 < self.fuel_two:
			self.fuel_two -= self.consumption_two
		elif 0 < self.fuel_three:
			self.fuel_three -= self.consumption_three
		return (self.massAddition())

	#Calculates the thrust of the rocket given atmospheric pressure and time (to get more accurate information)
	def calculateThrust(self, t, p):
		thrust = 0
		if self.fuel_one > 0:
			thrust = self.thrust_vac_one + (self.thrust_sea_one-self.thrust_vac_one)*p
		elif self.fuel_two > 0:
			thrust = self.thrust_vac_two + (self.thrust_sea_two-self.thrust_vac_two)*p
		elif self.fuel_three > 0:
			thrust = self.thrust_vac_three + (self.thrust_sea_three-self.thrust_vac_three)*p
		else:
			thrust = 0
		return thrust

def main():
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

if __name__ == "__main__":
    # execute only if run as a script
    main()
