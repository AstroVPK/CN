import math
import numpy as np
import pdb

# Value of alpha corresponds to Si
alpha_Si = 8.8e-5                                # thermal diffusivity (m^2/s)
# Value of alpha corresponds to Si3N4 with CNTs (Carbon Nano-Tubes)
alpha_Si3N4_CNTs = 9.142e-6                           # thermal diffusivity (m^2/s)
# Value of alpha corresponds to Si3N4 without CNTs (Carbon Nano-Tubes)
alpha_Si3N4 = 8.605e-6                           # thermal diffusivity (m^2/s)
# Value of alpha corresponds to 1% Carbon-Steel(Carbon Nano-Tubes)
alpha_1Steel = 1.172e-5                           # thermal diffusivity (m^2/s)
# Value of alpha corresponds to PVC
alpha_PVC = 8.0e-8                               # thermal diffusivity (m^2/s)
# Value of alpha corresponds to Pyrolytic graphite
alpha_PyroGraph = 1.22e-3                        # thermal diffusivity (m^2/s)
# Value of alpha corresponds to dry air
alpha_DryAir = 1.9e-5                            # thermal diffusivity (m^2/s)
# Value of alpha corresponds to wet air
alpha_WetAir = 2.0e-5                            # thermal diffusivity (m^2/s)
# Value of alpha corresponds to He
alpha_He = 1.9e-4                                # thermal diffusivity (m^2/s)
# Value of alpha corresponds to Silver
alpha_Ag = 1.6563e-4                             # thermal diffusivity (m^2/s)
# Value of alpha corresponds to Aluminium
alpha_Al = 9.7e-5                                # thermal diffusivity (m^2/s)
# Value of alpha corresponds to Aluminium 6061-T6 Alloy
alpha_Al_6061_T6_Alloy = 6.4e-5                  # thermal diffusivity (m^2/s)

delta_l = 1.0e-2                              # block size (m)
delta_t = 1.0e-1                              # timestep (s)
delta_l_sq = math.pow(delta_l, 2.0)

def dist(x, y, z, x0, y0, z0):
	return math.sqrt(math.pow((x - x0), 2.0) + math.pow((y - y0), 2.0) + math.pow((z - z0), 2.0))

def distSq(x, y, z, x0, y0, z0):
	return math.pow((x - x0), 2.0) + math.pow((y - y0), 2.0) + math.pow((z - z0), 2.0)

def r(x, y, z):                               # CFL number (<= 0.5 for accuracy)
	x0 = 1.0
	y0 = 1.0
	z0 = 0.0
	if x <= x0:
		alpha = alpha_Si3N4
	else:
		alpha = alpha_1Steel
	rval = (alpha*delta_t)/(2.0*delta_l_sq)
	if rval > 0.5:
		warnings.warn("r = %4.3e exceeds 0.5!")
	return rval

def T(x, y, z):                               # Initial temperature (K)
	x0 = 1.0
	y0 = 1.0
	z0 = 0.0
	l = dist(x, y, z, x0, y0, z0)
	temp = 10.0*math.exp(-dist(x, y, z, x0, y0, z0)/0.1)
	return temp

def resetT(x, y, z):
	x0 = 1.0
	y0 = 1.0
	z0 = 0.0
	l = dist(x, y, z, x0, y0, z0)
	if l < 0.5:
		temp = 100.0
	else:
		temp = None
	return temp