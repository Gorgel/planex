from __future__ import division
import numpy as np
import matplotlib as plt
from astropy import constants as con


# Exercise set 2

#1. 
#a)

d = 50*10**-3	#grain radius
k_abs = 100	
k_em = np.pi*d**2 		#emissivity cross-section
T_d = 1400
T_Sun = 5777
T_W = 150000
R_Sun = con.R_sun.cgs.value
R_W = 0.01*R_Sun

r_Sun = np.sqrt(k_abs/k_em)*R_Sun*T_Sun**2 /(2*T_d**2)
r_au = r_Sun/con.au.cgs.value
Sun_radii = r_Sun/R_Sun

r_W = np.sqrt(k_abs/k_em)*R_W*T_W**2 /(2*T_d**2)
r_Wau = r_W/con.au.cgs.value
W_radii = r_W/R_W

#Dependent of opacities? Yes, very much e.g. k_em = 0.34 yields r_au = 0.69
#A varying cross-section would be good for heating at optical wavelengths but not as efficient for cooling at higher wavelengths e.g. IR. Cooling would be more efficient at higher wavelengths. See figures from lectures (HA-HA!) product of absorption/emission with opacity gives higher value for cooling at higher wavelengths.

#b) 
#Obvious reasons say that the first particle should be hotter because of higher absorption in the optical. The other particle does not absorb anything in the optical and therefor does not even heat up. Nothing is known of the cooling and the absorption cross-section for higher wavelengths is the same for both particles.

#c)

