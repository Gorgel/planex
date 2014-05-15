from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
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


#3
#a)
L_Sun = con.L_sun.cgs.value
M_Sun = con.M_sun.cgs.value
G = con.G.cgs.value
c = con.c.cgs.value
au = con.au.cgs.value
k_abs3 = 200
s_3 = 10**-4		#Assumption from lectures
rho_3 = 3*10**6/(10**6)	#Assumption from lectures
m_3 = rho_3*4/3 * np.pi*s_3**3

#F_rad = L_Sun*k_abs3 /(4*np.pi*c)
#F_grav = G*M_Sun
#beta = F_rad/F_grav

m = np.linspace(10**-2, 0.1, 2000)
F_rad = L_Sun*k_abs3 /(4*np.pi*c)
F_grav = G*M_Sun*m
beta = F_rad/F_grav

plt.figure(1)
plt.plot(m, beta)
plt.xlabel('mass [g]')
plt.ylabel(r'$\beta$')


#Large m -> small F_rad. High size with low mass gives efficient F_rad (grains absorb more)

#b)
rho = 1
a = 50*10**-4
rho_MgFeSiO = 3.7
m_grain = rho_MgFeSiO*(4/3)*np.pi*a**3

#Assume distance r = 1 au
dist = np.linspace(0.5*au, 6*au, 2000)
F_grav_b = G*M_Sun*m_grain/dist**2
F_rad_b = L_Sun*k_abs3 /(4*np.pi*dist**2*c)

beta_3 = F_rad_b/F_grav_b

v = ((F_grav_b + F_rad_b)/(np.pi*a**2*rho))**(1/2)	#sign convention?

plt.figure(2)
plt.plot(dist/au, v, color='red')
plt.ylabel('v')
plt.xlabel('distance in au')
plt.show()


