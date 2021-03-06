from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
#from astropy import constants as con

#exercise2_2.py
#a)

c = 29979245800		#con.c.cgs.value
h = 6.6260755*10**-27	#con.h.cgs.value

wavelength = 1*10**-4 		#1 micro m

#Comments:  	E = h * nu
#		wavelength = c/nu
#		p = m*v 
#		delta_p = h v_b /c + (2m h v_k)**(1/2) * cos(theta)
#		theta = 0 corresponds to particle being sent back into the direction from which the photon came.
#		hv_b = unbind energy, hv_k = kinetic energy
#		forward direction, theta = 180?
#		Photon energy , E = p*c, E = h*c/wavelength
#		=> p = h/wavelength

nu = c/wavelength
p_photon = h/wavelength
theta = np.pi
a = 10**-4			#assumption
rho = 3				#assumption
m = (4/3)*np.pi*a**3 * rho	#assumption
#delta_p = h*nu_b/c + (2*m*h*nu_k)**(1/2) * np.cos(theta)
#assume no binding energy freeing required?
delta_p = (2*m*h*nu)**(1/2) * np.cos(theta)

#b)
delta_p_back = (2*m*h*nu)**(1/2) * np.cos(0)

#c)
# P(theta) dtheta = (1-np.cos(theat))dtheta/2
# P = probability
# P = (1-cos(theat))/2
thetas = np.linspace(0, np.pi*2, 1000)
P = (1-np.cos(thetas))/(2*np.pi)
#P = (thetas - np.sin(thetas))/(2*np.pi) *(2/1000)

plt.figure(1)
plt.plot(thetas, P)
plt.show()
plt.xlabel(r'Angle $\theta$', fontsize=18)
plt.ylabel('Probability', fontsize=18)

#Easy to see that we get a net acceleration in the backward direction. (Scattering forward direction, however, momentum change negative resulting in an opposite direction acceleration).


