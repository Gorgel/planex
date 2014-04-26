from __future__ import division
from scipy import integrate
from scipy.optimize import fsolve
from scipy import constants
from astropy import constants	#ap.constants.h.cgs
import scipy as sp
import matplotlib.pyplot as plt
import astropy as ap
from scipy.integrate import odeint
import numpy as np

mu = 2
T = 25
R = 8.31*10**7

cs = sp.sqrt(R*T/mu)		#32229 cm/s

#b)
m = np.array([])         #create an empty array for storing dimensionless mass values
rho_ext_frac = np.array([])         #create an empty array for storing values of the fraction rho_c/rho_ext
PextList = np.linspace(4*10**-10,10**-15,3000)         #make a list of external pressures

##########################
#
# this for-loop loopes over the list of external pressures creating a
# different endpoint on the bonner-evert sphere in each loop.
#
#########################
for P_ext in PextList:
    #P_ext = 8*10**-12		# -13 in Pa, -12 in cgs units (bary 'barr-ee')
    rho_ext = P_ext/(cs**2)		#7.7*10^-22 g cm^-3
    rho_c = 7.5*10**-20
    rho_ratio = rho_ext/rho_c	#0.1
    psi_ext = -sp.log(P_ext/(rho_c*cs**2))

    #defining a function that takes input argument y? explanation needed!
    def f(y, t):
        y0 = y[0]
        y1 = y[1]

        f0 = y1
        f1 = sp.exp(-y0) - (2/t)*y1
        return [f0, f1]



    y_init = [0, 0]         #create a list with initial conditions
    t = sp.linspace(0.0001, 500, 5000)         #create a list with equidistant integration "steps" (start,stop, nr of steps)

    #Python ODE solver. Takes the function f, inital conditions y_init and steps t as input
    solu = sp.integrate.odeint(f, y_init, t)
    psi = solu[:,0]          #saves the psi values in list psi
    psi_diff = solu[:,1]     #saves the derivative of psi in list psi_diff


    rho_frac = sp.exp(-psi)     #calculates the fraction rho/rho_c

    xi_ext_index_1ist = sp.where(rho_frac < rho_ratio)     #selects the indexes from the list of rho/rho_c values that are < rho_ratio
    xi_ext_index = xi_ext_index_1ist[0][0]         #select first element in list.
    xi_ext = t[xi_ext_index]         #uses the selected index to pick out the radius where the bonnor-evert sphere ends

    m_calc = np.sqrt(rho_ext/(4*np.pi*rho_c))*((xi_ext**2)*psi_diff[xi_ext_index])     #calculates the dimensionless mass
    rho_ext_frac = np.append(rho_ext_frac, np.log10((rho_c/rho_ext)))     #calculates rho_c/rho_ext and appends it to an array
    m = np.append(m,m_calc)     #append calculated dimensionless mass to array


# plot the results
plt.figure(1)
plt.plot(t, psi, color='blue', linewidth = 2)
plt.plot(t, rho_frac, color='red', linewidth = 2)

plt.figure(2)
plt.plot(rho_ext_frac, m, color='blue', linewidth = 2)
plt.show()
