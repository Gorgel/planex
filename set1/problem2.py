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
m = np.array([])
rho_ext_frac = np.array([])
PextList = np.linspace(4*10**-10,10**-15,3000)
for P_ext in PextList:
    #P_ext = 8*10**-12		# -13 in Pa, -12 in cgs units (bary 'barr-ee')
    rho_ext = P_ext/(cs**2)		#7.7*10^-22 g cm^-3
    rho_c = 7.5*10**-20
    rho_ratio = rho_ext/rho_c	#0.1

    #e)



    psi_ext = -sp.log(P_ext/(rho_c*cs**2))


    def f(y, t):
        y0 = y[0]
        y1 = y[1]

        f0 = y1
        f1 = sp.exp(-y0) - (2/t)*y1
        return [f0, f1]

    y_init = [0, 0]
    t = sp.linspace(0.0001, 500, 5000)

    solu = sp.integrate.odeint(f, y_init, t)
    psi = solu[:,0]
    psi_diff = solu[:,1]

    rho_frac = sp.exp(-psi)

    xi_ext_index_1 = sp.where(rho_frac < rho_ratio)
    xi_ext_index = xi_ext_index_1[0][0]
    xi_ext = t[xi_ext_index]

    m_temp = np.sqrt(rho_ext/(4*np.pi*rho_c))*((xi_ext**2)*psi_diff[xi_ext_index])
    rho_ext_frac = np.append(rho_ext_frac, np.log10((rho_c/rho_ext)))
    m = np.append(m,m_temp)

plt.figure(1)
plt.plot(t, psi, color='blue', linewidth = 2)
plt.plot(t, rho_frac, color='red', linewidth = 2)

plt.figure(2)
plt.plot(rho_ext_frac, m, color='blue', linewidth = 2)
plt.show()
