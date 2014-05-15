from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
#from astropy import constants as con

#exercise2_4.py

A_V = 2.90	#A_J
A_I = 1.83	#A_J
A_K = 0.37	#A_J
I_a = 10.53
I_bc = 11.23
J_a = 8.71
J_b = 10.05
J_c = 10.37
K_a = 6.53
K_b = 8.42
K_c = 8.82
K5_IJ = 0.80
K7_IJ = 0.92

###################################
###         Comments	##########
##################################

#E_ij = colour excess	Eij = (I- J)_obs - (I -J)_intrinsic
#A_ij = Extinction
#E_ij = A_i - A_j 
#A_J = A_I - E_IJ

#a)
E_IJ_a = (I_a - J_a) - K5_IJ
E_IJ_b = (I_bc - J_b) - K7_IJ
E_IJ_c = (I_bc - J_c) - K7_IJ
#A_V = 2.9*A_J = 2.9*(A_I - E_IJ) = 2.9*(1.83*A_J - E_IJ)
#E_IJ_a = A_I - A_J = 1.83A_J - AJ = 0.83 A_J
A_J_a = E_IJ_a/0.83
A_V_a = 2.9*A_J_a
I_b = A_V_a*0.83 + J_b+K7_IJ 
I_c = A_V_a*0.83 + J_c+K7_IJ

#b)
#E_JK_a = (J_a - K_a)_obs - (J_a - K_a)_intrinsic
#E_JK_a = A_J - A_K = A_J - 0.37*A_J = 0.63*A_J
#E_JK_a = (J_a - K_a)_obs - ()
E_JK_a = 0.63*A_J_a

#c)
BCJ_a = 1.41
BCJ_bc = 1.37
mu = 5.5		#distance modulus
M_Sun = 4.5

#################################33
#		Comments
############################3

#M_bol = M_J + BC(J)
#M = m - mu

M_bol_a = J_a - mu + BCJ_a
M_bol_b = J_b - mu + BCJ_bc
M_bol_c = J_c - mu + BCJ_bc
#M_bol - M_Sun = -2.5*log10(L_star / L_Sun)
#L_star = 10**(-(M_bol-M_Sun)/2.5) * L_sun

L_Sun = 3.846*10**33	#ap.con.L_sun.cgs.value

L_star_a = 10**(-(M_bol_a-M_Sun)/2.5) * L_Sun
L_star_b = 10**(-(M_bol_b-M_Sun)/2.5) * L_Sun
L_star_c = 10**(-(M_bol_c-M_Sun)/2.5) * L_Sun

#d)
T_a = 4340
T_bc = 4040

#a ~ 0.002 Gyear = 2 Myear 
#b ~  0.00636 Gyear = 6.36 Myear
#c ~  0.02005 Gyear = 20.05 Myear


