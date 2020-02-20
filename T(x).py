#VERIFY THIS!!!!
#INCLUDE q(x,z)!!!!

#coordinate system in the chord line
#at the wing's leading edge

#Trial inputs:

P = 20.6 * 1000 #N

import numpy as np
from Aero_Load_Functions import aero_load_interpolation, aero_load_1_int_function

def T_x(R1y,R2y,R3y,R_I,z_hat,x):
    
    coeff, q, x_dom = aero_load_interpolation('torque')
    val = aero_load_1_int_function(x, coeff, 'x')

    #Defining all the x locations:
    x1 = 0.174 #m
    x2 = 1.051
    x3 = 2.512
    xa = 0.3
    theta = 0.436332 #[rad], also, we assume constant theta, since initial
    #deflection angle is way larger than deflection
    #due twist
    h_a = 0.248 #[m], airfoil height
    #z_hat is the z location of the shear center
    #Doing Macaulay step function for the torque T(x):
    
    if x<x1:
        return val
    
    elif (x2-0.5*xa)>x>=x1:
        return -R1y * (-0.5* h_a-z_hat) + val 
    
    elif x2>x>=(x2-0.5*xa):
        return -R1y * (-0.5* h_a-z_hat) - R_I * np.sin(theta)* (-z_hat) + R_I * np.cos(theta) * 0.5 *h_a + val
    
    elif (x2+0.5*xa)>x>=x2:
        return -R1y * (-0.5* h_a-z_hat) - R_I * np.sin(theta)* (-z_hat) + R_I * np.cos(theta) * 0.5 *h_a - R2y * (-0.5*h_a-z_hat) + val
    
    elif x3>x>=(x2+0.5*xa):
        return -R1y * (-0.5* h_a-z_hat) - R_I * np.sin(theta)* (-z_hat) + R_I * np.cos(theta) * 0.5 *h_a - R2y * (-0.5*h_a-z_hat) \
                + P *np.sin(theta) * (-z_hat) - P * np.cos(theta) * (0.5*(h_a) ) + val
    
    elif x >= x3 :
        return -R1y * (-0.5* h_a-z_hat) - R_I * np.sin(theta)* (-z_hat) + R_I * np.cos(theta) * 0.5 *h_a - R2y * (-0.5*h_a-z_hat) \
                + P *np.sin(theta) * (-z_hat) - P * np.cos(theta) * (0.5*h_a) -R3y * (-0.5*h_a-z_hat) + val
    
