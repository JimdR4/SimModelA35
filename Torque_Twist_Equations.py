# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 15:24:16 2020
@author: Sietse de Vilder
"""

#VERIFY THIS!!!!
#INCLUDE q(x,z)!!!!

#coordinate system in the chord line
#at the wing's leading edge

#Trial inputs:

P = 20.6 * 1000 #N

import numpy as np
from Aero_Load_Functions import aero_load_interpolation, aero_load_1_int_function, aero_load_2_int_function

def T_x(R1y,R2y,R3y,R_I,z_hat,x,P):
    
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
    
def Twist_x(G,J,R1y,z_hat,R2y,R3y,P,R_I,x,C5):
   
    coeff, q, z_dom = aero_load_interpolation('torque')
    val = aero_load_2_int_function(x, coeff, 'x')
    
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
        return (1/(G*J))*val + C5
    
    elif (x2-0.5*xa)>x>=x1:
        return (1/(G*J))*(-R1y*(x-x1)*(-0.5* h_a-z_hat) + val) + C5 
    
    elif x2>x>=(x2-0.5*xa):
        return (1/(G*J))*(-R1y*(x-x1)*(-0.5* h_a-z_hat) - R_I*np.sin(theta)*(x-x2+xa/2)*(-z_hat) + R_I*np.cos(theta)*(x-x2+xa/2)*0.5*h_a + val) + C5
    
    elif (x2+0.5*xa)>x>=x2:
        return (1/(G*J))*(-R1y*(x-x1)*(-0.5* h_a-z_hat) - R_I*np.sin(theta)*(x-x2+xa/2)*(-z_hat) + R_I*np.cos(theta)*(x-x2+xa/2)*0.5*h_a - R2y*(x-x2)*(-0.5*h_a-z_hat) + val)  + C5
    
    elif x3>x>=(x2+0.5*xa):
        return (1/(G*J))*(-R1y*(x-x1)*(-0.5* h_a-z_hat) - R_I*np.sin(theta)*(x-x2+xa/2)*(-z_hat) + R_I*np.cos(theta)*(x-x2+xa/2)*0.5*h_a - R2y*(x-x2)*(-0.5*h_a-z_hat) \
                + P*np.sin(theta)*(x-x2-xa/2)*(-z_hat) - P*np.cos(theta)*(x-x2-xa/2)*(0.5*(h_a)) + val) + C5
    
    elif x >= x3 :
        return (1/(G*J))*(-R1y*(x-x1)*(-0.5* h_a-z_hat) - R_I*np.sin(theta)*(x-x2+xa/2)*(-z_hat) + R_I*np.cos(theta)*(x-x2+xa/2)*0.5*h_a - R2y*(x-x2)*(-0.5*h_a-z_hat) \
                + P*np.sin(theta)*(x-x2-xa/2)*(-z_hat) - P*np.cos(theta)*(x-x2-xa/2)*(0.5*(h_a)) -R3y(x-x3)*(-0.5*h_a-z_hat) + val) + C5
