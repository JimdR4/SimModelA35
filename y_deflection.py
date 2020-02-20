# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 11:41:23 2020
@author: Sietse de Vilder
"""

#coordinate system in the chord line!!!
import numpy as np
from Aero_Load_Functions import aero_load_interpolation, aero_load_4_int_function

def V_x(R1y,RIy,R2y,P,R3y,theta,x,C1,C2,E,Izz):
    #This function Returns the y-deflection as a function of x
    
    coeff, q, x_dom = aero_load_interpolation('torque')
    val = aero_load_4_int_function(x, coeff, 'x')

    #Defining all the x locations:
    x1 = 0.174 #m
    x2 = 1.051
    x3 = 2.512
    xa = 0.3
    theta = 0.436332 #[rad], also, we assume constant theta, since initial
    #deflection angle is way larger than deflection
    #due twist
    #Doing Macaulay step function for the y-deflection along x:
    if x<x1:
        return ((-1/(E*Izz))*val) + C1*x + C2
    
    #until actuator I:
    elif (x2-0.5*xa)>x>=x1:
        return (-1/(E*Izz))*(-(R1y/6)*(x-x1)**3 + val) + C1*x + C2 
    
    #until hinge 2:
    elif x2>x>=(x2-0.5*xa):
        return (-1/(E*Izz))*(-(RIy/6)*(x-(x2-0.5*xa))**3 - (R1y/6)*(x-x1)**3 + val) + C1*x + C2 
                
    #until actuator II:
    elif (x2+0.5*xa) > x >= x2:
        return (-1/(E*Izz))*(-(R2y/6)*(x-x2)**3 - (RIy/6)*(x-(x2-0.5*xa))**3 - (R1y/6)*(x-x1)**3 + val) + C1*x + C2 
                     
    #until hinge 3:
    elif x3>x>=(x2+0.5*xa):
        return  (-1/(E*Izz))*((P/6)*np.cos(theta)*(x-(x2+0.5*xa))**3 - (R2y/6)*(x-x2)**3 - (RIy/6)*(x-(x2-0.5*xa))**3 - (R1y/6)*(x-x1)**3 + val) + C1*x + C2 
    
    #for the rest of the aileron
    elif x>x3:
        return (-1/(E*Izz))*(-(R3y/6)*(x-x3)**3 + (P/6)*np.cos(theta)*(x-(x2+0.5*xa))**3 -(R2y/6)*(x-x2)**3 -(RIy/6)*(x-(x2-0.5*xa))**3 - (R1y/6)*(x-x1)**3 + val) + C1*x + C2 
    