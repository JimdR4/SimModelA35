# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 20:34:58 2020

@author: Sietse de Vilder
"""

#coordinate system in the chord line!!!

from Aero_Load_Functions import aero_load_interpolation, aero_load_function, aero_load_1_int_function, aero_load_2_int_function, aero_load_3_int_function, simpsons_rule_integration


#SINGLE INTEGRAL
def aero_load_1_int_function(coordinates, coeff, direction):
    ''' input:  coordinates =   x values  
                coeff =         coefficients from the interpolation
                direction =     choose between 'x' or 'z'
        output: val =           function values for the coordinates used as input
    '''
    
    Nb = 41     # number of spanwise segments
    Nc = 81

    if direction == 'x':
        N = Nb
    if direction == 'z':
        N = Nc

    if np.size(coordinates) >1: # when there is more then 1 coordinate
        val = []
        coeff1 = np.zeros(len(coeff))
        for co in coordinates:
            vect = np.zeros(N)
            for i in range(N):
                vect[i]  = co**(i+1)
                coeff1[i] = coeff[i]/(i+1)
            sol = np.dot(vect,coeff1)
            val.append(sol)

    else: # in case there is only 1 coordinate
        vect = np.zeros(N)
        coeff1 = np.zeros(len(coeff))
        for i in range(N):
            vect[i]  = coordinates**(i+1)
            coeff1[i] = coeff[i]/(i+1)
        val = np.dot(vect,coeff1)

    return val

#DOUBLE INTEGRAL
def aero_load_2_int_function(coordinates, coeff, direction):
    ''' input:  coordinates =   x values  
                coeff =         coefficients from the interpolation
                direction =     choose between 'x' or 'z'
        output: val =           function values for the coordinates used as input
    '''
    
    Nb = 41     # number of spanwise segments
    Nc = 81

    if direction == 'x':
        N = Nb
    if direction == 'z':
        N = Nc

    if np.size(coordinates) >1: # when there is more then 1 coordinate
        val = []
        coeff1 = np.zeros(len(coeff))
        coeff2 = np.zeros(len(coeff))
        for co in coordinates:
            vect = np.zeros(N)
            for i in range(N):
                vect[i]   = co**(i+2)
                coeff1[i] = coeff[i]/(i+1)
                coeff2[i] = coeff1[i]/(i+2)
            sol = np.dot(vect,coeff2)
            val.append(sol)

    else: # in case there is only 1 coordinate
        vect = np.zeros(N)
        coeff1 = np.zeros(len(coeff))
        coeff2 = np.zeros(len(coeff))
        for i in range(N):
            vect[i]   = coordinates**(i+2)
            coeff1[i] = coeff[i]/(i+1)
            coeff2[i] = coeff[i]/(i+2)
        val = np.dot(vect,coeff2)

    return val

import numpy as np
def My_x(R1z,RIz,R2z,P,R3z,x):
    #This function Returns My as a function of x
    
    #Defining all the x locations:
    x1 = 0.174 #m
    x2 = 1.051
    x3 = 2.512
    xa = 0.3
    theta = 0.436332 #[rad], also, we assume constant theta, since initial
    #deflection angle is way larger than deflection
    #due twist
    #Doing Macaulay step function for the moment My(x):
    if x<x1:
        return 0
    
    #until actuator I:
    elif (x2-0.5*xa)>x>=x1:
        return -R1z*(x-x1)
    
    #until hinge 2:
    elif x2>x>=(x2-0.5*xa):
        return -RIz*(x-(x2-0.5*xa)) - R1z*(x-x1)
                
    #until actuator II:
    elif (x2+0.5*xa) > x >= x2:
        return -R2z*(x-x2) -RIz*(x-(x2-0.5*xa)) - R1z*(x-x1)
                     
    #until hinge 3:
    elif x3>x>=(x2+0.5*xa):
        return  P*np.cos(theta)*(x-(x2+0.5*xa)) -R2z*(x-x2) -RIz*(x-(x2-0.5*xa)) - R1z*(x-x1)
    
    #for the rest of the aileron
    elif x>x3:
        return -R3z*(x-x3)+ P*np.cos(theta)*(x-(x2+0.5*xa)) -R2z*(x-x2) -RIz*(x-(x2-0.5*xa)) - R1z*(x-x1)

import numpy as np
def Mz_x(R1y,RIy,R2y,P,R3y,Py,x):
    #This function Returns Mz as a function of x
    
    #Defining all the x locations:
    coeff,q,x_dom = aero_load_interpolation('force')
    val = aero_load_2_int_function(x, coeff, 'x')
    x1 = 0.174 #m
    x2 = 1.051
    x3 = 2.512
    xa = 0.3
    theta = 0.436332 #[rad], also, we assume constant theta, since initial
    #deflection angle is way larger than deflection
    #due twist
    #Doing Macaulay step function for the moment Mz(x):
    if x<x1:
        return val
    
    #until actuator I:
    elif (x2-0.5*xa)>x>=x1:
        return -R1y*(x-x1) + val 
    
    #until hinge 2:
    elif x2>x>=(x2-0.5*xa):
        return -RIy*(x-(x2-0.5*xa)) - R1y*(x-x1) + val
    
    #until actuator II:
    elif (x2+0.5*xa) > x >= x2:
        return -R2y*(x-x2) -RIy*(x-(x2-0.5*xa)) - R1y*(x-x1) + val
                     
    #until hinge 3:
    elif x3>x>=(x2+0.5*xa):
        return  P*np.sin(theta)*(x-(x2+0.5*xa)) -R2y*(x-x2) -RIy*(x-(x2-0.5*xa)) - R1y*(x-x1) + val
    
    #for the rest of the aileron
    elif x>x3:
        return -R3y*(x-x3)+ P*np.sin(theta)*(x-(x2+0.5*xa)) -R2y*(x-x2) -RIy*(x-(x2-0.5*xa)) - R1y*(x-x1) + val
    
import numpy as np
def Sz_x(R1z,RIz,R2z,P,R3z,x):
    #This function Returns the shear in z-direction as a function of x
    
    #Defining all the x locations:
    x1 = 0.174 #m
    x2 = 1.051
    x3 = 2.512
    xa = 0.3
    theta = 0.436332 #[rad], also, we assume constant theta, since initial
    #deflection angle is way larger than deflection
    #due twist
    #Doing Macaulay step function for the Shear in z-direction:
    if x<x1:
        return 0
    
    #until actuator I:
    elif (x2-0.5*xa)>x>=x1:
        return -R1z
    
    #until hinge 2:
    elif x2>x>=(x2-0.5*xa):
        return -RIz - R1z
                
    #until actuator II:
    elif (x2+0.5*xa) > x >= x2:
        return -R2z -RIz - R1z
                     
    #until hinge 3:
    elif x3>x>=(x2+0.5*xa):
        return  P*np.cos(theta) - R2z - RIz - R1z
    
    #for the rest of the aileron
    elif x>x3:
        return -R3z + P*np.cos(theta) - R2z - RIz - R1z

import numpy as np
def Sy_x(R1y,RIy,R2y,P,R3y,x):
    #This function Returns the shear in y-direction:
    
    #Defining all the x locations:
    coeff,q,x_dom = aero_load_interpolation('force')
    val = aero_load_1_int_function(x, coeff, 'x')
    x1 = 0.174 #m
    x2 = 1.051
    x3 = 2.512
    xa = 0.3
    theta = 0.436332 #[rad], also, we assume constant theta, since initial
    #deflection angle is way larger than deflection
    #due twist
    #Doing Macaulay step function for the shear in y-direction:
    if x<x1:
        return val
    
    #until actuator I:
    elif (x2-0.5*xa)>x>=x1:
        return -R1y + val
    
    #until hinge 2:
    elif x2>x>=(x2-0.5*xa):
        return -RIy - R1y + val
                
    #until actuator II:
    elif (x2+0.5*xa) > x >= x2:
        return -R2y - RIy - R1y + val
                     
    #until hinge 3:
    elif x3>x>=(x2+0.5*xa):
        return  P*np.sin(theta) - R2y - RIy - R1y + val
    
    #for the rest of the aileron
    elif x>x3:
        return -R3y + P*np.sin(theta) - R2y - RIy - R1y + val

def Sigma_x(z,y,Mz_x,My_x,Iyy,Izy,Izz):
    return ((Iyy*y - Izy*z)/(Izz*Iyy - Izy**2))*Mz_x(R1y,RIy,R2y,P,R3y,Py,x) + ((Izz*z - Izy*y)/(Izz*Iyy - Izy**2))*My_x(R1z,RIz,R2z,P,R3z,x) #general equation for stress. Not sure about it.

    