# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 20:34:58 2020

@author: Sietse de Vilder
"""

#coordinate system in the chord line!!!

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
    x1 = 0.174 #m
    x2 = 1.051
    x3 = 2.512
    xa = 0.3
    theta = 0.436332 #[rad], also, we assume constant theta, since initial
    #deflection angle is way larger than deflection
    #due twist
    #Doing Macaulay step function for the moment Mz(x):
    if x<x1:
        return doubleintegral
    
    #until actuator I:
    elif (x2-0.5*xa)>x>=x1:
        return -R1y*(x-x1) + doubleintegral
    
    #until hinge 2:
    elif x2>x>=(x2-0.5*xa):
        return -RIy*(x-(x2-0.5*xa)) - R1y*(x-x1) + doubleintegral
                
    #until actuator II:
    elif (x2+0.5*xa) > x >= x2:
        return -R2y*(x-x2) -RIy*(x-(x2-0.5*xa)) - R1y*(x-x1) + doubleintegral
                     
    #until hinge 3:
    elif x3>x>=(x2+0.5*xa):
        return  P*np.sin(theta)*(x-(x2+0.5*xa)) -R2y*(x-x2) -RIy*(x-(x2-0.5*xa)) - R1y*(x-x1) + doubleintegral
    
    #for the rest of the aileron
    elif x>x3:
        return -R3y*(x-x3)+ P*np.sin(theta)*(x-(x2+0.5*xa)) -R2y*(x-x2) -RIy*(x-(x2-0.5*xa)) - R1y*(x-x1) + doubleintegral
    
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
    x1 = 0.174 #m
    x2 = 1.051
    x3 = 2.512
    xa = 0.3
    theta = 0.436332 #[rad], also, we assume constant theta, since initial
    #deflection angle is way larger than deflection
    #due twist
    #Doing Macaulay step function for the shear in y-direction:
    if x<x1:
        return singleintegral 
    
    #until actuator I:
    elif (x2-0.5*xa)>x>=x1:
        return -R1y + singleintegral
    
    #until hinge 2:
    elif x2>x>=(x2-0.5*xa):
        return -RIz - R1z + singleintegral
                
    #until actuator II:
    elif (x2+0.5*xa) > x >= x2:
        return -R2z - RIz - R1z + singleintegral
                     
    #until hinge 3:
    elif x3>x>=(x2+0.5*xa):
        return  P*np.sin(theta) - R2z - RIz - R1z + singleintegral
    
    #for the rest of the aileron
    elif x>x3:
        return -R3z + P*np.sin(theta) - R2z - RIz - R1z + singleintegral

def Sigma_x(z,y,Mz_x,My_x):
    return ((Iyy*y - Izy*z)/(Izz*Iyy - Ixy**2))*Mz_x(R1y,RIy,R2y,P,R3y,Py,x) + ((Izz*z - Izy*y)/(Ixx*Iyy - Ixy**2))*My_x(R1z,RIz,R2z,P,R3z,x) #general equation for stress. Not sure about it.

    