


#Verify this function!!!

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
        return  P*np.cos(theta) -R2z*(x-x2) -RIz*(x-(x2-0.5*xa)) - R1z*(x-x1)
    
    #for the rest of the aileron
    elif x>x3:
        return -R3z*(x-x3)+ P*np.cos(theta) -R2z*(x-x2) -RIz*(x-(x2-0.5*xa)) - R1z*(x-x1)
        