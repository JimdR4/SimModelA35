import numpy as np
from Aero_Load_Functions import aero_load_interpolation, aero_load_function, aero_load_1_int_function, aero_load_2_int_function, aero_load_3_int_function

def My_x(R1z,RIz,R2z,P,R3z,x,x1,x2,x3,xa,theta):
    #This function Returns My as a function of x
    
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

def Mz_x(R1y,RIy,R2y,P,R3y,Py,x,x1,x2,x3,xa,theta):
    #This function Returns Mz as a function of x
    
    #Defining all the x locations:
    coeff,q,x_dom = aero_load_interpolation('force')
    val = aero_load_2_int_function(x, coeff, 'x')
    
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

def Sz_x(R1z,RIz,R2z,P,R3z,x,x1,x2,x3,xa,theta):
    #This function Returns the shear in z-direction as a function of x
    
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

def Sy_x(R1y,RIy,R2y,P,R3y,x,x1,x2,x3,xa,theta):
    #This function Returns the shear in y-direction:
    
    #Defining all the x locations:
    coeff,q,x_dom = aero_load_interpolation('force')
    val = aero_load_1_int_function(x, coeff, 'x')
    
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

def T_x(R1y,R2y,R3y,R_I,zsc,x,P,x1,x2,x3,xa,theta,ha):
    
    coeff, q, x_dom = aero_load_interpolation('torque')
    val = aero_load_1_int_function(x, coeff, 'x')
    
    #Doing Macaulay step function for the torque T(x):
    
    if x<x1:
        return val
    
    elif (x2-0.5*xa)>x>=x1:
        return -R1y * (-0.5* ha-zsc) + val
    
    elif x2>x>=(x2-0.5*xa):
        return -R1y * (-0.5* ha-zsc) - R_I * np.sin(theta)* (-zsc) + R_I * np.cos(theta) * 0.5 *ha + val
    
    elif (x2+0.5*xa)>x>=x2:
        return -R1y * (-0.5* ha-zsc) - R_I * np.sin(theta)* (-zsc) + R_I * np.cos(theta) * 0.5 *ha - R2y * (-0.5*ha-zsc) + val
    
    elif x3>x>=(x2+0.5*xa):
        return -R1y * (-0.5* ha-zsc) - R_I * np.sin(theta)* (-zsc) + R_I * np.cos(theta) * 0.5 *ha - R2y * (-0.5*ha-zsc) \
                + P *np.sin(theta) * (-zsc) - P * np.cos(theta) * (0.5*(ha) ) + val
    
    elif x >= x3 :
        return -R1y * (-0.5* ha-zsc) - R_I * np.sin(theta)* (-zsc) + R_I * np.cos(theta) * 0.5 *ha - R2y * (-0.5*ha-zsc) \
                + P *np.sin(theta) * (-zsc) - P * np.cos(theta) * (0.5*ha) -R3y * (-0.5*ha-zsc) + val

