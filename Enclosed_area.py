

#Input parameters
#DELETE THESE PARAMETERS AFTER VERIFICATION!!!!!



import numpy as np

#This code gives you the enclosed area of the aileron cross section
def Am(Ca,ha):
    
    Semicircular_section  = np.pi *0.5* (0.5*ha)**2
    
    Triangular_section = 0.5 * ha *(Ca-0.5*ha)
    
    #this array returns the total area, the semicircular area
    #and the triangular area in [m^2]
    
    return np.array([Semicircular_section+Triangular_section,Semicircular_section,Triangular_section])







