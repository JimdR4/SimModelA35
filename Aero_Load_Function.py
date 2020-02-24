'''
Created on Mon Feb 24 2020
@author: Frans Vossen
'''

import numpy as np 
from CubSplineInterpolation import aero_nodes, interpolation, funct_cub_spline
from Integration import val_simpsons_rule_integration, funct_simpsons_rule_integration
from matplotlib import pyplot as plt # for testing

def aero_load_interpolation(force_or_torque, n, Ca, la, Nb, Nc, sc_z):
    ''' Input:  force_or_torque =   Choose between 'force' or 'torque
                n =                 The number of nodes used when setting an axis at which the cubic splines will be evaluated
        Output: coeff =             These are the coefficients obtained from the interpolation
    '''
    
    aero_data = '/Users/frans/Desktop/AE/AE3212_II-Simulation_Varification_and_Validation/Python/Local/aerodynamicloaddo228.dat'
    data = np.genfromtxt(aero_data,delimiter=',') # column is spanwise // row is chordwise

    ''' Start program '''
    # Constants
    z_nodes, x_nodes = aero_nodes(Ca, la, Nc, Nb)

    #-----
    int_along_chord = np.zeros(len(x_nodes))
    for bi in range(Nb): # repeating process for each spanwise location
        # Evaluating along the chord & iterating along the span
        q_cs = data[:,bi]               # array containing all point loads along the chord
        q_cs_times_x = q_cs*z_nodes     # array containing all point loads multiplied with its z coordinate (for torque calculation)

        if force_or_torque == 'force':
            z_dom = np.linspace(z_nodes[0],z_nodes[-1],n)                           # create a z-axis on which to evaluate the created splines

            coeff_along_chord = interpolation(q_cs,z_nodes)                         # get the coefficients of the splines that interpolate the point loads
            f_along_chord = funct_cub_spline(z_dom, z_nodes, coeff_along_chord)     # function values obtained from the splines
            q_c = val_simpsons_rule_integration(z_dom, n, f_along_chord)            # integrate function values over the chord
            int_along_chord[bi] = q_c

        if force_or_torque == 'torque':
            z_dom = np.linspace(z_nodes[0], z_nodes[-1], n)                         # create a z-axis on which to evaluate the created splines

            coeff_q_cs_x = interpolation(q_cs_times_x, z_nodes)                     # find the coefficients of the splines for the q values multiplied by its z value
            f_q_cs_x = funct_cub_spline(z_dom, z_nodes, coeff_q_cs_x)               # get the function values using the created spline
            int_q_cs_x = val_simpsons_rule_integration(z_dom, n, f_q_cs_x)          # integrate these function values along the chord

            coeff_along_chord = interpolation(q_cs, z_nodes)                        # get the coefficients of the splines that interpolate the point loads
            f_along_chord = funct_cub_spline(z_dom, z_nodes, coeff_along_chord)     # function values obtained from the splines
            q_c = val_simpsons_rule_integration(z_dom, n, f_along_chord)            # integrate function values over the chord

            z_point_of_applic = int_q_cs_x/q_c                                      # find the z location where the integrated point load acts on
            int_along_chord[bi] = q_c*(z_point_of_applic - sc_z)                    # multiply this point load by the arm to the shear center to get the torque due to the aerodynamic load

    coeff = interpolation(int_along_chord, x_nodes)                                 # get the coefficients for the splines that interpolates along the span

    return coeff
    
''' For testing '''

# Only change these values !!!
Ca = 0.515  # [m] / chord length 
la = 2.691  # [m] / span length
Nz = 81
Nx = 41
sc_z = -0.2
n  = 1001
force_or_torque = 'force' # either 'torque' or 'force'


#----------------
z_nodes, x_nodes = aero_nodes(Ca, la, Nz, Nx)

coeff = aero_load_interpolation(force_or_torque,n, Ca, la, Nx, Nz, sc_z)

x_dom = np.linspace(x_nodes[0],x_nodes[-1],n)
y   = funct_cub_spline(x_dom, x_nodes, coeff)

int1, x_dom_int1 = funct_simpsons_rule_integration(x_dom, n, y)
int2, x_dom_int2 = funct_simpsons_rule_integration(x_dom_int1, n, int1)
int3, x_dom_int3 = funct_simpsons_rule_integration(x_dom_int2, n, int2)
int4, x_dom_int4 = funct_simpsons_rule_integration(x_dom_int3, n, int3)


plt.plot(x_dom, y)
plt.plot(x_dom_int1, int1)
plt.plot(x_dom_int2, int2)
plt.plot(x_dom_int3, int3)
plt.plot(x_dom_int4, int4)
plt.show()
