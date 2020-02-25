'''
Created on Tue Feb 25 2020
@author: Frans Vossen
'''

import numpy as np 
from CubSplineInterpolation import aero_nodes, interpolation 
from Integration import simp_val, simp_func

def aero_interpolation(force_or_torque, Ca, la, Nx, Nz, sc_z, x):
    ''' Input:  force_or_torque =   Choose between 'force' or 'torque
                Ca =                length of the chord
                la =                length of the span
                Nx =                number of aero nodes in x direction (span)
                Nz =                number of aero nodes in z direction (chord)
                sc_z =              z coordinate of shear center
                x =                 either a single value or array that will be used as function variable
        Output: qx_dom =    the function value(s) at the evaluated x

        Remark: changing the value 'n' will lead to a change in accuray of the numerical model. An increase in 'n' will lead to finer grid that is used to evaluate the cubic splines. 
                This results in a more accurate representation of the created cubic splines. This finer grid used in the cubic splines is also used for the numerical integration, 
                where three grid points are used to evaluate the area, which results in one grid point in the integrated function. This integrated function will therefore have less grid points then 
                the cubic splines. To deal with this problem and to represent the integral as a continous function, the integral function will be interpolated using the created grid points from the simpson's rule.
    '''
    n = 1001 # Change this value to obtain high enough accuracy

    aero_data = '/Users/frans/Desktop/AE/AE3212_II-Simulation_Varification_and_Validation/Python/Local/aerodynamicloaddo228.dat'
    data = np.genfromtxt(aero_data,delimiter=',') # column is spanwise // row is chordwise

    ''' Start program '''
    # Constants
    z_nodes, x_nodes = aero_nodes(Ca, la, Nz, Nx)

    #-----
    qx_disc = np.zeros(len(x_nodes))
    for bi in range(Nx): # repeating process for each spanwise location
    
        z_dom = np.linspace(z_nodes[0],z_nodes[-1],n)                               # create a z-axis on which to evaluate the created splines
        
        if force_or_torque == 'force':
            qz_disc = data[:,bi]                                                    # array containing all discretized point loads along the chord
            qz_dom  = interpolation(z_dom, qz_disc, z_nodes)                        # array containing all point loads on the defined z_dom, using cubic splines
            qz_int  = simp_val(0, Ca, n, qz_dom)                                    # numerical integration over the chord (from 0 to Ca) which results in a single value

            qx_disc[bi] = qz_int

        if force_or_torque == 'torque':
            qzt_disc = data[:,bi]*z_nodes                                           # array containing all discretized point loads along the chord multiplied by its z coordinate
            qzt_dom  = interpolation(z_dom, qzt_disc, z_nodes)                      # array containing all point loads along the chord multiplied by its z coordinate on the defined z_dom, using cubic splines
            qzt_int  = simp_val(0, Ca, n, qzt_dom)                                  # numerical integration over the chord (from 0 to Ca) which results in a single value

            qz_disc = data[:,bi]                                                    # array containing all discretized point loads along the chord
            qz_dom  = interpolation(z_dom, qz_disc, z_nodes)                        # array containing all point loads on the defined z_dom, using cubic splines
            qz_int  = simp_val(0, Ca, n, qz_dom)                                    # numerical integration over the chord (from 0 to Ca) which results in a single value

            z_applic = qzt_int / qz_int
            t_z      = qz_int * z_applic

            qx_disc[bi] = t_z

    qx_dom = interpolation(x, qx_disc, x_nodes) * 10**3                             # going from kN... to N...

    return qx_dom

def aero_int_1(force_or_torque, Ca, la, Nx, Nz, sc_z, x):
    ''' Input:  force_or_torque =   Choose between 'force' or 'torque
                Ca =                length of the chord
                la =                length of the span
                Nx =                number of aero nodes in x direction (span)
                Nz =                number of aero nodes in z direction (chord)
                sc_z =              z coordinate of shear center
                x =                 either a single value or array that will be used as function variable
        Output: int1 =  the function value(s) at the evaluated x
    '''
    n = 1001

    x_nodes = aero_nodes(Ca, la, Nz, Nx)[1]
    x_dom   = np.linspace(x_nodes[0], x_nodes[-1], n)                               # create an x-axis on which to evaluate the created splines
                                 
    values = aero_interpolation(force_or_torque, Ca, la, Nx, Nz, sc_z, x_dom)       # get the aerodynamic load qx on the created x_dom
    int1    = simp_func(x, n, x_dom, values)                                         # get the integral of the aerodynamic load int(qx) value(s) 

    return int1

def aero_int_2(force_or_torque, Ca, la, Nx, Nz, sc_z, x):
    ''' Input:  force_or_torque =   Choose between 'force' or 'torque
                Ca =                length of the chord
                la =                length of the span
                Nx =                number of aero nodes in x direction (span)
                Nz =                number of aero nodes in z direction (chord)
                sc_z =              z coordinate of shear center
                x =                 either a single value or array that will be used as function variable
        Output: int2 =  the function value(s) at the evaluated x
    '''
    n = 1001

    x_nodes = aero_nodes(Ca, la, Nz, Nx)[1]
    x_dom   = np.linspace(x_nodes[0], x_nodes[-1], n)                               # create an x-axis on which to evaluate the created splines
                                 
    int1_dom = aero_int_1(force_or_torque, Ca, la, Nx, Nz, sc_z, x_dom)             # get the integral of the aerodynamic load values
    int2    = simp_func(x, n, x_dom, int1_dom)                                      # get the second integral of the aerodynamic load value(s) 

    return int2

def aero_int_3(force_or_torque, Ca, la, Nx, Nz, sc_z, x):
    ''' Input:  force_or_torque =   Choose between 'force' or 'torque
                Ca =                length of the chord
                la =                length of the span
                Nx =                number of aero nodes in x direction (span)
                Nz =                number of aero nodes in z direction (chord)
                sc_z =              z coordinate of shear center
                x =                 either a single value or array that will be used as function variable
        Output: int3 =  the function value(s) at the evaluated x
    '''
    n = 1001

    x_nodes = aero_nodes(Ca, la, Nz, Nx)[1]
    x_dom   = np.linspace(x_nodes[0], x_nodes[-1], n)                               # create an x-axis on which to evaluate the created splines
                                 
    int2_dom = aero_int_2(force_or_torque, Ca, la, Nx, Nz, sc_z, x_dom)             # get the secon integral of the aerodynamic load values 
    int3     = simp_func(x, n, x_dom, int2_dom)                                     # get the third integral of the aerodynamic load value(s)

    return int3

def aero_int_4(force_or_torque, Ca, la, Nx, Nz, sc_z, x):
    ''' Input:  force_or_torque =   Choose between 'force' or 'torque
                Ca =                length of the chord
                la =                length of the span
                Nx =                number of aero nodes in x direction (span)
                Nz =                number of aero nodes in z direction (chord)
                sc_z =              z coordinate of shear center
                x =                 either a single value or array that will be used as function variable
        Output: int4 =  the function value(s) at the evaluated x
    '''
    n = 1001

    x_nodes = aero_nodes(Ca, la, Nz, Nx)[1]
    x_dom   = np.linspace(x_nodes[0], x_nodes[-1], n)                               # create an x-axis on which to evaluate the created splines
                                 
    int3_dom = aero_int_3(force_or_torque, Ca, la, Nx, Nz, sc_z, x_dom)             # get the third integral of the aerodynamic load values
    int4     = simp_func(x, n, x_dom, int3_dom)                                     # get the fourth integral of the aerodynamic load value(s)

    return int4
