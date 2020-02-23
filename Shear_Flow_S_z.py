'''
Created on Fri Feb 21 2020
@author: Frans Vossen
'''

import numpy as np 
from matplotlib import pyplot as plt

from Z_Coordinates_Stringers import z_coordinates_stringers
from MOI_and_Centroid import MOI_and_centroid

def q_b12(S_z, I_yy, n_stringer, Ca, h, t_sk, z_ce, B, omega1):
    ''' Input:      S_z = shear force in z direction
                    I_yy = moment of inertia around the y axis
                    n_stringer = the number of stringers used
                    Ca = lenth of chord
                    h = height of aileron
                    t_sk = skin thickness
                    z_sc = z coordinate of shear center
                    B = boom area
                    omega1 = defined as clockwise positive // works for both single value and array of values
        Output:     values = shear flow at omega1 // singel value or array of values
                    if wanted remove #: z_glob_co = the z-coordinate in the global coordinate frame corresponding to omega 1 // works for both single value of omega1 and array of omega1
                    if wanted remove #: y_glob_co = the y-coordinate in the global coordinate frame corresponding to omega 1 // works for both single value of omega1 and array of omega1
        Remarks:    - This is the shear flow in the circular part starting at the leading edge going clockwise
                    - omega 1 runs from 0 to pi/2
                    - Global coordinate system is defined at leading edge with y-axis pointing upwards and z-axis pointing away from the leading edge
                    - The coordinate system used for the integrations has the same orientation as the global coordinate system, but has its origin the the centroid
                    - Z-axis is assumed to be axis of symmetry
    '''
    if np.size(omega1) > 1:
        values = []
        # z_glob_co = []
        # y_glob_co = []
        for omega in omega1:
            if omega < 0 or omega > np.pi/2:
                raise ValueError('Omega1 should have a value between 0 and pi/2')

            R = h/2
            z = -(1-np.cos(omega))*R - z_ce
            z_co, str_loc = z_coordinates_stringers(Ca, h, n_stringer)[:2]
            val = -S_z/I_yy * (-t_sk*R**2*omega - t_sk*R**2*np.sin(omega) + t_sk*z_ce*R*omega + 0.5*B*z_co[0])

            z_glob_c = z + z_ce
            # y_glob_c = R*np.sin(omega)

            for i in range(len(z_co)):
                if str_loc[i] == 'c':
                    if abs(z_glob_c) >= abs(z_co[i]): 
                        val = val + -S_z/I_yy* B* z_co[i]
            
            values.append(val)
            # z_glob_co.append(z_glob_c)
            # y_glob_co.append(y_glob_c)

    else:
        if omega1 < 0 or omega1 > np.pi/2:
            raise ValueError('Omega1 should have a value between 0 and pi/2')

        R = h/2
        z = -(1-np.cos(omega1))*R - z_ce
        z_co, str_loc = z_coordinates_stringers(Ca, h, n_stringer)[:2]
        values = -S_z/I_yy * (-t_sk*R**2*omega1 - t_sk*R**2*np.sin(omega1) + t_sk*z_ce*R*omega1 + 0.5*B*z_co[0])

        z_glob_co = z + z_ce
        # y_glob_co = R*np.sin(omega1)

        for i in range(len(z_co)):
            if str_loc[i] == 'c':
                if abs(z_glob_co) >= abs(z_co[i]): 
                    values = values + -S_z/I_yy* B* z_co[i]
    

    return values # , z_glob_co, y_glob_co

def q_b42(S_z, I_yy, h, t_sp, z_ce, y1):
    ''' Input:      S_z = shear force in z direction
                    I_yy = moment of inertia around the y axis
                    h = height of aileron
                    t_sp = spar thickness
                    z_sc = z coordinate of shear center
                    y1 = defined as upwards positive // works for both single value and array of values
        Output:     values = shear flow at y // singel value or array of values
                    if wanted remove #: z_glob_co = the z-coordinate in the global coordinate frame corresponding to y1 // works for both single value of y1 and array of y1
                    if wanted remove #: y_glob_co = the y-coordinate in the global coordinate frame corresponding to y1 // works for both single value of y1 and array of y1
        Remarks:    - This is the shear flow in the spar starting at the chord going upwards
                    - y1 runs from 0 to R
                    - Global coordinate system is defined at leading edge with y-axis pointing upwards and z-axis pointing away from the leading edge
                    - The coordinate system used for the integrations has the same orientation as the global coordinate system, but has its origin the the centroid
                    - Z-axis is assumed to be axis of symmetry
    '''
    if np.size(y1) > 1:
        values = []
        # z_glob_co = []
        # y_glob_co = []
        for y in y1:
            R = h/2
            if y < 0 or y > R:
                raise ValueError('y should have a vale between 0 and R')

            val = -S_z/I_yy * (-R*t_sp*y - z_ce*t_sp*y)

            # z_glob_c = -R
            # y_glob_c = y

            values.append(val)
            # z_glob_co.append(z_glob_c)
            # y_glob_co.append(y_glob_c)

    else:
        R = h/2
        if y1 < 0 or y1 > R:
            raise ValueError('y should have a vale between 0 and R')

        values = -S_z/I_yy * (-R*t_sp*y1 - z_ce*t_sp*y1)

        # z_glob_co = -R
        # y_glob_co = y1

    return values # , z_glob_co, y_glob_co

def q_b23(S_z, I_yy, n_stringer, Ca, h, t_sk, z_ce, B, s1):
    ''' Input:      S_z = shear force in z direction
                    I_yy = moment of inertia around the y axis
                    n_stringer = the number of stringers used
                    Ca = lenth of chord
                    h = height of aileron
                    t_sk = skin thickness
                    z_sc = z coordinate of shear center
                    B = boom area
                    s1 = defined as clockwise positive // works for both single value and array of values
        Output:     values = shear flow at s1 // singel value or array of values
                    if wanted remove #: z_glob_co = the z-coordinate in the global coordinate frame corresponding to s1 // works for both single value of s1 and array of s1
                    if wanted remove #: y_glob_co = the y-coordinate in the global coordinate frame corresponding to s1 // works for both single value of s1 and array of s1
        Remarks:    - This is the shear flow in the straight part of the skin starting on top of the aileron going downwards to the right
                    - s1 runs from 0 to the length of the straight part, l_sk
                    - Global coordinate system is defined at leading edge with y-axis pointing upwards and z-axis pointing away from the leading edge
                    - The coordinate system used for the integrations has the same orientation as the global coordinate system, but has its origin the the centroid
                    - Z-axis is assumed to be axis of symmetry
    '''
    if np.size(s1) > 1:
        values = []
        # z_glob_co = []
        # y_glob_co = []
        for s in s1:
            R = h/2
            z_co, str_loc, l_sk = z_coordinates_stringers(Ca, h, n_stringer)

            if s < 0 or s > l_sk:
                raise ValueError('s1 should have a value between 0 and length of the straight skin')

            z = -R-z_ce - (Ca-R)/l_sk * s
            val = -S_z/I_yy * (-t_sk*R*s - t_sk*z_ce*s - t_sk*(Ca-R)/(2*l_sk)*s**2) + q_b12(S_z, I_yy, n_stringer, Ca, h, t_sk, z_ce, B, np.pi/2) + q_b42(S_z, I_yy, h, t_sp, z_ce, R)#[0] (also behind q_b12()), when evaluating z_global_co and y_glob_co[0]

            z_glob_c = z+z_ce
            # y_glob_c = R - R/l_sk * s

            for i in range(len(z_co)):
                if str_loc[i] == 's':
                    if abs(z_glob_c) >= abs(z_co[i]): 
                        val = val + -S_z/I_yy*B*z_co[i]
            
            values.append(val)
            # z_glob_co.append(z_glob_c)
            # y_glob_co.append(y_glob_c)

    else:
        R = h/2
        z_co, str_loc, l_sk = z_coordinates_stringers(Ca, h, n_stringer)

        if s1 < 0 or s1 > l_sk:
            raise ValueError('s1 should have a value between 0 and length of the straight skin')

        z = -R-z_ce - (Ca-R)/l_sk * s1
        values = -S_z/I_yy * (-t_sk*R*s1 - t_sk*z_ce*s1 - t_sk*(Ca-R)/(2*l_sk)*s1**2) + q_b12(S_z, I_yy, n_stringer, Ca, h, t_sk, z_ce, B, np.pi/2)[0] + q_b42(S_z, I_yy, h, t_sp, z_ce, R)#[0], when evaluating z_global_co and y_glob_co

        z_glob_co = z+z_ce
        # y_glob_co = R - R/l_sk * s1

        for i in range(len(z_co)):
            if str_loc[i] == 's':
                if abs(z_glob_co) >= abs(z_co[i]): 
                    values = values + -S_z/I_yy*B*z_co[i]
            else:
                continue

    return values # , z_glob_co, y_glob_co

def q_b35(S_z, I_yy, n_stringer, Ca, h, t_sk, z_ce, B, s2):
    ''' Input:      S_z = shear force in z direction
                    I_yy = moment of inertia around the y axis
                    n_stringer = the number of stringers used
                    Ca = lenth of chord
                    h = height of aileron
                    t_sk = skin thickness
                    z_sc = z coordinate of shear center
                    B = boom area
                    s1 = defined as clockwise positive // works for both single value and array of values
        Output:     val = shear flow at s2 // singel value or array of values
                    if wanted remove #: z_glob_co = the z-coordinate in the global coordinate frame corresponding to s2 // works for both single value of s2 and array of s2
                    if wanted remove #: y_glob_co = the y-coordinate in the global coordinate frame corresponding to s2 // works for both single value of s2 and array of s2
        Remarks:    - This is the shear flow in the straight part of the skin starting at the trailing edge going downwards to the left
                    - s2 runs from 0 to the length of the straight part, l_sk
                    - Global coordinate system is defined at leading edge with y-axis pointing upwards and z-axis pointing away from the leading edge
                    - The coordinate system used for the integrations has the same orientation as the global coordinate system, but has its origin the the centroid
                    - Z-axis is assumed to be axis of symmetry
    '''
    l_sk = z_coordinates_stringers(Ca, h, n_stringer)[2]

    if np.size(s2) > 1:
        for s in s2:
            if s < 0 or s> l_sk:
                raise ValueError('s1 should have a value between 0 and length of the straight skin') 

    else: 
        if s2 < 0 or s2> l_sk:
            raise ValueError('s1 should have a value between 0 and length of the straight skin') 
    
    values = q_b23(S_z, I_yy, n_stringer, Ca, h, t_sk, z_ce, B, l_sk-s2)
    # values, z, y = q_b23(S_z, I_yy, n_stringer, Ca, h, t_sk, z_ce, B, l_sk-s2) // when evaluating z_global_co and y_glob_co
    values = np.asarray(values)
    values = values*-1

    # z_glob_co = z
    # y_glob_co = y*-1

    return values # , z_glob_co, y_glob_co

def q_b54(S_z, I_yy, h, t_sp, z_ce, y2):
    ''' Input:      S_z = shear force in z direction
                    I_yy = moment of inertia around the y axis
                    h = height of aileron
                    t_sp = spar thickness
                    z_sc = z coordinate of shear center
                    y = defined as upwards positive // works for both single value and array of values
        Output:     val = shear flow at y // singel value or array of values
                    if wanted remove #: z_glob_co = the z-coordinate in the global coordinate frame corresponding to y2 // works for both single value of y2 and array of y2
                    if wanted remove #: y_glob_co = the y-coordinate in the global coordinate frame corresponding to y2 // works for both single value of y2 and array of y2
        Remarks:    - This is the shear flow in the spar starting at the bottom of the aileron going upwards
                    - y runs from -R to 0
                    - Global coordinate system is defined at leading edge with y-axis pointing upwards and z-axis pointing away from the leading edge
                    - The coordinate system used for the integrations has the same orientation as the global coordinate system, but has its origin the the centroid
                    - Z-axis is assumed to be axis of symmetry
    '''
    R = h/2

    if np.size(y2) > 1:
        for y in y2:
            if y > 0 or y < -R:
                raise ValueError('y should have a vale between -R and 0')
        
    else:
        if y2 > 0 or y2 < -R:
            raise ValueError('y should have a vale between -R and 0')

    values = q_b42(S_z, I_yy, h, t_sp, z_ce, -1*y2)
    # values, z, y = q_b42(S_z, I_yy, h, t_sp, z_ce, -1*y2) // when evaluating z_global_co and y_glob_co
    values = np.asarray(values)
    values = values*-1

    # z_glob_co = z
    # y_glob_co = y*-1 

    return values # , z_glob_co, y_glob_co

def q_b51(S_z, I_yy, n_stringer, Ca, h, t_sk, z_ce, B, omega2):
    ''' Input:      S_z = shear force in z direction
                    I_yy = moment of inertia around the y axis
                    n_stringer = the number of stringers used
                    Ca = lenth of chord
                    h = height of aileron
                    t_sk = skin thickness
                    z_sc = z coordinate of shear center
                    B = boom area
                    omega2 = defined as clockwise positive
        Output:     val = shear flow at omega1
        Remarks:    - This is the shear flow in the circular part starting at the bottom of the aileron going clokwise
                    - omega 1 runs from -pi/2 to 0
                    - Global coordinate system is defined at leading edge with y-axis pointing upwards and z-axis pointing away from the leading edge
                    - The coordinate system used for the integrations has the same orientation as the global coordinate system, but has its origin the the centroid
                    - Z-axis is assumed to be axis of symmetry
    '''
    if np.size(omega2) > 1:
        for omega in omega2:
            if omega > 0 or omega < -np.pi/2:
                raise ValueError('Omega1 should have a value between -pi/2 and 0')

    else:
        if omega2 > 0 or omega2 < -np.pi/2:
            raise ValueError('Omega1 should have a value between -pi/2 and 0')

    values = q_b12(S_z, I_yy, n_stringer, Ca, h, t_sk, z_ce, B, -1*omega2)
    # values, z, y = q_b12(S_z, I_yy, n_stringer, Ca, h, t_sk, z_ce, B, -1*omega2) // when evaluating z_global_co and y_glob_co
    values = np.asarray(values)
    values = values*-1

    # z_glob_co = z
    # y_glob_co = y*-1

    return values # , z_glob_co, y_glob_co

''''''''' FOR TESTING '''''''''
''' Remarks:    Coordinate system at leading edge with y-axis pointing upwards and z-axis pointing away from the aileron
                These shearflows are caused only due to a horizontal shear force S_z(x) which acts along the chord and therefore through the shear center '''

# Input values
S_z = 10
B   = 5 *10**(-4) # going from dm**2 to m**2
I_yy, I_zz, z_ce = MOI_and_centroid()   # correct values // z_ce is a negative value!!
t_sk = 1.1 * 10**(-3)                   # correct value
t_sp = 2.2 * 10**(-3)                   # correct value
h    = 0.248                            # correct value
R    = h/2                              # correct value
Ca   = 0.515                            # correct value
n_stringer = 15                         # correct value

number_of_steps_on_domain = 1000

# Program
l_straight_skin = np.sqrt(R**2 + (Ca-R)**2)

# omega1 = np.pi/2        # value between 0 and pi/2
omega1 = np.linspace(0,np.pi/2,number_of_steps_on_domain)
shear_flow_1 = q_b12(S_z, I_yy, n_stringer, Ca, h, t_sk, z_ce, B, omega1)

# y_1 = R                  # value between 0 and R
y_1 = np.linspace(0,R,number_of_steps_on_domain)
shear_flow_2 = q_b42(S_z, I_yy, h, t_sp, z_ce, y_1)

# s1 = l_straight_skin    # value between 0 and l_straight_skin
s1 = np.linspace(0,l_straight_skin,number_of_steps_on_domain)
shear_flow_3 = q_b23(S_z, I_yy, n_stringer, Ca, h , t_sk, z_ce, B, s1)

# s2 = 0                  # value between 0 and l_straight_skin
s2 = np.linspace(0,l_straight_skin,number_of_steps_on_domain)
shear_flow_4 = q_b35(S_z, I_yy, n_stringer, Ca, h, t_sk, z_ce, B, s2)

# y_2 = -R                 # value between -R and 0
y_2 = np.linspace(-R,0,number_of_steps_on_domain)
shear_flow_5 = q_b54(S_z, I_yy, h, t_sp, z_ce, y_2)

# omega2 = -np.pi/2       # value between -pi/2 and 0
omega2 = np.linspace(-np.pi/2,0,number_of_steps_on_domain)
shear_flow_6 = q_b51(S_z, I_yy, n_stringer, Ca, h, t_sk, z_ce, B, omega2)

# print('\nshear flow circ_up: ', shear_flow_1,'\nshear flow spar_up: ', shear_flow_2,'\nshear flow straight_up: ', shear_flow_3,'\nshear flow straight_lo: ', shear_flow_4,'\nshear flow spar_lo: ', shear_flow_5,'\nshear flow circ_lo: ', shear_flow_6)
plt.figure('circ part 1, pos: clock wise')
plt.plot(omega1, shear_flow_1)
plt.show()

plt.figure('spar 1, pos: upwards')
plt.plot(y_1, shear_flow_2)
plt.show()

plt.figure('straight part 1, pos: down to the right')
plt.plot(s1,shear_flow_3)
plt.show()

plt.figure('straight part 2, pos: down to the left')
plt.plot(s2, shear_flow_4)
plt.show()

plt.figure('spar 2, pos: upwards')
plt.plot(y_2, shear_flow_5)
plt.show()

plt.figure('circ part 2, pos: clock wise')
plt.plot(omega2, shear_flow_6)
plt.show()

print()