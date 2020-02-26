# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 14:30:55 2020
@author: cysch
"""

import numpy as np
from matplotlib import pyplot as plt

from Y_Coordinates_Stringers import y_coordinates_stringers
from Z_Coordinates_Stringers import z_coordinates_stringers
from Equations_needed_for_stress import Sz_x, Sy_x


def q_b12y(R1y,RIy,R2y,P,R3y,x, I_zz, nst, Ca, ha, tsk, zsc, B, omega1):
    if np.size(omega1) > 1:
        values = []
        for omega in omega1:
            if omega < 0 or omega > np.pi/2:
                raise ValueError('Omega1 should have a value between 0 and pi/2')
            r = ha/2
            z = -(1-np.cos(omega))*r - zsc
            # y = h*np.sin(omega)
            z_co, str_loc = z_coordinates_stringers(Ca, ha, nst)[:2]
            y_co, str_loc = y_coordinates_stringers(Ca, ha, nst)[:2]

            val = -Sy_x(R1y,RIy,R2y,P,R3y,x)/I_zz * (-tsk*r**2*np.cos(omega) + tsk*r**2 + 0.5*B*y_co[0]) 

            z_glob_c = z + zsc

            for i in range(len(z_co)):
                if str_loc[i] == 'c':
                    if abs(z_glob_c) >= abs(z_co[i]): 
                        val = val + -Sy_x(R1y,RIy,R2y,P,R3y,x)/I_zz*B*y_co[i]
            
            values.append(val)

    else:
        if omega1 < 0 or omega1 > np.pi/2:
            raise ValueError('Omega1 should have a value between 0 and pi/2')
            
        r = h/2
        z = -(1-np.cos(omega1))*r - zsc
        # y = h*np.sin(omega1)
        z_co, str_loc = z_coordinates_stringers(Ca, ha, nst)[:2]
        y_co, str_loc = y_coordinates_stringers(Ca, ha, nst)[:2]

        values = -Sy_x(R1y,RIy,R2y,P,R3y,x)/I_zz * (-tsk*r**2*np.cos(omega1) + tsk*r**2 + 0.5*B*y_co[0]) 

        z_glob_co = z + zsc

        for i in range(len(z_co)):
            if str_loc[i] == 'c':
                if abs(z_glob_co) >= abs(z_co[i]): 
                    values = values + -Sy_x(R1y,RIy,R2y,P,R3y,x)/I_zz* B* (y_co[i])

    return values 

def q_b42y(R1y,RIy,R2y,P,R3y,x, I_zz, ha, tsp, zsc, y1):
    ''' Input:      R1y,RIy,R2y,P,R3y,x = shear force in y direction inputs
                    I_zz = moment of inertia around the y axis
                    ha = height of aileron
                    tsp = spar thickness
                    zsc = z coordinate of shear center
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
        for y in y1:
            r = ha/2
            if y < 0 or y > r:
                raise ValueError('y should have a vale between 0 and R')

            val = -Sy_x(R1y,RIy,R2y,P,R3y,x)/I_zz * ((tsp*y**2)/2)

            values.append(val)

    else:
        r = ha/2
        if y1 < 0 or y1 > r:
            raise ValueError('y should have a vale between 0 and R')

        values = -Sy_x(R1y,RIy,R2y,P,R3y,x)/I_zz * ((tsp*y1**2)/2)

    return values 

def q_b23y(R1y,RIy,R2y,P,R3y,x, I_zz, nst, Ca, ha, tsk, tsp, zsc, B, s1):
    ''' Input:      R1y,RIy,R2y,P,R3y,x = shear force in y direction inputs
                    I_zz = moment of inertia around the y axis
                    nst = the number of stringers used
                    Ca = lenth of chord
                    ha = height of aileron
                    tsk = skin thickness
                    zsc = z coordinate of shear center
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

        for s in s1:
            r = ha/2
            y_co, str_loc, l_sk = y_coordinates_stringers(Ca, ha, nst)
            z_co, str_loc, l_sk = z_coordinates_stringers(Ca, ha, nst)


            if s < 0 or s > l_sk:
                raise ValueError('s1 should have a value between 0 and length of the straight skin')
            z = -r-zsc - (Ca-r)/l_sk * s
            y = r-(r/l_sk)*s
            val = -Sy_x(R1y,RIy,R2y,P,R3y,x)/I_zz * (tsk*r*s - tsk*r/(2*l_sk)*s**2) + q_b12y(R1y,RIy,R2y,P,R3y,x, I_zz, nst, Ca, ha, tsk, zsc, B, np.pi/2) + q_b42y(R1y,RIy,R2y,P,R3y,x, I_zz, ha, tsp, zsc, r)

            z_glob_c = z+zsc

            for i in range(len(z_co)):
                if str_loc[i] == 's':
                    if abs(z_glob_c) >= abs(z_co[i]): 
                        val = val - Sy_x(R1y,RIy,R2y,P,R3y,x)/I_zz*B*(y_co[i])
            
            values.append(val)

    else:
        r = ha/2
        z_co, str_loc, l_sk = z_coordinates_stringers(Ca, ha, nst)
        y_co, str_loc, l_sk = y_coordinates_stringers(Ca, ha, nst)


        if s1 < 0 or s1 > l_sk:
            raise ValueError('s1 should have a value between 0 and length of the straight skin')
        
        y = r-(r/l_sk)*s1
        z = -r-zsc - (Ca-r)/l_sk * s1

        values = -Sy_x(R1y,RIy,R2y,P,R3y,x)/I_zz * (tsk*r*s1 - tsk*r/(2*l_sk)*s1**2) + q_b12y(R1y,RIy,R2y,P,R3y,x, I_zz, nst, Ca, ha, tsk, zsc, B, np.pi/2) + q_b42y(R1y,RIy,R2y,P,R3y,x, I_zz, ha, tsp, zsc, r)

        z_glob_co = z+zsc

        for i in range(len(z_co)):
            if str_loc[i] == 's':
                if abs(z_glob_co) >= abs(z_co[i]): 
                    values = values + -Sy_x(R1y,RIy,R2y,P,R3y,x)/I_zz*B*(y_co[i])

    return values

def q_b35y(R1y,RIy,R2y,P,R3y,x, I_zz, nst, Ca, ha, tsk, zsc, B, s2):
    ''' Input:      R1y,RIy,R2y,P,R3y,x = shear force in y direction inputs
                    I_zz = moment of inertia around the y axis
                    nst = the number of stringers used
                    Ca = lenth of chord
                    ha = height of aileron
                    tsk = skin thickness
                    zsc = z coordinate of shear center
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
    if np.size(s2) > 1:
        values = []

        for s in s2:
            r = ha/2
            y_co, str_loc, l_sk = y_coordinates_stringers(Ca, ha, nst)
            z_co, str_loc, l_sk = z_coordinates_stringers(Ca, ha, nst)

            y_co = y_co*-1

            if s < 0 or s > l_sk:
                raise ValueError('s1 should have a value between 0 and length of the straight skin')

            # y = -(R/l_sk)*s
            z = -Ca-zsc + (Ca-r)/l_sk * s

            val = -Sy_x(R1y,RIy,R2y,P,R3y,x)/I_zz * (-tsk*r/(2*l_sk)*s**2) + q_b23y(R1y,RIy,R2y,P,R3y,x, I_zz, nst, Ca, ha, tsk, zsc, B, l_sk)

            z_glob_c = z+zsc

            for i in range(len(z_co)):
                if str_loc[i] == 's':
                    if abs(z_glob_c) <= abs(z_co[i]): 
                        val = val + -Sy_x(R1y,RIy,R2y,P,R3y,x)/I_zz*B*(y_co[i])
            
            values.append(val)

    else:
        r = ha/2
        z_co, str_loc, l_sk = z_coordinates_stringers(Ca, ha, nst)
        y_co, str_loc, l_sk = y_coordinates_stringers(Ca, ha, nst)

        y_co = -1*y_co

        if s2 < 0 or s2 > l_sk:
            raise ValueError('s2 should have a value between 0 and length of the straight skin')
        
        # y = -(R/l_sk)*s2
        z = -Ca-zsc + (Ca-r)/l_sk * s2

        values = -Sy_x(R1y,RIy,R2y,P,R3y,x)/I_zz * (-tsk*r/(2*l_sk)*s2**2) + q_b23y(R1y,RIy,R2y,P,R3y,x, I_zz, nst, Ca, ha, tsk, zsc, B, l_sk)

        z_glob_co = z+zsc

        for i in range(len(z_co)):
            if str_loc[i] == 's':
                if abs(z_glob_co) <= abs(z_co[i]): 
                    values = values + -Sy_x(R1y,RIy,R2y,P,R3y,x)/I_zz*B*(y_co[i])

    return values

def q_b54y(R1y,RIy,R2y,P,R3y,x, I_zz, ha, tsp, zsc, y2):
    ''' Input:      R1y,RIy,R2y,P,R3y,x = shear force in y direction inputs
                    I_zz = moment of inertia around the y axis
                    ha = height of aileron
                    tsp = spar thickness
                    zsc = z coordinate of shear center
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
    r = ha/2
    if np.size(y2) > 1:
        values=[]
        for y in y2:
            if y > 0 or y < -r:
                raise ValueError('y should have a vale between -R and 0')
                
            val = -Sy_x(R1y,RIy,R2y,P,R3y,x)/I_zz * ((tsp*y**2)/2 )
            
            values.append(val)
        
    else:
        if y2 > 0 or y2 < -r:
            raise ValueError('y should have a vale between -R and 0')
        
        values = -Sy_x(R1y,RIy,R2y,P,R3y,x)/I_zz * ((tsp*y2**2)/2)


    return values

def q_b51y(R1y,RIy,R2y,P,R3y,x, I_zz, nst, Ca, ha, tsk, tsp, zsc, l, y_ce, B, omega2):
    if np.size(omega2) > 1:
        values = []
        
        for omega in omega2:
            if omega > 0 or omega < -np.pi/2:
                raise ValueError('Omega2 should have a value between 0 and pi/2')
            r = ha/2
            z = -(1-np.cos(omega))*r - zsc
            # y = h*np.sin(omega)
            z_co, str_loc = z_coordinates_stringers(Ca, ha, nst)[:2]
            y_co, str_loc = y_coordinates_stringers(Ca, ha, nst)[:2]

            y_co = y_co*-1

            val = -Sy_x(R1y,RIy,R2y,P,R3y,x)/I_zz * (-tsk*r**2*np.cos(omega)) + q_b35y(R1y,RIy,R2y,P,R3y,x, I_zz, nst, Ca, ha, tsk, zsc, B, l) - q_b54y(R1y,RIy,R2y,P,R3y,x, I_zz, ha, tsp, zsc, -r)

            z_glob_c = z + zsc

            for i in range(len(z_co)):
                if str_loc[i] == 'c':
                    if abs(z_glob_c) <= abs(z_co[i]): 
                        val = val + -Sy_x(R1y,RIy,R2y,P,R3y,x)/I_zz* B* (y_co[i])
            
            values.append(val)

    else:
        if omega2 < 0 or omega2 > np.pi/2:
            raise ValueError('Omega1 should have a value between 0 and pi/2')

        r = ha/2
        y = ha*np.sin(omega)
        z_co, str_loc = z_coordinates_stringers(Ca, ha, nst)[:2]
        y_co, str_loc = y_coordinates_stringers(Ca, ha, nst)[:2]

        y_co = y_co**-1

        values = -Sy_x(R1y,RIy,R2y,P,R3y,x)/I_zz * (-tsk*r**2*np.cos(omega))

        z_glob_co = z + zsc

        for i in range(len(z_co)):
            if str_loc[i] == 'c':
                if abs(z_glob_co) <= abs(z_co[i]): 
                    values = values + -Sy_x(R1y,RIy,R2y,P,R3y,x)/I_zz* B* (y_co[i])

    return values 

"""
''''''''' FOR TESTING '''''''''
''' Remarks:    Coordinate system at leading edge with y-axis pointing upwards and z-axis pointing away from the aileron
                These shearflows are caused only due to a horizontal shear force S_z(x) which acts along the chord and therefore through the shear center '''

# Only change these input values !!!
S_y = 1
B   = 5.4*10**(-5) # going from dm**2 to m**2
I_yy = 0.000014222
z_ce = -0.207287029
t_sk = 1.1 * 10**(-3)                   
t_sp = 2.2 * 10**(-3)                   
h    = 0.248                            
R    = h/2                              
Ca   = 0.515                            
n_stringer = 11                         
l_straight_skin = np.sqrt(R**2 + (Ca-R)**2)

number_of_steps_on_domain = 1000


''' Program for single values 
omega1 = np.pi/2        # value between 0 and pi/2
shear_flow_1 = q_b12(S_y, I_zz, n_stringer, Ca, h, t_sk, z_ce, B, omega1)
y_1 = R                  # value between 0 and R
shear_flow_2 = q_b42(S_y, I_zz, h, t_sp, z_ce, y_1)
s1 = l_straight_skin    # value between 0 and l_straight_skin
shear_flow_3 = q_b23(S_y, I_zz, n_stringer, Ca, h , t_sk, z_ce, B, s1)
s2 = 0                  # value between 0 and l_straight_skin
shear_flow_4 = q_b35(S_y, I_zz, n_stringer, Ca, h, t_sk, z_ce, B, s2)
y_2 = -R                 # value between -R and 0
shear_flow_5 = q_b54(S_y, I_zz, h, t_sp, z_ce, y_2)
omega2 = -np.pi/2       # value between -pi/2 and 0
shear_flow_6 = q_b51(S_y, I_zz, n_stringer, Ca, h, t_sk, z_ce, B, omega2)
print('\nshear flow circ_up: ', shear_flow_1,'\nshear flow spar_up: ', shear_flow_2,'\nshear flow straight_up: ', shear_flow_3,'\nshear flow straight_lo: ', shear_flow_4,'\nshear flow spar_lo: ', shear_flow_5,'\nshear flow circ_lo: ', shear_flow_6)
'''


''' Program for multiple values '''

omega1 = np.linspace(0,np.pi/2,number_of_steps_on_domain)
shear_flow_1 = q_b12(S_y, I_yy, n_stringer, Ca, h, t_sk, z_ce, B, omega1)

y_1 = np.linspace(0,R,number_of_steps_on_domain)
shear_flow_2 = q_b42(S_y, I_yy, h, t_sp, z_ce, y_1)
s1 = np.linspace(0,l_straight_skin,number_of_steps_on_domain)
shear_flow_3 = q_b23(S_y, I_yy, n_stringer, Ca, h , t_sk, z_ce, B, s1)

s2 = np.linspace(0,l_straight_skin,number_of_steps_on_domain)
shear_flow_4 = q_b35(S_y, I_yy, n_stringer, Ca, h, t_sk, z_ce, B, s2)

y_2 = np.linspace(0,-R,number_of_steps_on_domain)
shear_flow_5 = q_b54(S_y, I_yy, h, t_sp, z_ce, y_2)

omega2 = np.linspace(-np.pi/2,0,number_of_steps_on_domain)
shear_flow_6 = q_b51(S_y, I_yy, n_stringer, Ca, h, t_sk, z_ce, B, omega2)

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
"""
