'''
Created on Fri Feb 21 2020
@author: Frans Vossen
'''
import numpy as np 

def z_coordinates_stringers(Ca, h, n_stringer):
    """ Input:      Ca = chord length
                    h  = aileron height
                    n_stringer = number of stringers used (ONLY ODD NUMBERS)
        Output:     strloc = the z coordinates for the stringers on the leading edge and top/bottom
                    str_on_circ = number of stringers on the circular part (except leading edge)
                    str_on_straight = number of stringers on the straight part
        Remarks:    length of strloc = (number of stringer on the circular part) + (number of stringers on the straight part) + (string on leading edge)
                    function works only for an odd number of stringers
    """
    R = h/2

    l_circ     = np.pi*R*0.5
    l_straight = np.sqrt(R**2 + (Ca - R)**2)
    l_skin     = l_straight*2 + l_circ*2
    d_stringer = l_skin/n_stringer

    str_on_circ     = int(l_circ/d_stringer)
    str_on_straight = int((n_stringer)/2 - str_on_circ)
    strloc          = np.zeros(str_on_circ+str_on_straight+1)
    strname         = np.empty(str_on_circ+str_on_straight+1,dtype=str)

    for i in range(1,str_on_circ+1):
        l1 = np.cos( (d_stringer*i)/l_circ * np.pi/2 )*R
        z = R - l1
        strloc[i] = -1*z
        strname[i] = 'c'

    for i in range(str_on_straight):
        z = Ca - (l_circ + l_straight - ( (str_on_circ+1)*d_stringer + i*d_stringer ))*(Ca-R)/l_straight
        strloc[i+str_on_circ+1] = -1*z
        strname[i+str_on_circ+1] = 's'

    return strloc, strname, l_straight
