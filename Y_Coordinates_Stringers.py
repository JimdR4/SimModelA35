import numpy as np 

def y_coordinates_stringers(Ca, ha, nst):
    """ Input:      Ca = chord length
                    ha  = aileron height
                    nst = number of stringers used (ONLY ODD NUMBERS)
        Output:     strloc = the y coordinates for the stringers on the leading edge and top/bottom
                    str_on_circ = number of stringers on the circular part (except leading edge)
                    str_on_straight = number of stringers on the straight part
        Remarks:    length of strloc = (number of stringer on the circular part) + (number of stringers on the straight part) + (string on leading edge)
                    function works only for an odd number of stringers
    """
    r = ha/2

    l_circ     = np.pi*r*0.5
    l_straight = np.sqrt(r**2 + (Ca - r)**2)
    l_skin     = l_straight*2 + l_circ*2
    d_stringer = l_skin/n_stringer

    str_on_circ     = int(l_circ/d_stringer)
    str_on_straight = int((n_stringer)/2 - str_on_circ)
    strloc          = np.zeros(str_on_circ+str_on_straight+1)
    strlocc         = np.zeros(nst)
    strname         = np.empty(str_on_circ+str_on_straight+1,dtype=str)

    for i in range(1,str_on_circ+1):
        l1 = np.sin( (d_stringer*i)/l_circ * np.pi/2 )*R
        y = l1
        strloc[i] = -1*y
        strlocc[i] = -1*y
        strname[i] = 'c'

    for i in range(str_on_straight):
        y = (l_circ + l_straight - ( (str_on_circ+1)*d_stringer + i*d_stringer ))*(R)/l_straight        
        strloc[i+str_on_circ+1] = -1*y
        strlocc[i+str_on_circ+1] = -1*y
        strname[i+str_on_circ+1] = 's'
    
    for i in range(1,(int(nst/2)+1)):
        strlocc [nst-i] = strloc [i]

    return strloc, strname, l_straight, strlocc
