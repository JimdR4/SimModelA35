''' Still needs to be verified '''

import numpy as np 
# from Z_Coordinates_Stringers import z_coordinates_stringers
# from MOI_and_Centroid import MOI_and_centroid

def MOI_and_centroid():
    import math as m
    h = 24.8    #cm
    c = 51.5    #cm
    r = h/2     #cm
    nstring = 11
    phi = m.atan(r / (c-r))     #rad
    l = 269.1   #cm
    #print(phi)
    #stringer locations!!!

    lskin = m.pi*r + 2*m.sqrt(r*r + (c-r)*(c-r))    #cm
    #print(lskin)

    bstring = lskin/nstring     #cm (stringer spacing)
    a = 11      #number of stringers
    b = 2       #y and z coordinates both need a column
    strloc = np.zeros((a, b))


    for i in range(1, 6):
        #stringer location on circular part
        if i == 1:
            z = r-(r * m.cos((bstring / (m.pi * r/2))*m.pi/2))  
            y = (r * m.sin((bstring / (m.pi * r/2))*m.pi/2))
        #stringer locations on straight part    
        if i >= 2:
            z = r + (2*bstring + bstring*(i-2) - (m.pi * r/2))*m.cos(phi)
            y = r - (2*bstring + bstring*(i-2) - (m.pi * r/2))*m.sin(phi)
            
        strloc[i,0] = z
        strloc[i,1] = y
        #print(strloc)
        
    for i  in range(6,11):
        #stringers on the bottom side, so z coordinate is the same and y coordinate is the same but negative
        strloc[i,0] = strloc[a-i,0]
        strloc[i,1] = -strloc[a-i,1]
    #print(strloc)

    #centroid z coordinate
    #thicknesses and spar area using thinwalled
    tskin = 0.11 #cm
    tspar = 0.22 #cm
    astr = 0.12*(1.5+3) #cm^2
    #area calculations of parts
    acirc = m.pi*r*tskin    #cm^2
    acone = 2* m.sqrt(r*r + (c-r)*(c-r))*tskin  #cm^2
    aspar = h*tspar     #cm^2
    #centroid calculation
    zbar = (2*r/m.pi*acirc+ r * aspar + (r+(c-r)/2)*acone + strloc[1,0]*2*astr + strloc[2,0]*2*astr + strloc[3,0]*2*astr + strloc[4,0]*2*astr + strloc[5,0]*2*astr) / (acirc+ aspar + acone + 11*astr)
    #print(zbar)

    #Iyy
    #component Iyy
    Iyycirc = (r**3)*tskin*m.pi/2 + acirc*(((zbar-(r-(2*r/m.pi)))**2)-((2*r/m.pi)**2))
    Iyyspar = ((zbar-r)**2)*aspar
    Iyycone = (tskin*(m.sqrt(r*r + (c-r)*(c-r)))**3)*(m.cos(phi))/12+(((r+(c-r)/2)-zbar)**2)*acone
    Iyystring = zbar*zbar*astr + ((zbar - strloc[1,0])**2)*2*astr + ((zbar - strloc[2,0])**2)*2*astr +((zbar - strloc[3,0])**2)*2*astr +((zbar - strloc[4,0])**2)*2*astr +((zbar - strloc[5,0])**2)*2*astr
    #total Iyy
    Iyy = Iyycirc + Iyyspar + Iyycone + Iyystring #unit = cm^4
    #print(Iyy)

    #Izz
    #component Izz
    Izzcirc = (r**3)*tskin*m.pi/2
    Izzcone = tskin*((2*m.sqrt(r*r + (c-r)*(c-r)))**3)*((m.sin(phi))**2)/12
    Izzspar = tspar * (h**3)/12
    Izzstringer = (((strloc[1,1])**2)*2*astr) + (((strloc[2,1])**2)*2*astr) + (((strloc[3,1])**2)*2*astr) + (((strloc[4,1])**2)*2*astr) + (((strloc[5,1])**2)*2*astr)
    #total Izz
    Izz = Izzcirc + Izzcone + Izzspar + Izzstringer
    #print(Izz)

    zbar = -1*zbar/100
    Iyy  = Iyy * 10**(-12)
    Izz  = Izz * 10**(-12)

    return Iyy, Izz, zbar

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

def q_b12(S_z, I_yy, n_stringer, Ca, h, t_sk, z_ce, B, omega1):
    ''' Input:      S_z = shear force in z direction
                    I_yy = moment of inertia around the y axis
                    n_stringer = the number of stringers used
                    Ca = lenth of chord
                    h = height of aileron
                    t_sk = skin thickness
                    z_sc = z coordinate of shear center
                    B = boom area
                    omega1 = defined as clockwise positive
        Output:     val = shear flow at omega1
        Remarks:    - This is the shear flow in the circular part starting at the leading edge going clockwise
                    - omega 1 runs from 0 to pi/2
                    - Global coordinate system is defined at leading edge with y-axis pointing upwards and z-axis pointing away from the leading edge
                    - The coordinate system used for the integrations has the same orientation as the global coordinate system, but has its origin the the centroid
                    - Z-axis is assumed to be axis of symmetry
    '''
    if omega1 < 0 or omega1 > np.pi/2:
        raise ValueError('Omega1 should have a value between 0 and pi/2')

    R = h/2
    z = -(1-np.cos(omega1))*R - z_ce
    z_co, str_loc, l_sk = z_coordinates_stringers(Ca, h, n_stringer)
    val = -S_z/I_yy * (-t_sk*R**2*omega1 - t_sk*R**2*np.sin(omega1) + t_sk*z_ce*R*omega1 + 0.5*B*z_co[0])

    for i in range(len(z_co)):
        if str_loc[i] == 'c':
            if abs(z+z_ce) >= abs(z_co[i]): 
                val = val + -S_z/I_yy* B* z_co[i]
        
    return val

def q_b42(S_z, I_yy, h, t_sp, z_ce, y):
    ''' Input:      S_z = shear force in z direction
                    I_yy = moment of inertia around the y axis
                    h = height of aileron
                    t_sp = spar thickness
                    z_sc = z coordinate of shear center
                    y = defined as upwards positive
        Output:     val = shear flow at y
        Remarks:    - This is the shear flow in the spar starting at the chord going upwards
                    - y runs from 0 to R
                    - Global coordinate system is defined at leading edge with y-axis pointing upwards and z-axis pointing away from the leading edge
                    - The coordinate system used for the integrations has the same orientation as the global coordinate system, but has its origin the the centroid
                    - Z-axis is assumed to be axis of symmetry
    '''
    R = h/2
    if y < 0 or y > R:
        raise ValueError('y should have a vale between 0 and R')

    val = -S_z/I_yy * (-R*t_sp*y - z_ce*t_sp*y)
    return val

def q_b23(S_z, I_yy, n_stringer, Ca, h, t_sk, z_ce, B, s1):
    ''' Input:      S_z = shear force in z direction
                    I_yy = moment of inertia around the y axis
                    n_stringer = the number of stringers used
                    Ca = lenth of chord
                    h = height of aileron
                    t_sk = skin thickness
                    z_sc = z coordinate of shear center
                    B = boom area
                    s1 = defined as clockwise positive
        Output:     val = shear flow at s1
        Remarks:    - This is the shear flow in the straight part of the skin starting on top of the aileron going downwards to the right
                    - s1 runs from 0 to the length of the straight part, l_sk
                    - Global coordinate system is defined at leading edge with y-axis pointing upwards and z-axis pointing away from the leading edge
                    - The coordinate system used for the integrations has the same orientation as the global coordinate system, but has its origin the the centroid
                    - Z-axis is assumed to be axis of symmetry
    '''
    R = h/2
    z_co, str_loc, l_sk = z_coordinates_stringers(Ca, h, n_stringer)

    if s1 < 0 or s1 > l_sk:
        raise ValueError('s1 should have a value between 0 and length of the straight skin')

    z = -R-z_ce - (Ca-R)/l_sk * s1
    val = -S_z/I_yy * (-t_sk*R*s1 - t_sk*z_ce*s1 - t_sk*(Ca-R)/(2*l_sk)*s1**2) + q_b12(S_z, I_yy, n_stringer, Ca, h, t_sk, z_ce, B, np.pi/2) + q_b42(S_z, I_yy, h, t_sp, z_ce, R)

    for i in range(len(z_co)):
        if str_loc[i] == 's':
            if abs(z+z_ce) >= abs(z_co[i]): 
                val = val + -S_z/I_yy*B*z_co[i]
        else:
            continue

    return val

def q_b35(S_z, I_yy, n_stringer, Ca, h, t_sk, z_ce, B, s2):
    ''' Input:      S_z = shear force in z direction
                    I_yy = moment of inertia around the y axis
                    n_stringer = the number of stringers used
                    Ca = lenth of chord
                    h = height of aileron
                    t_sk = skin thickness
                    z_sc = z coordinate of shear center
                    B = boom area
                    s1 = defined as clockwise positive
        Output:     val = shear flow at s2
        Remarks:    - This is the shear flow in the straight part of the skin starting at the trailing edge going downwards to the left
                    - s2 runs from 0 to the length of the straight part, l_sk
                    - Global coordinate system is defined at leading edge with y-axis pointing upwards and z-axis pointing away from the leading edge
                    - The coordinate system used for the integrations has the same orientation as the global coordinate system, but has its origin the the centroid
                    - Z-axis is assumed to be axis of symmetry
    '''
    z_co, str_loc, l_sk = z_coordinates_stringers(Ca, h, n_stringer)

    if s2 < 0 or s2 > l_sk:
        raise ValueError('s1 should have a value between 0 and length of the straight skin') 

    return q_b23(S_z, I_yy, n_stringer, Ca, h, t_sk, z_ce, B, l_sk-s2) * -1

def q_b54(S_z, I_yy, h, t_sp, z_ce, y):
    ''' Input:      S_z = shear force in z direction
                    I_yy = moment of inertia around the y axis
                    h = height of aileron
                    t_sp = spar thickness
                    z_sc = z coordinate of shear center
                    y = defined as upwards positive
        Output:     val = shear flow at y
        Remarks:    - This is the shear flow in the spar starting at the bottom of the aileron going upwards
                    - y runs from -R to 0
                    - Global coordinate system is defined at leading edge with y-axis pointing upwards and z-axis pointing away from the leading edge
                    - The coordinate system used for the integrations has the same orientation as the global coordinate system, but has its origin the the centroid
                    - Z-axis is assumed to be axis of symmetry
    '''
    R = h/2
    if y > 0 or y < -R:
        raise ValueError('y should have a vale between -R and 0')

    return q_b42(S_z, I_yy, h, t_sp, z_ce, -1*y) * -1

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
    if omega2 > 0 or omega1 < -np.pi/2:
        raise ValueError('Omega1 should have a value between -pi/2 and 0')

    return q_b12(S_z, I_yy, n_stringer, Ca, h, t_sk, z_ce, B, -1*omega2) * -1

''' FOR TESTING '''
''' Remarks:    Coordinate system at leading edge with y-axis pointing upwards and z-axis pointing away from the aileron
                These shearflows are caused only due to a horizontal shear force S_z(x) which acts along the chord and therefore through the shear center '''

# Input values
S_z = 100
B   = 1
I_yy, I_zz, z_ce = MOI_and_centroid()   # correct values // z_ce is a negative value!!
t_sk = 1.1 * 10**(-3)                   # correct value
t_sp = 2.2 * 10**(-3)                   # correct value
h    = 0.248                            # correct value
R    = h/2                              # correct value
Ca   = 0.515                            # correct value
n_stringer = 11                         # correct value

l_straight_skin = np.sqrt(R**2 + (Ca-R)**2)

omega1 = np.pi/2        # value between 0 and pi/2
shear_flow_1 = q_b12(S_z, I_yy, n_stringer, Ca, h, t_sk, z_ce, B, omega1)

y1 = R                  # value between 0 and R
shear_flow_2 = q_b42(S_z, I_yy, h, t_sp, z_ce, y1)

s1 = l_straight_skin    # value between 0 and l_straight_skin
shear_flow_3 = q_b23(S_z, I_yy, n_stringer, Ca, h , t_sk, z_ce, B, s1)

s2 = 0                  # value between 0 and l_straight_skin
shear_flow_4 = q_b35(S_z, I_yy, n_stringer, Ca, h, t_sk, z_ce, B, s2)

y2 = -R                 # value between -R and 0
shear_flow_5 = q_b54(S_z, I_yy, h, t_sp, z_ce, y2)

omega2 = -np.pi/2       # value between -pi/2 and 0
shear_flow_6 = q_b51(S_z, I_yy, n_stringer, Ca, h, t_sk, z_ce, B, omega2)

print('\nshear flow circ_up: ',shear_flow_1,'\nshear flow spar_up: ', shear_flow_2,'\nshear flow straight_up: ', shear_flow_3,'\nshear flow straight_lo: ', shear_flow_4,'\nshear flow spar_lo: ', shear_flow_5,'\nshear flow circ_lo: ', shear_flow_6)
