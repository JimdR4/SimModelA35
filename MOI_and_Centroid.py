import math as m
import numpy as np

def MOI_and_centroid():
    h = 24.8    #cm
    c = 51.5    #cm
    r = h/2     #cm
    nstring = 11
    phi = m.atan(r / (c-r))     #rad
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