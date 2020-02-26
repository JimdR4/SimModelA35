from Z_Coordinates_Stringers import z_coordinates_stringers
from Y_Coordinates_Stringers import y_coordinates_stringers
from mooie_centroid_calcs import Centroid_z
import math as m

def I_yy(Ca,ha,nst,tst,tsk,tsp,hst,wst):
    r = ha/2
    l = m.sqrt(r*r + (Ca-r)*(Ca-r))
    phi = m.atan(r / (Ca-r))
    b_string = tst * (hst + wst)
    acirc = m.pi*r*tsk
    aspar = ha*tsp
    adia = l*tsk
    strlocz = z_coordinates_stringers(Ca, ha, nst)[3]
    #strlocy = y_coordinates_stringers(Ca, ha, nst)[3]
    zc = Centroid_z(Ca,ha,nst,tst,tsk,tsp,hst,wst)
    #contributions from the stringers
    Iyy = 0
    for i in range(0,nst):
        zdist = strlocz[i]-zc
        Iyy = Iyy + b_string*zdist*zdist
    #contributions of other parts
    Iyycirc = (r**3)*tsk*m.pi/2 + acirc*(((-zc-(r-(2*r/m.pi)))**2)-((2*r/m.pi)**2))
    Iyyspar = ((-zc-r)**2)*aspar
    Iyycone = 2*(tsk*(l**3)*(m.cos(phi)))/12+(((r+(Ca-r)/2)+zc)**2)*adia*2
    #add them together
    Iyytot = Iyy + Iyycirc + Iyyspar + Iyycone
    
    return Iyytot

def I_zz(Ca,ha,nst,tst,tsk,tsp,hst,wst):
    r = ha/2
    l = m.sqrt(r*r + (Ca-r)*(Ca-r))
    phi = m.atan(r / (Ca-r))
    b_string = tst * (hst + wst)
    #strlocz = z_coordinates_stringers(Ca, ha, nst)[3]
    strlocy = y_coordinates_stringers(Ca, ha, nst)[3]
    #zc = Centroid_z(Ca,ha,nst,tst,tsk,tsp,hst,wst)
    #contributions from the stringers
    Izz = 0
    for i in range(0,nst):
        ydist = strlocy[i]
        Izz = Izz + b_string*ydist*ydist
    #contributions of other parts
    Izzcirc = (r**3)*tsk*m.pi/2
    Izzcone = tsk*((2*l)**3)*((m.sin(phi))**2)/12
    Izzspar = tsp * (ha**3)/12
    #add them together
    Izztot = Izz + Izzcirc + Izzspar + Izzcone
    
    return Izztot
