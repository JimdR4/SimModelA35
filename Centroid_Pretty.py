import math as m
from Z_Coordinates_Stringers import z_coordinates_stringers

def Centroid_z(Ca,ha,nst,tst,tsk,tsp,hst,wst):
    r = ha/2
    l = m.sqrt(r*r + (Ca-r)*(Ca-r))
    b_string = tst * (hst + wst)
    acirc = m.pi*r*tsk
    aspar = ha*tsp
    adia = l*tsk
    strlocc = z_coordinates_stringers(Ca, ha, nst)[3]
    sumAz = 0
    sumA = 0
    #stringer contributions
    for i in range(0,nst):
        sumAz = sumAz + (strlocc[i]*b_string)
        sumA = sumA + b_string
    #contributions of other parts
    #in order, circular part, spar, cone
    sumAz = -sumAz + ((r-2*r/m.pi)*acirc) + r*aspar + (((Ca-r)/2)+r)*adia*2
    sumA = sumA + acirc + aspar + adia*2
    
    zc = -sumAz/sumA
    #print(zc)
    return zc

#zc = Centroid_z(Ca,ha,nst,tst,tsk,tsp,hst,wst)
#print(zc)
