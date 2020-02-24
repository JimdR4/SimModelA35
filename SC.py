#import all the variables +numpy here first

#shear centre
#first calculate the closed section shear flows
T=1 
#changing T has an influence on q01 and q02, which is logical, but it does not 
#have an influence on J, which is also logical as J is only dependent on the geometry
def SC(r,tskin,tspar,strloc,c,h,T):
    a1 = Am(c,h)[1]
    a2 = Am(c,h)[2]
    #T = 2*a1*q01 + 2*a2*q02
    #dtheta1 = (1/(tskin*Izz))*(2*tskin*r*r*r) + (1/Izz)*(r*r*r/3) + q01*(m.pi*r/tskin + 2*r/tspar) - q02*2*r/tspar
    #dtheta2 = (-2/(tskin*Izz))*(l*l*r*tskin/3 + tskin*r*r*l*3/2) - (1/Izz)*(r*r*r*3) + q02*(2*l/tskin + 2*r/tspar) - q01*2*r/tspar
    matone = np.array([[(m.pi*r/tskin + 2*r/tspar + 2*r/tspar),(-2*r/tspar - (2*l/tskin + 2*r/tspar))],[2*a1, 2*a2]])
    mattwo = np.array([(-2*T/(tskin*Izz))*(l*l*r*tskin/3 + tskin*r*r*l*3/2) - (T/Izz)*(r*r*r*3) - ((T/(tskin*Izz))*(2*tskin*r*r*r) + (T/Izz)*(r*r*r/3)),T])
    qs = np.linalg.solve(matone, mattwo)
    return (qs)

q01 = SC(r,tskin,tspar,strloc,c,h,T)[0]
q02 = SC(r,tskin,tspar,strloc,c,h,T)[1]
#print(q01, q02)

Gdthetadz = (1/(2*a1))*(m.pi*r*tskin*q01 + 2*r*tspar*(q01-q02))
J = T/ Gdthetadz
#print(J)

###### calculate all the open section shear flows for Vy(x)######
Vy = 1
#changing Vy does have an influence on the location of the SC, but that is due 
#to the fact that q01 and q02 are also in there, and they do not change with Vy. 
#Maybe we should leave them out?
#up to boom 1, so 0<s<bstring
q11 = -(Vy/Izz)*tskin*r*r*-m.cos(s/r)
#from boom 1 to top of spar, so bstring<s<pi*r/2
q1 = -(Vy/Izz)*(tskin*r*r*-m.cos(s/r) + astr*strloc[1,1])+ q11
#from centre of spar to top of spar, so 0<s<r
q2 = -(Vy/Izz)*tspar*s*s/2
#from top of spar to boom 2, so 0<s<(2*bstring-(pi*r/2))
q32 = -(Vy/Izz)*(tskin*r*s-tskin*r*s*s/(2*l))+q1+q2
#from boom 2 to boom 3, so (2*bstring-(pi*r/2))<s<(3*bstring-(pi*r/2))
q33 = -(Vy/Izz)*(tskin*r*s-tskin*r*s*s/(2*l) + astr*strloc[2,1])+q32
#from boom 3 to boom 4, so (3*bstring-(pi*r/2))<s<(4*bstring-(pi*r/2))
q34 = -(Vy/Izz)*(tskin*r*s-tskin*r*s*s/(2*l) + astr*strloc[3,1])+q33
#from boom 4 to boom 5, so (4*bstring-(pi*r/2))<s<(5*bstring-(pi*r/2))
q35 = -(Vy/Izz)*(tskin*r*s-tskin*r*s*s/(2*l) + astr*strloc[4,1])+q34
#from boom 4 to trailing edge, so (5*bstring-(pi*r/2))<s<l
q3 = -(Vy/Izz)*(tskin*r*s-tskin*r*s*s/(2*l) + astr*strloc[5,1])+q35
#from trailing edge to boom 6, so 0<s<(6*bstring-(pi*r/2)-l)
q46 = -(Vy/Izz)*(-tskin*r*s*s/(2*l))+q3
#from boom 6 to boom 7, so (6*bstring-(pi*r/2)-l)<s<(7*bstring-(pi*r/2)-l)
q47 = -(Vy/Izz)*(-tskin*r*s*s/(2*l) + astr*strloc[6,1])+q46
#from boom 7 to boom 8, so (7*bstring-(pi*r/2)-l)<s<(8*bstring-(pi*r/2)-l)
q48 = -(Vy/Izz)*(-tskin*r*s*s/(2*l) + astr*strloc[7,1])+q47
#from boom 8 to boom 9, so (8*bstring-(pi*r/2)-l)<s<(9*bstring-(pi*r/2)-l)
q49 = -(Vy/Izz)*(-tskin*r*s*s/(2*l) + astr*strloc[8,1])+q48
#from boom 9 to bottom of spar, so (9*bstring-(pi*r/2)-l)<s<l
q4 = -(Vy/Izz)*(-tskin*r*s*s/(2*l) + astr*strloc[9,1])+q49
#from centre of spar to bottom of spar, so 0<s<r
q5 = (Vy/Izz)*(tspar*s*s/2)
#from top of spar to boom 10, so 0<s<((pi*r/2)-bstring)
#and from boom 10 to end, so ((pi*r/2)-bstring)<s<(pi*r/2)
#this part is mirrored from the first two, so q610=-q1, and q6=-q11
q610 = -q1
q6 = -q11

#now for the SC:
#up to boom 1, so 0<s<bstring
sb = 0
se = bstring
S11 = -(Vy/Izz)*(-tskin*r*r*r*m.sin(se/r)) - (-(Vy/Izz)*(-tskin*r*r*r*m.sin(sb/r)))
#from boom 1 to top of spar, so bstring<s<pi*r/2
sb = bstring
se = m.pi*r/2
S1 = -(Vy/Izz)*(tskin*r*r*r*-m.sin(se/r) + astr*strloc[1,1]*se)+ q11*se - (-(Vy/Izz)*(tskin*r*r*r*-m.sin(sb/r) + astr*strloc[1,1]*sb)+ q11*sb)
#from top of spar to boom 2, so 0<s<(2*bstring-(pi*r/2))
sb = 0
se = (2*bstring-(m.pi*r/2))
S32 = -(Vy/Izz)*(tskin*r*se*se/2-tskin*r*se*se*se/(6*l))+(q1+q2)*se - (-(Vy/Izz)*(tskin*r*sb*sb/2-tskin*r*sb*sb*sb/(6*l))+(q1+q2)*sb)
#from boom 2 to boom 3, so (2*bstring-(pi*r/2))<s<(3*bstring-(pi*r/2))
sb = (2*bstring-(m.pi*r/2))
se = (3*bstring-(m.pi*r/2))
S33 = -(Vy/Izz)*((tskin*r*se*se/2-tskin*r*se*se*se/(6*l)) + astr*strloc[2,1]*se)+q32*se - (-(Vy/Izz)*((tskin*r*sb*sb/2-tskin*r*sb*sb*sb/(6*l)) + astr*strloc[2,1]*sb)+q32*sb)
#from boom 3 to boom 4, so (3*bstring-(pi*r/2))<s<(4*bstring-(pi*r/2))
sb = (3*bstring-(m.pi*r/2))
se = (4*bstring-(m.pi*r/2))
S34 = -(Vy/Izz)*((tskin*r*se*se/2-tskin*r*se*se*se/(6*l)) + astr*strloc[3,1]*se)+q33*se - (-(Vy/Izz)*((tskin*r*sb*sb/2-tskin*r*sb*sb*sb/(6*l)) + astr*strloc[3,1]*sb)+q33*sb)
#from boom 4 to boom 5, so (4*bstring-(pi*r/2))<s<(5*bstring-(pi*r/2))
sb = (4*bstring-(m.pi*r/2))
se = (5*bstring-(m.pi*r/2))
S35 = -(Vy/Izz)*((tskin*r*se*se/2-tskin*r*se*se*se/(6*l)) + astr*strloc[4,1]*se)+q34*se - (-(Vy/Izz)*((tskin*r*sb*sb/2-tskin*r*sb*sb*sb/(6*l)) + astr*strloc[4,1]*sb)+q34*sb)
#from boom 4 to trailing edge, so (5*bstring-(pi*r/2))<s<l
sb = (2*bstring-(m.pi*r/2))
se = l
S3 = -(Vy/Izz)*((tskin*r*se*se/2-tskin*r*se*se*se/(6*l)) + astr*strloc[5,1]*se)+q35*se - (-(Vy/Izz)*((tskin*r*sb*sb/2-tskin*r*sb*sb*sb/(6*l)) + astr*strloc[5,1]*sb)+q35*sb)
#from trailing edge to boom 6, so 0<s<(6*bstring-(pi*r/2)-l)
sb = 0
se = (6*bstring-(m.pi*r/2)-l)
S46 = -(Vy/Izz)*(-tskin*r*se*se*se/(6*l))+q3*se - (-(Vy/Izz)*(-tskin*r*sb*sb*sb/(6*l))+q3*sb)
#from boom 6 to boom 7, so (6*bstring-(pi*r/2)-l)<s<(7*bstring-(pi*r/2)-l)
sb = (6*bstring-(m.pi*r/2)-l)
se = (7*bstring-(m.pi*r/2)-l)
S47 = -(Vy/Izz)*(-tskin*r*se*se*se/(6*l) + astr*strloc[6,1]*se)+q46*se - (-(Vy/Izz)*(-tskin*r*sb*sb*sb/(6*l) + astr*strloc[6,1]*sb)+q46*sb)
#from boom 7 to boom 8, so (7*bstring-(pi*r/2)-l)<s<(8*bstring-(pi*r/2)-l)
sb = (7*bstring-(m.pi*r/2)-l)
se = (8*bstring-(m.pi*r/2)-l)
S48 = -(Vy/Izz)*(-tskin*r*se*se*se/(6*l) + astr*strloc[7,1]*se)+q47*se - (-(Vy/Izz)*(-tskin*r*sb*sb*sb/(6*l) + astr*strloc[7,1]*sb)+q47*sb)
#from boom 8 to boom 9, so (8*bstring-(pi*r/2)-l)<s<(9*bstring-(pi*r/2)-l)
sb = (8*bstring-(m.pi*r/2)-l)
se = (9*bstring-(m.pi*r/2)-l)
S49 = -(Vy/Izz)*(-tskin*r*se*se*se/(6*l) + astr*strloc[8,1]*se)+q48*se - (-(Vy/Izz)*(-tskin*r*sb*sb*sb/(6*l) + astr*strloc[8,1]*sb)+q48*sb)
#from boom 9 to bottom of spar, so (9*bstring-(pi*r/2)-l)<s<l
sb = (9*bstring-(m.pi*r/2)-l)
se = l
S4 = -(Vy/Izz)*(-tskin*r*se*se*se/(6*l) + astr*strloc[9,1]*se)+q49*se - (-(Vy/Izz)*(-tskin*r*sb*sb*sb/(6*l) + astr*strloc[9,1]*sb)+q49*sb)
#from top of spar to boom 10, so 0<s<((pi*r/2)-bstring)
#and from boom 10 to end, so ((pi*r/2)-bstring)<s<(pi*r/2)
#this part is equal the first two, so S610=S1, and S6=S11
S610 = S1
S6 = S11
#zsc is the distance from the centre of the spar to the SC, so if the 
#SC wrt the origin is desired the distance r has to be added to it
#zsc with q01 and q02
#zsc = ((S11+S1+S610+S6+(q01*m.pi*r))*r+(S32+S33+S34+S35+S3+S46+S47+S48+S49+S4+(q02*m.pi*r))*(r*m.sin(m.pi-phi)))/Vy
#zsc without q01 and q02
zsc = ((S11+S1+S610+S6)*r+(S32+S33+S34+S35+S3+S46+S47+S48+S49+S4)*(r*m.sin(m.pi-phi)))/Vy
print(zsc)
