import numpy as np
import math as m
from matplotlib import pyplot as plt
from mooie_I import I_zz,I_yy

def yz_12(ha, omega1):
    r = ha/2
    #make empty array
    y = np.zeros((len(omega1)))
    #fill the empty array with y coordinates
    for i in range(len(omega1)):
        y[i] = -r*m.sin(omega1[i])
    
    z = np.zeros((len(omega1)))
    for i in range(len(omega1)):
        z[i] = -r+r*m.cos(omega1[i])
    return y, z
        
def yz_42(ha, y_1):
    r = ha/2
    #make empty array
    y = np.zeros((len(y_1)))
    #fill the empty array with y coordinates
    for i in range(len(y_1)):
        y[i] = -y_1[i]
    
    z = np.zeros((len(y_1)))
    for i in range(len(y_1)):
        z[i] = -r
    return y, z

def yz_23(ha, l, Ca, s1):
    r = ha/2
    #make empty array
    y = np.zeros((len(s1)))
    #fill the empty array with y coordinates
    for i in range(len(s1)):
        y[i] = -r+s1[i]*(r/l)
    
    z = np.zeros((len(s1)))
    for i in range(len(s1)):
        z[i] = -r-s1[i]*(Ca-r)/l
    return y, z

def yz_35(ha, l, Ca, s2):
    r = ha/2
    #make empty array
    y = np.zeros((len(s2)))
    #fill the empty array with y coordinates
    for i in range(len(s2)):
        y[i] = s2[i]*(r/l)
    
    z = np.zeros((len(s2)))
    for i in range(len(s2)):
        z[i] = -Ca+s2[i]*(Ca-r)/l
    return y, z

def yz_54(ha, y_2):
    r = ha/2
    #make empty array
    y = np.zeros((len(y_2)))
    #fill the empty array with y coordinates
    for i in range(len(y_2)):
        y[i] = y_2[i]+r
    
    z = np.zeros((len(y_2)))
    for i in range(len(y_2)):
        z[i] = -r    
    return y, z

def yz_51(ha, omega2):
    r = ha/2
    #make empty array
    y = np.zeros((len(omega2)))
    #fill the empty array with y coordinates
    for i in range(len(omega2)):
        y[i] = -r*m.sin(omega2[i])
    
    z = np.zeros((len(omega2)))
    for i in range(len(omega2)):
        z[i] = -r+r*m.cos(omega2[i])   
    return y, z
    
#input values(some are not needed)    
Ca = 0.515  # m
la = 2.691  # m
x1 = 0.174  # m
x2 = 1.051  # m
x3 = 2.512  # m
xa = 0.30   # m
ha = 0.248  # m
tsk = 1.1/1000  # m
tsp = 2.2/1000  # m
tst = 1.2/1000  # m
hst = 15./1000   # m
wst = 30./1000   # m
nst = 11  # -
d1 = 0.01034  # m
d3 = 0.02066  # m
theta = m.radians(25)  # rad
P = 20.6*1000  # N

l = m.sqrt(r**2 + (Ca-r)**2)
number_of_steps_on_domain = 100

#define domains
omega1 = np.linspace(0,np.pi/2,number_of_steps_on_domain)
y_1 = np.linspace(0,r,number_of_steps_on_domain)
s1 = np.linspace(0,l,number_of_steps_on_domain)
s2 = np.linspace(0,l_,number_of_steps_on_domain)
y_2 = np.linspace(-r,0,number_of_steps_on_domain)
omega2 = np.linspace(-np.pi/2,0,number_of_steps_on_domain)

#plot to check
plt.plot(yz_12(ha, omega1)[1],yz_12(ha, omega1)[0])
plt.plot(yz_42(ha, y_1)[1],yz_42(ha, y_1)[0])
plt.plot(yz_23(ha, l, Ca, s1)[1],yz_23(ha, l, Ca, s1)[0])
plt.plot(yz_35(ha, l, Ca, s2)[1],yz_35(ha, l, Ca, s2)[0])
plt.plot(yz_54(ha, y_2)[1],yz_54(ha, y_2)[0])
plt.plot(yz_51(ha, omega2)[1],yz_51(ha, omega2)[0])
plt.xlim(0.01,-1.1*Ca)
plt.ylim(1.1*r,-1.1*r)
plt.show()
