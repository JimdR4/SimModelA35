# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 16:31:52 2020

@author: adria
"""



ha = 0.248
Ca = 0.515

import numpy as np
from matplotlib import pyplot as plt


#Plotting the airfoil profile

z = np.arange(-0.5*ha,0,0.00001) #start,stop,step

#First Semicircle
y1 = np.sqrt((0.5*ha)**2-(z+0.5*ha)**2)

#Second semicircle
y2 = -y1

#Top triangular part
z1 = np.arange(-Ca,-0.5*ha,0.00001) #start,stop,step
y3 = 0.5*ha*z1*(Ca-0.5*ha)**(-1) + 0.5*ha*Ca*(Ca-0.5*ha)**(-1)

#Lower triangular part
y4 = -y3

#Spar 
y_coord = np.arange(-0.5*ha,0.5*ha,0.001)
z4 = -0.5*ha*np.ones(len(y_coord))



#Plotting stuff

plt.plot(z,y1,'c', linewidth=3.3)
plt.plot(z,y2,'c', linewidth=3.3)
plt.plot(z1,y3,'c', linewidth=3.3)
plt.plot(z1,y4,'c', linewidth=3.3)
plt.plot(z4,y_coord,'c', linewidth=3.3)

plt.title('Airfoil')
plt.xlim(0.01,-1.1*Ca)
plt.xlabel('z [m]')
plt.ylabel('y [m]')
plt.grid()
plt.show()

