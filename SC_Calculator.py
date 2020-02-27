# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 18:09:42 2020

@author: cysch
"""

from Intergration import val_simpsons_rule_integration
from ShearflowSy import *
#Geometries

alpha = np.arcsin(R/l_straight_skin)



F12=val_simpsons_rule_integration(omega1,3001,shear_flow_1)
F42=val_simpsons_rule_integration(y_1,3001,shear_flow_1)
F23=val_simpsons_rule_integration(s1,3001,shear_flow_1)
F35=val_simpsons_rule_integration(s2,3001,shear_flow_1)
F54=val_simpsons_rule_integration(y_2,3001,shear_flow_1)
F61=val_simpsons_rule_integration(omega2,3001,shear_flow_1)


A = np.array([[(np.pi*R/t_sk)+2*R/t_sp,-2*R/t_sp],
             [-2*R/t_sk,l_straight_skin/t_sk+2*R/t_sp]])

B = np.array([[-1*((F12/t_sk)-(F42/t_sp)-(F54/t_sp)+(F61/t_sk))],
               [-1*((F23/t_sk)+(F35/t_sk)+(F54/t_sp)+(F42/t_sp))]])

                
C = np.linalg.solve(A,B)

q01 = C[0]
q02 = C[1]

n=(F12+q01*np.pi/2*R)*R+(F61+q01*np.pi/2*R)*R+(F23+q02*l_straight_skin)*R*np.sin(np.pi/2-alpha)+(F35+q02*l_straight_skin)*R*np.sin(np.pi/2-alpha)







