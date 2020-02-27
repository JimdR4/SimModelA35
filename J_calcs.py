import math as m
import numpy as np

def J(ha,a1,a2,tsp,tsk,l):
    TT = 1
    matrix1 = np.array([[2*a1,2*a2],
                    [(m.pi*(0.5*ha)/(2*a1*tsk)+((0.5*ha)/(a1*tsp))+((0.5*ha)/(a2*tsp))),(-((0.5*ha)/(a1*tsp))-l/(a2*tsk)-(0.5*ha)/(a2*tsp))]])
    m2 = np.array([[TT],
                [0]])

    solu = np.linalg.solve(matrix1,m2)

    q01 = solu[0]
    q02 = solu[1]
  
    Gdthetadz1 = (1/(2*a1))*(m.pi*(0.5*ha)/tsk*q01 + 2*(0.5*ha)/tsp*(q01-q02))
    Gdthetadz2 = (1/(2*a2))*(2*l/tsk*q02+2*(0.5*ha)/tsp*(q02-q01))

    J = (TT/ Gdthetadz2)[0]
    return J
