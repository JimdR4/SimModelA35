import numpy as np
import math

#import data set
dataset = np.loadtxt(r'C:\Users\jimri\Downloads\Studie\VerificationNValidation\aerodynamicloaddo228.dat',dtype='float', delimiter=',')

''' Input variables '''
Ca = 0.515  # [m] / chord length 
la = 2.691  # [m] / span length
Nx = 41     # number of spanwise segments
Nz = 81     # number of chordwise segments
Force_z = dataset.T #transpose is Airodynamic force in Z direction
Force_x = dataset 

#Create to store the unequally spaced nodes
Xnodes = np.zeros(Nx)
Znodes = np.zeros(Nz)


def angle(i,N):
    th =  (i-1)/N*math.pi
    return th

#### Implementing the node spacing formulas in the lists
for i in range(1, Nz+1):
    Znodes[i-1] =  -Ca/4*(2 - math.cos(angle(i, Nz)) - math.cos(angle(i+1, Nz)))
for i in range(1, Nx+1):
    Xnodes[i-1] = la/4*(2 - math.cos(angle(i, Nx)) - math.cos(angle(i+1, Nx)))


#All data is structured, lets start with those cubic splinessss :), for this the second sammery for Applied Numerical Analysis is used#
def interpolation(F, axis):
    Sol = np.zeros((len(F)-1, 4, len(F.T)-1)) #all splines are stored in here

    for r in range(len(F)-1):

        n = len(axis)-1 #number of needed splines 
        h = np.zeros(n) #the stepsize between nodes will be stored in here
        a = np.zeros(n) #a,b,c need to be solved to create our spline
        b = np.zeros(n)
        c = np.zeros(n)
        d = F[r,:] #d is always equal to the force
        
        
        for i in range(0,n):           
            h[i] = (axis[i+1] - axis[i])
        
        S = np.zeros((n+1, n+1)) #the matrix to store the equations
        S[0,0] = 1                  #top right and bottom left corner are  non zero
        S[n,n] = 1
        for j in range(1,n):        #sammary used to solve all values
            S[j,j] = (h[j-1] + h[j])/3
            S[j,j-1] = h[j-1]/6
            S[j,j+1] = h[j]/6
        
            
        #To be able to solve S*c = v for c, we create v
        v = np.zeros(n+1) #add one because we substracted one earlier  
        
        for k in range(0,n-1):
            v[k+1] = 3*((d[k+2]-d[k+1])/h[k+1] - (d[k+1]-d[k])/h[k])
        
        #solve S*b = v, v are the M values, the solutions to the derivatives of the first and second order
        b = np.linalg.solve(S, v)
        
        #Solve for a and c   
        for l in range(0, n):
            c[l] = (d[l+1]-d[l])/h[l] - h[l]*(2*b[l]+b[l+1])/3        #c variable solved
            a[l] = (b[l+1] - b[l])/(3*h[l])                             #a variable solved
        #Remove last element to have 80 intervals
        b = b[:-1]
        d = d[:-1]
        
        #store the solutions to Sol in the proper place
        Sol[r, 0, :] = a
        Sol[r, 1, :] = b
        Sol[r, 2, :] = c
        Sol[r, 3, :] = d         
        
    return Sol



Sol_z = interpolation(Force_z,Znodes )
Sol_x = interpolation(Force_x, Xnodes)






    