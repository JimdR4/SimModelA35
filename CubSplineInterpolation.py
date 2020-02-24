import numpy as np

def aero_nodes(Ca, la, Nz, Nx):

    Xnodes = np.zeros(Nx)
    Znodes = np.zeros(Nz)

    def angle(i,N):
        th =  (i-1)/N*np.pi
        return th

    for i in range(1, Nz+1):
        Znodes[i-1] =  -Ca/4*(2 - np.cos(angle(i, Nz)) - np.cos(angle(i+1, Nz)))
    for i in range(1, Nx+1):
        Xnodes[i-1] = la/4*(2 - np.cos(angle(i, Nx)) - np.cos(angle(i+1, Nx)))

    return Znodes, Xnodes

def interpolation(F, axis):
    ''' Input:  F = 1-D array containing all function values of the data points that needs to be interpolated
                axis = 1-D array containing all grid points
        Output: Coeff = 4xlen(axis) array containing all coefficients for the splines
        Remarks:      s_i = a*(x-x_i)**3 + b*(x-x_i)**2 + c*(x-x_i) + d
                      Coeff[0,:] all 'a' coefficients
                      Coeff[1,:] all 'b' coefficients
                      Coeff[2,:] all 'c' coefficients
                      Coeff[3,:] all 'd' coefficients
    '''
    Coeff = np.zeros((4, len(F.T)-1)) #all splines are stored in here

    n = len(axis)-1 # number of needed splines 
    h = np.zeros(n) # the stepsize between nodes will be stored in here
    a = np.zeros(n) # a,b,c need to be solved to create our spline
    b = np.zeros(n)
    c = np.zeros(n)
    d = F        # d is always equal to the force
        
    for i in range(n):           
        h[i] = (axis[i+1] - axis[i])
    
    S = np.zeros((n+1, n+1)) #the matrix to store the equations
    S[0,0] = 1               #top right and bottom left corner are  non zero
    S[n,n] = 1
    for j in range(1,n):     #sammary used to solve all values
        S[j,j] = (h[j-1] + h[j])/3
        S[j,j-1] = h[j-1]/6
        S[j,j+1] = h[j]/6
    
    #To be able to solve S*c = v for c, we create v
    v = np.zeros(n+1) #add one because we substracted one earlier  
    
    for k in range(len(axis)-2):
        v[k] = (d[k+2] - d[k+1])/h[k+1] - (d[k+1] - d[k])/h[k]

    # for k in range(n-1):
    #     v[k+1] = 3*((d[k+2]-d[k+1])/h[k+1] - (d[k+1]-d[k])/h[k])
    
    # solve S*b = v, v are the M values, the solutions to the derivatives of the first and second order
    b = np.linalg.solve(S, v)
    
    #Solve for a and c   
    for l in range(n):
        # c[l] = (d[l+1]-d[l])/h[l] - h[l]*(2*b[l]+b[l+1])/3          #c variable solved
        c[l] = (d[l+1]-d[l])/h[l] - h[l]*b[l]/3 - h[l]*b[l+1]/6
        # a[l] = (b[l+1] - b[l])/(3*h[l])                             #a variable solved
        a[l] = (d[l+1] - d[l])/(6*h[l])

    #Remove last element to have 80 intervals
    b = b[:-1]
    d = d[:-1]
    
    #store the solutions to Sol in the proper place
    Coeff[0,:] = a
    Coeff[1,:] = b
    Coeff[2,:] = c
    Coeff[3,:] = d         
        
    return Coeff

def funct_cub_spline(domain, nodes_co, coeff):
    ''' Input:  domain = array containing all x-values at which the function will be evaluated
                nodes_co = array containing all the nodes that were used to setup the cubic spline
                coeff = the coefficient of the cubic splines, obtained from the interpolation function
        Output: q = the function values at all predefined x values (domain)
    '''
    q = []
    for co in domain:
        if co == nodes_co[0]:
            index = 0
        if co == nodes_co[-1]:
            index = len(coeff[0]) -1
        else:
            index = np.where(np.abs(nodes_co) > abs(co))[0][0] -1

        val = coeff[0][index] * (co - nodes_co[index])**3 + coeff[1][index] * (co - nodes_co[index])**2 + coeff[2][index] * (co - nodes_co[index]) + coeff[3][index] 
        q.append(val)

    q = np.asarray(q)

    return q
