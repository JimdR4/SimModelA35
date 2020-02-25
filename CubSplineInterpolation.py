import numpy as np

def aero_nodes(Ca, la, Nz, Nx):
    ''' Input:  Ca = chord length
                la = span length
                Nz = number of segments in chordwise direction
                Nx = number of segments in spanwise direction
        Output: Znodes = array containing all z-coordinates of the point loads 
                Xnodes = array containing all x-coordinates of the point loads
    '''

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

def interpolation(x, values, nodes_co):
    ''' Input:  x =         either a value or an array as function variable
                values =    1-D array containing all function values of the data points that needs to be interpolated
                nodes_co =  1-D array containing all grid points         
        Output: val =       function values at the predefined function variable(s) (=input x)

        Remarks:    the 'values' input are the function values at the 'nodes_co' input

                    s_i = a*(x-x_i)**3 + b*(x-x_i)**2 + c*(x-x_i) + d
                    Coeff[0,:] all 'a' coefficients
                    Coeff[1,:] all 'b' coefficients
                    Coeff[2,:] all 'c' coefficients
                    Coeff[3,:] all 'd' coefficients
    '''
    coeff = np.zeros((4, len(values.T)-1))                      # all coefficients of the splines will be stored here

    n = len(nodes_co)-1                                         # number of needed splines 
    h = np.zeros(n)                                             # the stepsize between nodes will be stored in here
    a = np.zeros(n)                                             # a,b,c need to be solved to create our spline
    b = np.zeros(n)
    c = np.zeros(n)
    d = values                                                  # d is always equal to the force
        
    for i in range(n):           
        h[i] = (nodes_co[i+1] - nodes_co[i])                    # stepsize calculation
    
    S = np.zeros((n+1, n+1))                                    # the matrix to store the equations
    S[0,0] = 1                                                  # top right and bottom left corner are  non zero
    S[n,n] = 1
    for j in range(1,n):                                        # storing the values in the matrix S
        S[j,j] = (h[j-1] + h[j])/3
        S[j,j-1] = h[j-1]/6
        S[j,j+1] = h[j]/6
    
    # we first have to solve the matrix equation S*b=v for b
    v = np.zeros(n+1)                                           # array to store values in v
    
    for k in range(len(nodes_co)-2):                            
        v[k] = (d[k+2] - d[k+1])/h[k+1] - (d[k+1] - d[k])/h[k]  # storing the values in the v array
    
    b = np.linalg.solve(S, v)                                   # solving for b
     
    for l in range(n):
        # c[l] = (d[l+1]-d[l])/h[l] - h[l]*(2*b[l]+b[l+1])/3        
        c[l] = (d[l+1]-d[l])/h[l] - h[l]*b[l]/3 - h[l]*b[l+1]/6 # calculating all 'c' coefficients
        # a[l] = (b[l+1] - b[l])/(3*h[l])                      
        a[l] = (d[l+1] - d[l])/(6*h[l])                         # calculating all 'a' coefficients

    b = b[:-1]                                                  # Remove last element to have 80 intervals
    d = d[:-1]
    
    coeff[0,:] = a                                              # store the coefficients in the proper place
    coeff[1,:] = b
    coeff[2,:] = c
    coeff[3,:] = d         
    
    # -------- at this point we have all coefficients for the cubic splines --------

    if np.size(x) > 1:
        val = []
        for co in x:
            if co == nodes_co[0]:
                index = 0
            if co == nodes_co[-1]:
                index = len(coeff[0]) -1
            else:
                index = np.where(np.abs(nodes_co) > abs(co))[0][0] -1

            value = coeff[0][index] * (co - nodes_co[index])**3 + coeff[1][index] * (co - nodes_co[index])**2 + coeff[2][index] * (co - nodes_co[index]) + coeff[3][index] 
            val.append(value)

        val = np.asarray(val)

    else:
        if x == nodes_co[0]:
            index = 0
        if x == nodes_co[-1]:
            index = len(coeff[0]) -1
        else:
            index = np.where(np.abs(nodes_co) > abs(x))[0][0] -1

        val = coeff[0][index] * (x - nodes_co[index])**3 + coeff[1][index] * (x - nodes_co[index])**2 + coeff[2][index] * (x - nodes_co[index]) + coeff[3][index] 
    
    return val
