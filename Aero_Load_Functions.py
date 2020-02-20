import numpy as np 

def aero_load_interpolation(force_or_torque):
    ''' Input:  force_or_torque =   Choose if you want to find the coefficients for the interpolant q(x) ['force'] or the coefficients for the interpolant t(x) ['shear']
        Output: coeff =             These are the coefficients a0, a1, a2, ..., a40 which correspond to q(x) or t(x)
                q_c or t_c =        These are the loads obtained after integrating along the chord for each spanwise location
                x_dom =             This is the x domain over which the function is interpolated
    '''
    
    aero_data = '/Users/frans/Desktop/AE/AE3212_II-Simulation_Varification_and_Validation/Python/Local/aerodynamicloaddo228.dat'
    data = np.genfromtxt(aero_data,delimiter=',') # column is spanwise // row is chordwise

    ''' Input variables '''
    Ca = 0.515  # [m] / chord length 
    la = 2.691  # [m] / span length
    Nb = 41     # number of spanwise segments
    Nc = 81     # number of chordwise segments
    sc_z = -0.1 # assumption for now / chord runs from 0 to -0.515

    ''' Start program '''
    # Constants
    dc = Ca/Nc                          # segment width in chordwise direction
    db = la/Nb                          # segment width in spanwise direction

    x_dom, z_dom = aero_co()
    shear_dist = z_dom - sc_z

    #---------- Make choice between computing   q(x) needed for Mz      --> 'force'
    #                                           t(z) needed for Torque  --> 'torque'
    if force_or_torque == 'force':

        q_c = np.zeros(Nb)                  # empty array that will contain integrated q values along the chord for each span segment location
        for bi in range(Nb):                # repeating process for each spanwise location
            q_cs = data[:,bi]               # array containing all point loads along the chord for a spanwise location
            # Integrating q(z)
            q_c[bi] = simpsons_rule_integration(z_dom,Nc,q_cs)

        # Interpolation to obtain q(x)
        # Assume as solution,     q(x) = a0 + a1*x + a2*x^2 + a3*x^3 + ... + a40*x^40
        # Coefficients can be found by solving Ax=q_c
        # with A = [ 1   , x0  , x0^2  , x0^3  , ... , x0^40  ]
        #          [ 1   , x1  , x1^2  , x1^3  , ... , x1^40  ]
        #          [ ... ,                           , ...    ]
        #          [ 1   , x40 , x40^2 , x40^3 , ... , x40^40 ]

        A = np.ones((Nb,Nb))
        for column in range(1,Nb):
            A[:,column] = x_dom**column
        coeff = np.linalg.solve(A,q_c) # coefficients of the function // q(z) = a0 + a1*z + a2*z^2 + ... + a60*z^60

        return coeff, q_c, x_dom
    
    #---------
    if force_or_torque == 'torque':

        t_c = np.zeros(Nb)
        for bi in range(Nb):
            t_cs = data[:,bi] * shear_dist
            # Integrating t(z)
            t_c[bi] = simpsons_rule_integration(z_dom, Nc, t_cs)

        # Interpolation to obtain t(x)
        # Assume as solution,     t(x) = a0 + a1*x + a2*x^2 + a3*x^3 + ... + a40*x^40
        # Coefficients can be found by solving Ax=q_c
        # with A = [ 1   , x0  , x0^2  , x0^3  , ... , x0^40  ]
        #          [ 1   , x1  , x1^2  , x1^3  , ... , x1^40  ]
        #          [ ... ,                           , ...    ]
        #          [ 1   , x40 , x40^2 , x40^3 , ... , x40^40 ]

        A = np.ones((Nb,Nb))
        for column in range(1,Nb):
            A[:,column] = x_dom**column
        coeff = np.linalg.solve(A,t_c) # coefficients of the function // q(z) = a0 + a1*z + a2*z^2 + ... + a60*z^60

        return coeff, t_c, x_dom

def aero_load_function(coordinates, coeff, direction):
    ''' input:  coordinates =   x values for q(x) or t(x)  
                coeff =         coefficients from the interpolation
                direction =     choose between 'x' or 'z'
        output: val =           function values for the coordinates used as input
    '''
    
    Nb = 41     # number of spanwise segments
    Nc = 81

    if direction == 'x':
        N = Nb
    if direction == 'z':
        N = Nc

    if np.size(coordinates) >1: # when there is more then 1 coordinate
        val = []
        for co in coordinates:
            vect = np.ones(N)
            for i in range(1,N):
                vect[i] = co**i
            sol = np.dot(vect,coeff)
            val.append(sol)

    else: # in case there is only 1 coordinate
        vect = np.ones(N)
        for i in range(1,N):
            vect[i] = coordinates**i
        val = np.dot(vect,coeff)

    return val

def aero_load_1_int_function(coordinates, coeff, direction):
    ''' input:  coordinates =   x values  
                coeff =         coefficients from the interpolation
                direction =     choose between 'x' or 'z'
        output: val =           function values for the coordinates used as input
    '''
    
    Nb = 41     # number of spanwise segments
    Nc = 81

    if direction == 'x':
        N = Nb
    if direction == 'z':
        N = Nc

    if np.size(coordinates) >1: # when there is more then 1 coordinate
        val = []
        coeff1 = np.zeros(len(coeff))
        for co in coordinates:
            vect = np.zeros(N)
            for i in range(N):
                vect[i]  = co**(i+1)
                coeff1[i] = coeff[i]/(i+1)
            sol = np.dot(vect,coeff1)
            val.append(sol)

    else: # in case there is only 1 coordinate
        vect = np.zeros(N)
        coeff1 = np.zeros(len(coeff))
        for i in range(N):
            vect[i]  = coordinates**(i+1)
            coeff1[i] = coeff[i]/(i+1)
        val = np.dot(vect,coeff1)

    return val

def aero_load_2_int_function(coordinates, coeff, direction):
    ''' input:  coordinates =   x values  
                coeff =         coefficients from the interpolation
                direction =     choose between 'x' or 'z'
        output: val =           function values for the coordinates used as input
    '''
    
    Nb = 41     # number of spanwise segments
    Nc = 81

    if direction == 'x':
        N = Nb
    if direction == 'z':
        N = Nc

    if np.size(coordinates) >1: # when there is more then 1 coordinate
        val = []
        coeff1 = np.zeros(len(coeff))
        coeff2 = np.zeros(len(coeff))
        for co in coordinates:
            vect = np.zeros(N)
            for i in range(N):
                vect[i]   = co**(i+2)
                coeff1[i] = coeff[i]/(i+1)
                coeff2[i] = coeff1[i]/(i+2)
            sol = np.dot(vect,coeff2)
            val.append(sol)

    else: # in case there is only 1 coordinate
        vect = np.zeros(N)
        coeff1 = np.zeros(len(coeff))
        coeff2 = np.zeros(len(coeff))
        for i in range(N):
            vect[i]   = coordinates**(i+2)
            coeff1[i] = coeff[i]/(i+1)
            coeff2[i] = coeff[i]/(i+2)
        val = np.dot(vect,coeff2)

    return val

def aero_load_3_int_function(coordinates, coeff, direction):
    ''' input:  coordinates =   x values  
                coeff =         coefficients from the interpolation
                direction =     choose between 'x' or 'z'
        output: val =           function values for the coordinates used as input
    '''
    
    Nb = 41     # number of spanwise segments
    Nc = 81

    if direction == 'x':
        N = Nb
    if direction == 'z':
        N = Nc

    if np.size(coordinates) >1: # when there is more then 1 coordinate
        val = []
        coeff1 = np.zeros(len(coeff))
        coeff2 = np.zeros(len(coeff))
        coeff3 = np.zeros(len(coeff))
        for co in coordinates:
            vect = np.zeros(N)
            for i in range(N):
                vect[i]   = co**(i+2)
                coeff1[i] = coeff[i]/(i+1)
                coeff2[i] = coeff1[i]/(i+2)
                coeff3[i] = coeff2[i]/(i+3)
            sol = np.dot(vect,coeff3)
            val.append(sol)

    else: # in case there is only 1 coordinate
        vect = np.zeros(N)
        coeff1 = np.zeros(len(coeff))
        coeff2 = np.zeros(len(coeff))
        coeff3 = np.zeros(len(coeff))
        for i in range(N):
            vect[i]   = coordinates**(i+2)
            coeff1[i] = coeff[i]/(i+1)
            coeff2[i] = coeff[i]/(i+2)
            coeff3[i] = coeff[i]/(i+3)
        val = np.dot(vect,coeff3)

    return val

def simpsons_rule_integration(coordinates, n, f):
    ''' Input:      coordinates =   the domain over which to evaluate the integral
                    n =             number of nodes used
                    f =             function values 
        Output:     int =           integrated value

        Remarks:    this method can only be used when n is an odd number
    '''
    if n%2 == 0:
        raise ValueError('n (the number of nodes used) must be an odd number')

    if coordinates[1] < 0:
        a = coordinates[-1]
        b = coordinates[0]
    else:
        a = coordinates[0]
        b = coordinates[-1]

    h = (b-a)/n
    int = h/3 * np.sum(f[0:-1:2] + 4*f[1::2] + f[2::2])
    
    return int

def aero_co():
    #x and z discretization
    Nz = 81
    Nx = 41
    la = 2.691
    Ca = 0.515

    xdisc = np.zeros(Nx)
    #determine the x locations of the distcretization of q(x,z)
    for i in range(1,Nx):
        #calculate theta_(x,i)
        thetaxi = ((i-1)/Nx)*np.pi
        #calculate theta_(x,i+1)
        thetaxi2 = ((i)/Nx)*np.pi
        #calculate x_i
        xi = 0.5 * ((la/2)*(1-np.cos(thetaxi))+(la/2)*(1-np.cos(thetaxi2)))
        xdisc[i] = xi
    xdisc = xdisc
    
    zdisc = np.zeros(Nz)
    #determine the y locations of the distcretization of q(x,z)
    for i in range(1,Nz):
        #calculate theta_(z,i)
        thetazi = ((i-1)/Nz)*np.pi
        #calculate theta_(z,i+1)
        thetazi2 = ((i)/Nz)*np.pi
        #calculate z_i
        zi = 0.5 * ((Ca/2)*(1-np.cos(thetazi))+(Ca/2)*(1-np.cos(thetazi2)))
        zdisc[i] = zi
    zdisc = -1 * zdisc

    return xdisc, zdisc
