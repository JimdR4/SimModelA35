import numpy as np 
from matplotlib import pyplot as plt

def aero_load_interpolation(direction):
    ''' Input:  direction:  choose between 'x' or 'x_torque' 
        Output: 1. coeff:   coefficients of polynomial; 41 coefficients
                2. q_x or q_x_torque:   vector of all point forces in spanwise direction at constant chord segment location 
                                        OR vector of all point torques in spanwise direction at constant chord segment location
                3. x_co_q: x coordinate where point loads act, which is in the center of each segment 
    '''

    aero_data = '/Users/frans/Desktop/AE/AE3212_II-Simulation_Varification_and_Validation/Python/Aero_load/aerodynamicloaddo228.dat'
    data = np.genfromtxt(aero_data,delimiter=',') # column is spanwise // row is chordwise
    # b stands for span
    # c stands for chord
    # q stands for distributed point force 

    Ca = 0.515  # [m] / chord length 
    la = 2.691  # [m] / span length
    Nb = 41     # number of spanwise segments
    Nc = 81     # number of chordwise segments

    dc = Ca/Nc  # segment width in chordwise direction
    db = la/Nb  # segment width in spanwise direction

    sc_z = -0.2 # assumption for now / chord runs from 0 to -0.515

    if direction == 'x':
        # integrate along chordline for all spanwise segment
        q_x = np.zeros(Nb)

        for b_i in range(Nb):
            q_along_c = data[:,b_i]         # vector of all point forces in chordwise direction at same spanwise segment location 
            q_x[b_i]  = sum(q_along_c * dc) # numerical integration

        x_co_q = np.arange(db/2,la,db)      # x coordinate where point loads act, which is in the center of each segment

        # *** interpolation ***
        #    Assume as solution,     w(x) = a0 + a1*x + a2*x^2 + a3*x^3 + ... + a40*x^40
        #    Coefficients can be found by solving,   A /cdot coeff = q_x
        #           with A = [ 1   , x0  , x0^2  , x0^3  , ... , x0^40  ]
        #                    [ 1   , x1  , x1^2  , x1^3  , ... , x1^40  ]
        #                    [ ... ,                           , ...    ]
        #                    [ 1   , x40 , x40^2 , x40^3 , ... , x40^40 ]

        A = np.zeros((Nb,Nb))            
        A[:,0] = 1
        for column in range(1,Nb):
            A[:,column] = x_co_q**column

        coeff = np.linalg.solve(A,q_x)

        return coeff, q_x, x_co_q

    if direction == 'x_torque':
        # integrate along chordline for all spanwise segment
        q_x_torque = np.zeros(Nb)

        z_co_q = -1 * np.arange(dc/2,Ca,dc)      # negative z coordinate where point loads act, which is in the center of each segment along the chord

        for b_i in range(Nb):
            q_along_c  = data[:,b_i]                      # vector of all point forces in chordwise direction at same spanwise segment location 
            shear_dist = z_co_q - sc_z                    # array with all distances from segment location to shear center
            q_x_torque[b_i]   = sum(q_along_c * dc * shear_dist) # numerical integration of torque along the chord 

        x_co_q = np.arange(db/2,la,db)      # x coordinate where point loads act, which is in the center of each segment

        # *** interpolation ***
        #    Assume as solution,     w(x) = a0 + a1*x + a2*x^2 + a3*x^3 + ... + a40*x^40
        #    Coefficients can be found by solving,   A /cdot coeff = q_x
        #           with A = [ 1   , x0  , x0^2  , x0^3  , ... , x0^40  ]
        #                    [ 1   , x1  , x1^2  , x1^3  , ... , x1^40  ]
        #                    [ ... ,                           , ...    ]
        #                    [ 1   , x40 , x40^2 , x40^3 , ... , x40^40 ]
                
        A = np.zeros((Nb,Nb))            
        A[:,0] = 1
        for column in range(1,Nb):
            A[:,column] = x_co_q**column

        coeff = np.linalg.solve(A,q_x_torque)

        return coeff, q_x_torque, x_co_q

def w(coordinates,coeff,direction):
    ''' input:  coordinates: values used in the interpolated function
                coeff: coefficients from the interpolation
                direction: choose between 'x' or 'z'
        output: sol: function values for the coordinates used as input
    '''
    
    Nb = 41     # number of spanwise segments

    if direction == 'x':
        N = Nb
    if direction == 'x_torque':
        N = Nb

    if np.size(coordinates) >1: # when there is more then 1 coordinate
        sol = []
        for co in coordinates:
            vect = np.zeros(N)
            vect[0] = 1
            for i in range(1,N):
                vect[i] = co**i
        
            sol.append( np.dot(vect,coeff) )

    else: # in case there is only 1 coordinate
        vect = np.zeros(N)
        vect[0] = 1
        for i in range(1,N):
            vect[i] = coordinates**i
        
        sol = np.dot(vect,coeff)

    return sol

direction = 'x_torque' # choose between 'x' OR 'x_torque'
coeff, q, co = aero_load_interpolation(direction)

x = co
y = w(x, coeff, direction)
plt.plot(x,y)    # interpolation plot
plt.scatter(x,q) # scatter plot from the data file
plt.grid()
plt.show()

print()