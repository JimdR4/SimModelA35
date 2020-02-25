'''
Created on Thu Feb 20 2020
@author: Frans Vossen
'''

import numpy as np
from CubSplineInterpolation import interpolation
# from matplotlib import pyplot as plt

def simp_val(a, b , n, values):
    ''' Input:      a = the start value of the integration
                    b = the end value of the integration
                    n = number of nodes used
                    f = function values 
        Output:     int = integrated value

        Remarks:    this method can only be used when n is an odd number
                    this method can only be used for a constant node spacing
    '''
    if n%2 == 0:
        raise ValueError('n (the number of nodes used) must be an odd number')

    h = (b-a)/n
    int = h/3 * np.sum(values[0:-1:2] + 4*values[1::2] + values[2::2])
    
    return int

def simp_func(x, n, nodes_co, values):
    ''' 
    '''
    if n%2 == 0:
        raise ValueError('n (the number of nodes used) must be an odd number')

    if nodes_co[1] < 0:
        a = nodes_co[-1]
        b = nodes_co[0]
    else:
        a = nodes_co[0]
        b = nodes_co[-1]

    h = (b-a)/n

    value = [values[0]]
    node_co = [nodes_co[0]]
    val = values[0]
    
    for i in range(1,len(values)-1,2):
        val = val + h/3 * (values[i-1] + 4*values[i] + values[i+1])
        value.append(val)
        node_co.append(nodes_co[i])

    tot_int = simp_val(nodes_co[0], nodes_co[-1], n, values)
    value.append(tot_int)
    node_co.append(nodes_co[-1])

    value = np.asarray(value)
    node_co = np.asarray(node_co)

    int = interpolation(x, value, node_co)

    return int
    
# error = []
# h = []
# for i in range(8001, 13002, 1000):
#     n = i
#     x_0 = 0
#     x_n = 10

#     x = np.linspace(x_0,x_n,n)
#     y = 3*x**4 + 5*x**3 + 2*x**2 + 8*x + 6
#     int_exa = (3/5*x**5 + 5/4*x**4 + 2/3*x**3+4*x**2+6*x)[-1]

#     int_num_fun, dom = funct_simpsons_rule_integration(x, n, y) # numerical values van de integraal
#     coeff = interpolation(int_num_fun, dom) # interpolation van de numerical values
#     val = funct_cub_spline(x, dom, coeff)[-1]   # function values van the splines

#     h.append( (x_n-x_0)/n )
#     error.append( abs(int_exa - val) )
#     print(i)

# plt.plot(h,error,'-c')
# plt.show()


# print()