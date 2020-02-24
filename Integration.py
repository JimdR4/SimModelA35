import numpy as np
from CubSplineInterpolation import interpolation, funct_cub_spline

def val_simpsons_rule_integration(coordinates, n, f):
    ''' Input:      coordinates = the domain over which to evaluate the integral
                    n = number of nodes used
                    f = function values 
        Output:     int = integrated value
        Remarks:    this method can only be used when n is an odd number
                    this method can only be used for a constant node spacing
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

def funct_simpsons_rule_integration(coordinates, n, f):
    if n%2 == 0:
        raise ValueError('n (the number of nodes used) must be an odd number')

    if coordinates[1] < 0:
        a = coordinates[-1]
        b = coordinates[0]
    else:
        a = coordinates[0]
        b = coordinates[-1]

    h = (b-a)/n

    values = [f[0]]
    coordin = [coordinates[0]]
    val = f[0]
    
    for i in range(1,len(f)-1,2):
        val = val + h/3 * (f[i-1] + 4*f[i] + f[i+1])
        values.append(val)
        coordin.append(coordinates[i])

    tot_int = val_simpsons_rule_integration(coordinates, n, f)
    values.append(tot_int)
    coordin.append(coordinates[-1])

    values = np.asarray(values)
    coordin = np.asarray(coordin)

    domain = np.linspace(coordinates[0],coordinates[-1],n)
    coeff = interpolation(values, coordin)
    f_values = funct_cub_spline(domain, coordin, coeff)

    # return values, coordin
    return f_values, domain
    