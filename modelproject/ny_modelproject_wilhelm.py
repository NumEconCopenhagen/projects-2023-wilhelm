import numpy as np
from scipy import optimize


def k_solow_equation(k, h, alpha, phi, n, g, s_k, s_h, delta):
    ''' find the solow equation for capital
    
    Args: 
        k (float): capital
        h (float): human capital
        alpha (float): capital share
        phi (float): human capital share
        n (float): population growth rate
        g (float): rate of technological development
        delta (float): depreciation rate
    
    Returns: 
    
        change in technology adjusted captal pr capita
        
    '''
    return (s_k*(k**alpha)*h**(phi)-(n+g+delta+n*g)*k)/((1+n)*(1+g))


def h_solow_equation(k, h, alpha, phi, n, g, s_k, s_h, delta):
    ''' find the solow equation for capital
    
    Args: 
        k (float): capital
        h (float): human capital
        alpha (float): capital share
        phi (float): human capital share
        n (float): population growth rate
        g (float): rate of technological development
        delta (float): depreciation rate
    
    Returns: 
    
        change in technology adjusted human captal pr capita
        
    '''
    return (s_h*(k**alpha)*(h**phi)-(n+g+delta+n*g)*h)/((1+n)*(1+g))


def solve(alpha, phi, n, g, s_k, s_h, delta):
    ''' find the solow equation for capital
    
    Args: 
        alpha (float): capital share
        phi (float): human capital share
        n (float): population growth rate
        g (float): rate of technological development
        delta (float): depreciation rate
    
    Returns: 
    
        nullclines for capital and human capital
        
    '''
    # a. create grids
    k_vec = np.linspace(1.0e-10, 500, 1000)
    h_vec_k= np.empty(1000)
    h_vec_h  = np.empty(1000)
    
    # b. loop over every k and solve the solow equations
    for i, k in enumerate(k_vec):
        
        # i. obj function
        def obj(k, x, alpha, phi, n, g, s_k, s_h, delta): 
            return h_solow_equation(x, k, alpha, phi, n, g, s_k, s_h, delta)
        
        # ii. find root
        result = optimize.root_scalar(obj, args=(k, alpha, phi, n, g, s_k, s_h, delta), method='brentq', bracket=[1.0e-10, 500])
 
        h_vec_k[i] = result.root

        # solve for delta_k=0
        
        # i. obj function
        def obj_1(x, h, alpha, phi, n, g, s_k, s_h, delta):
            return k_solow_equation(x, h, alpha, phi, n, g, s_k, s_h, delta)

        # ii. find root
        result = optimize.root_scalar(obj, args=(k, alpha, phi, n, g, s_k, s_h, delta), method='brentq', bracket=[1.0e-10, 500])
        
        h_vec_k[i] = result.root

    return k_vec, h_vec_k, h_vec_h