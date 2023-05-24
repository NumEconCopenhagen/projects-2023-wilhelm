import numpy as np
from scipy import optimize

def k_solow_equation(k, alpha, n, g, s, delta):
    ''' find the solow equation for capital
    
    Args: 
        k (float): capital
        alpha (float): capital share
        n (float): population growth rate
        g (float): rate of technological development
        s (float): savings rate
        delta (float): depreciation rate
    
    Returns: 
    
        change in technology adjusted captal pr capita
        
    '''

    return (s*(k**alpha)-(n+g+delta+n*g)*k)/((1+n)*(1+g))

def transistion_diagram(k, alpha, n, g, s, delta, t):
    ''' the transistion diagram
    
    Args: 
        k (float): capital
        alpha (float): capital share
        n (float): population growth rate
        g (float): rate of technological development
        s (float): savings rate
        delta (float): depreciation rate
        t (int): periods
    
    Returns: 
    
        diagram whith growth in captial that illustrates the growthpath towards steady state
        
    '''
    # lists for the two lines that make up the diagram
    fortyfive = [0]
    k_values = [k]

    # the 45 line values
    for a in range(1, t):
        x = (n+g+n*g+delta)*a
        fortyfive.append(x)
        
    # capital growth
    for a in range(1, t):
        k_t = s*a**alpha
        k_values.append(k_t)
        
    # make the plot
    plt.figure()
    plt.plot(fortyfive[:t], label=r'$(n+g+ng+\delta)k_t$')
    plt.plot(k_values[:t], label=r'$s\tilde{k}_t^{\alpha}$')
    plt.xlim(0, t)
    plt.ylim(0, fortyfive[-1])
    plt.xlabel('$k_t$')
    plt.ylabel('$k_{t+1}$')
    plt.legend()
    plt.title('Transition diagram')


def k_solow_equation_human(k, h, alpha, phi, n, g, s_k, s_h, delta):
    ''' find the solow equation for capital
    
    Args: 
        k (float): capital
        h (float): human capital
        alpha (float): capital share
        phi (float): human capital share
        n (float): population growth rate
        g (float): rate of technological development
        s_k (float): savings rate for capital
        s_h (float): savings rate for human capital
        delta (float): depreciation rate
    
    Returns: 
    
        change in technology adjusted captal pr capita
        
    '''
    
    return (s_k*(k**alpha)*h**(phi)-(n+g+delta+n*g)*k)/((1+n)*(1+g))


def h_solow_equation_human(k, h, alpha, phi, n, g, s_k, s_h, delta):
    ''' find the solow equation for capital
    
    Args: 
        k (float): capital
        h (float): human capital
        alpha (float): capital share
        phi (float): human capital share
        n (float): population growth rate
        g (float): rate of technological development
        s_k (float): savings rate for capital
        s_h (float): savings rate for human capital
        delta (float): depreciation rate
    
    Returns: 
    
        change in technology adjusted human captal pr capita
        
    '''
    return (s_h*(k**alpha)*(h**phi)-(n+g+delta+n*g)*h)/((1+n)*(1+g))


def solve_human(alpha, phi, n, g, s_k, s_h, delta):
    ''' find the solow equation for capital
    
    Args: 
        alpha (float): capital share
        phi (float): human capital share
        n (float): population growth rate
        g (float): rate of technological development
        s_k (float): savings rate for capital
        s_h (float): savings rate for human capital
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
            return h_solow_equation_human(x, k, alpha, phi, n, g, s_k, s_h, delta)
        
        # ii. find root
        result = optimize.root_scalar(obj, args=(k, alpha, phi, n, g, s_k, s_h, delta), method='brentq', bracket=[1.0e-10, 500])
 
        h_vec_k[i] = result.root

        # solve for delta_k=0
        
        # i. obj function
        def obj_1(x, h, alpha, phi, n, g, s_k, s_h, delta):
            return k_solow_equation_human(x, h, alpha, phi, n, g, s_k, s_h, delta)

        # ii. find root
        result = optimize.root_scalar(obj, args=(k, alpha, phi, n, g, s_k, s_h, delta), method='brentq', bracket=[1.0e-10, 500])
        
        h_vec_k[i] = result.root

    return k_vec, h_vec_k, h_vec_h