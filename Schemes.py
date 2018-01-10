import numpy as np
import matplotlib.pyplot as plt


N = 100.
dx = 1/N
dt = 1/(N)
T = .9
v0 = 1
rho0 = 1

p0 =  0.2 + 0.1*np.sin(np.arange(0, 2*np.pi, 2*np.pi*dx))
p = np.copy(p0)


def flux(shift, p):    
    """ Helper function for conservative schemes. """
    if shift == 0:
        return v0*(1 - p/rho0)*p
    else:
        return v0*(1 - np.roll(p, -shift)/rho0)*np.roll(p, -shift)

def velocity(shift = 0, array=p):
    """ Helper function for Godunov scheme. """
    if shift == 0:
        return v0*(1 - array/rho0)
    else:
        return v0*(1 - np.roll(array, -shift)/rho0)

def LF(p):
    """ Finite elements Lax-Friedrichs """    
    A = (np.roll(p, -1) + np.roll(p, 1))/2.
    B = -dt*(1 - 2*p)*(np.roll(p, -1) - np.roll(p, 1))/(2*dx)
    return (A + B)

def LW(p):
    """ Finite elements Lax-Wendroff """
    A = -(1 - 2*p)*(np.roll(p, -1) - np.roll(p, 1))/(2*dx)
    B = dt/2*(1 - 2*p)**2*(np.roll(p, -1) - 2*p + np.roll(p, 1))/(dx*dx)
    p += dt*(A + B)
    return p
     
def LF_cons(p):
    """ Lax-Friedrichs conservative scheme. """
    f_left = 0.5*(flux(-1, p) + flux(0, p)) - dt/(2*dx)*(p - np.roll(p, 1))
    f_right = 0.5*(flux(0, p) + flux(1, p)) - dt/(2*dx)*(np.roll(p, -1) - p)

    p += -dt/dx*(f_right - f_left)
    return p

def LW_cons(p):    
    """ Lax-Wendroff conservative scheme. """
    u_right = 0.5*(np.roll(p, -1) + p) - dt/(2*dx)*(flux(1, p) - flux(0, p))
    u_left  = 0.5*(p + np.roll(p, 1))  - dt/(2*dx)*(flux(0, p) - flux(-1, p))

    f_right = v0*(1 - u_right/rho0)*u_right
    f_left = v0*(1 - u_left/rho0)*u_left
    
    p += -dt/dx*(f_right - f_left)
    return p

def Godunov(p):
    """ Godunov scheme """
    vl = velocity(0, p)
    vr = velocity(1, p)
    ul = p
    ur = np.roll(p, -1)

    #Shock or rarefaction?
    Shock = vl - vr > 0 # 1 is shock, 0 = rarefaction.
    Rare = Shock != 1

    #Direction of shock. 
    Sl = (vl + vr)/2. < 0
    Sr = Sl != 1

    #Sign of velocity.
    v_right = (vl >= 0) == (vr > 0)
    v_centre = (vl < 0) == (vr > 0) #Not necessary?
    v_left = (vr <= 0) == (vl > 0)
    
    # Adding shock.
    u  = Shock*Sl*ur   
    u += Shock*Sr*ul

    #Adding rarefaction.
    u += Rare*v_right*ul
    u += Rare*v_left*ur

    f_right = velocity(0, u)*u
    f_left  = velocity(-1, u)*np.roll(u, 1)

    p += -dt/dx*(f_right - f_left)
    return p

def Lagrange(X, p0 = p0):
    """ Lagrangian solution. Returns x, not p!"""
    X = np.mod(X + (1-2*p0)*dt, 1)
    return X

def L2(A, B):
    return np.sum(np.sqrt(((A - B)/A)**2))/N


def dxdt(scheme, lower, upper, dr):
    time = np.arange(0, T, dt)
    
    for ratio in np.arange(lower, upper, dr):
        dx = ratio*dt

        p0 =  0.2 + 0.1*np.sin(np.arange(0, 2*np.pi, 2*np.pi*dx))
        p = np.copy(p0)
        
        x = np.arange(0, 1, dx) 
        l = np.arange(0, 1, dx) 

        difference = []

        for t in np.arange(0, T, dt):
            p = scheme(p)
            l = Lagrange(l, p0)
            lp = np.interp(x, l, p0, period=1)

            difference.append(L2(lp, p))
        
        plt.plot(time, difference, label=str(dx/dt))
    plt.legend()
    plt.xlabel("Time")
    plt.ylabel("Error")
    plt.title("dx dt comparison")
    
    plt.show()

if __name__ == "__main__":
    dxdt(LF, .8, 1.4, .1)






