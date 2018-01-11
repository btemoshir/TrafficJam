# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 15:22:10 2018

@author: Larry
"""

from Schemes import *

L_factor = 100.

schemes = [LF, LW, LF_cons, LW_cons, Godunov]
names = ["LF", "LW", "LF cons", "LW cons", "Godunov", "Lagrange"]

p_vec = [np.copy(p0) for p in range(len(schemes))]
x = np.arange(0, 1, dx/L_factor) 
p_vec.append(x)
pL = 0.2 + 0.1*np.sin(np.arange(0, 2*np.pi, 2*np.pi*dx/L_factor))
x0 = np.arange(0, 1, dx)

fig, ax = plt.subplots()


def L2(A, B):
    return np.sum(np.sqrt(((A - B)/A)**2))

def L_inf(A, B):
    return np.max((np.sqrt(((A - B)/A)**2)))

def check_collapse(x):
    test = (x - np.roll(x, -1)) > 0
    count= 0
    for value in test:
        if value:
            count += 1
        if count == 2:
            return True
    return False
    

def mom(p):
    return np.sum((p + np.roll(p, -1))/(2*p))
    return np.sum(np.sum(p)*(p + np.roll(p, -1))/2.)

time = []
def momentum(T):
    t_collapse = 0
    n = len(schemes)
    momenta = [[] for i in range(n+1)]
    
    

    for t in np.arange(0, T, dt):
        p_vec[n] = Lagrange(p_vec[n], pL)
        lagrange_p = np.interp(x0, p_vec[-1], pL, period=1)


        
        for i in range(n):
            p_vec[i] = schemes[i](p_vec[i])
            momenta[i].append(mom(p_vec[i]))
        
        momenta[n].append(mom(lagrange_p))
        
        time.append(t)
        if check_collapse(p_vec[-1]) and t_collapse == 0:
            #break
            t_collapse = t
        
    plt.title("Momenta")
    for i in range(n+1):
        plt.plot(time, momenta[i], label=names[i])
    if t_collapse != 0:
        plt.plot([t_collapse, t_collapse], [np.min(momenta), np.max(momenta)], '--')
    plt.legend()
    plt.xlabel("Time")
    plt.ylabel("Relative error")
    plt.show()
    
    print (np.shape(momenta))

momentum(.9)