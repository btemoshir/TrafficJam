# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 17:14:56 2018

@author: Larry
"""
import numpy as np
import matplotlib.pyplot as plt

u = .2

N = 100.
dx = 1/N
dt = 1/N
dk = 2*np.pi/N

def d0(k):
    return np.sin(k*dx)/(k*dx)

def lw(k):
    a = d0(k)
    b = 8j*(1-2*u)*np.sin(k*dx/2.)**2/(k*dx/2.)
    return np.abs(a + b)

k = np.arange(0, 2*np.pi, dk)

plt.plot(k, d0(k), label="d0")
plt.plot(k, (lw(k)), label="LW")
plt.legend()
plt.show()



#
#p0 =  0.2 + 0.1*np.sin(np.arange(0, 2*np.pi, 2*np.pi*dx))
#
#k0 = np.fft.fft(p0)
#plt.plot(k, k0)