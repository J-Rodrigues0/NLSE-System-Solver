from EQUATION import *
from SYSTEM import *
from TERM import *

import scipy as sc
import numpy as np
from scipy import special

def V(x,y,V0,R0):
	V_ = np.zeros((len(x),len(y))) + V0
	index = np.where(x**2 + y**2 < R0**2)
	V_[index] = 0
	return V_

def psi_0(X,Y):
	return np.exp(-(((X-0.2)**2 + (Y-0.1)**2)/0.1**2)**4)*np.exp(1.j*100.*X) +\
			np.exp(-(((X+0.2)**2 + (Y+0.1)**2)/0.1**2)**4)*np.exp(-1.j*100.*X)

#Initializing system and equation
system = SYSTEM()
eq1 = EQUATION()

system.add_equation(eq1) #Adding equation to the system

Z = 8
N = 2**Z

L = 1.
dh = 2.*L/(N-1)

X, Y = np.mgrid[-L:L+dh:dh,-L:L+dh:dh]

k_line = np.fft.fftfreq(N,dh)
k_x, k_y = np.meshgrid(k_line,k_line)
K2 = k_x**2 + k_y**2

V0 = 1e10
R0 = 0.4

dt = 1e-4
N_stride = 100
N_steps = 100

term1 = TERM(0.5*K2,'Momentum','P_squared')
eq1.add_term(term1)

term2 = TERM(V(X,Y,V0,R0),'Position','Potential')
eq1.add_term(term2)

eq1.solution = psi_0(X,Y)
eq1.solve(X,Y,dt,N_stride,N_steps)