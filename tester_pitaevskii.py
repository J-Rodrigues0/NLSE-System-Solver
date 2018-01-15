'''

Tests the solver for a single Schrodinger equation of the type H psi = E psi with H = P²/2 + V

'''

from EQUATION import *
from SYSTEM import *
from TERM import *

import scipy as sc
import numpy as np
from scipy import special

def V(x,y,V0,R0):
	V_ = np.zeros((len(x),len(y)))
	index = np.where(x**2 + y**2 < R0**2)
	V_[index] = V0
	return V_

def psi_0(X,Y):
	return np.exp(-(((X-0.2)**2 + (Y-0.2)**2)/0.05**2))*np.exp(-1.j*10000.*X)


#Initializing system and equation
system = SYSTEM()
eq1 = EQUATION('eq1')

system.add_equation(eq1) #Adding equation to the system

Z = 6
N = 2**Z

L = 1.
dh = 2.*L/(N-1)

X, Y = np.mgrid[-L:L+dh:dh,-L:L+dh:dh]

k_line = np.fft.fftfreq(N,dh)
k_x, k_y = np.meshgrid(k_line,k_line)
K2 = k_x**2 + k_y**2

V0 = 1e10
R0 = 0.8

dt = 1e-4
N_stride = 100
N_steps = 100

term1 = TERM(0.5*K2,'Momentum','P_squared')
eq1.add_term(term1)

term2 = TERM(V(X,Y,V0,R0),'Position','Binding Potential')
eq1.add_term(term2)

f = lambda eq: np.abs(eq.solution)**2
kwargs = {'Function': f,'Variables':{'eq':eq1}}

term3 = TERM(eq1.solution,'Position','Pitaevskii',True,**kwargs)
eq1.add_term(term3)

eq1.solution = psi_0(X,Y,K2)
system.solve(X,Y,dt,N_stride,N_steps)