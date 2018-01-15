'''

Tests the solver for a single Schrodinger equation of the type H psi = E psi with H = P²/2 + V 
with Absorbing Boundary Conditions (ABCs).

'''

'''

Tests the solver for two Schrodinger equations of the type: P²/2 psi + V psi + abs(psi)² psi + abs(phi)² psi = 0
															P²/2 phi + V phi + abs(phi)² phi + abs(psi)² phi = 0

'''

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
	return np.exp(-(((X+0.4)**2 + (Y-0.)**2)/0.1**2))*np.exp(1.j*10.*X)

Z = 6
N = 2**Z

L = 1.
dh = 2.*L/(N-1)

X, Y = np.mgrid[-L:L+dh:dh,-L:L+dh:dh]

k_line = np.fft.fftfreq(N,dh)
k_x, k_y = np.meshgrid(k_line,k_line)
K2 = k_x**2 + k_y**2

V0 = 1e10
R0 = 0.5

dt = 1e-4
N_stride = 100
N_steps = 100

#Initializing system and equation
system = SYSTEM()
eq1 = EQUATION('eq1')

system.add_equation(eq1) #Adding equation to the system

k_line = np.fft.fftshift(k_line)
k0 = k_line[0]
print(k0)

#Kinetic Energy Term
term1 = TERM(np.sqrt(2)*k0*k_x + k0**2,'Momentum','P_squared')
eq1.add_term(term1)

#Gross-Pitaevskii Term
f = lambda eq: np.abs(eq.solution)**2
kwargs1 = {'Function': f,'Variables':{'eq':eq1}} #(eq1.solution)² for eq1

term3_1 = TERM(eq1.solution,'Position','Gross-Pitaevskii Term',True,**kwargs1)
eq1.add_term(term3_1)


#Initializing solutions
eq1.solution = psi_0(X,Y)

#Solving system
system.solve(X,Y,dt,N_stride,N_steps)