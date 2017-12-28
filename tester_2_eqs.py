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
	return np.exp(-(((X-0.2)**2 + (Y-0.2)**2)/0.1**2))*np.exp(1.j*10.*X)


def psi_1(X,Y):
	return np.exp(-(((X+0.2)**2 + (Y+0.2)**2)/0.1**2))*np.exp(-1.j*10.*X)

#Initializing system and equation
system = SYSTEM()
eq1 = EQUATION('eq1')
eq2 = EQUATION('eq2')

system.add_equation(eq1) #Adding equation to the system
system.add_equation(eq2) #Adding equation to the system

Z = 6
N = 2**Z

L = 1.
dh = 2.*L/(N-1)

X, Y = np.mgrid[-L:L+dh:dh,-L:L+dh:dh]

k_line = np.fft.fftfreq(N,dh)
k_line = np.fft.fftshift(k_line)
k_x, k_y = np.meshgrid(k_line,k_line)
K2 = k_x**2 + k_y**2

V0 = 1e10
R0 = 0.5

dt = 1e-4
N_stride = 100
N_steps = 100

#Kinetic Energy Term
term1 = TERM(0.5*K2,'Momentum','P_squared')
eq1.add_term(term1)
eq2.add_term(term1)

#Binding Potential Term
term2 = TERM(V(X,Y,V0,R0),'Position','Binding Potential')
eq1.add_term(term2)
eq2.add_term(term2)

#Gross-Pitaevskii Term
f = lambda eq: np.abs(eq.solution)**2
kwargs1 = {'Function': f,'Variables':{'eq':eq1}} #(eq1.solution)² for eq1
kwargs2 = {'Function': f,'Variables':{'eq':eq2}} #(eq2.solution)² for eq2

term3_1 = TERM(eq1.solution,'Position','Gross-Pitaevskii Term',True,**kwargs1)
term3_2 = TERM(eq2.solution,'Position','Gross-Pitaevskii Term',True,**kwargs2)
eq1.add_term(term3_1)
eq2.add_term(term3_2)

#Cross Potential Term
kwargs1 = {'Function': f,'Variables':{'eq':eq2}} #(eq2.solution)² for eq1
kwargs2 = {'Function': f,'Variables':{'eq':eq1}} #(eq1.solution)² for eq2

term4_1 = TERM(eq1.solution,'Position','Cross Term',True,**kwargs1)
term4_2 = TERM(eq2.solution,'Position','Cross Term',True,**kwargs2)
eq1.add_term(term4_1)
eq2.add_term(term4_2)

#Initializing solutions
eq1.solution = psi_0(X,Y)
eq2.solution = psi_1(X,Y)

#Solving system
system.solve(X,Y,dt,N_stride,N_steps)