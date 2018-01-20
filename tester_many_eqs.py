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
import time
import pylab as pl

def V(x,y,V0,R0):
	V_ = np.zeros((len(x),len(y))) + V0
	index = np.where(x**2 + y**2 < R0**2)
	V_[index] = 0
	return V_

def psi_0(X,Y):
	return np.exp(-((X**2 + Y**2)/0.1**2))*np.exp(1.j*10000.*X)

Z = 8
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
N_stride = 10
N_steps = 10

N_eqs_vals = np.arange(1,101)
times = []
for N_eqs in N_eqs_vals:
	#Initializing system
	system = SYSTEM()

	eqs = {}

	for i in range(N_eqs):
		eqs[i] = EQUATION('eq%s'%i)

		eq = eqs[i]
		system.add_equation(eq) #Adding equation to the system

		#Kinetic Energy Term
		term1 = TERM(0.5*K2,'Momentum','P_squared')
		eq.add_term(term1)

		#Binding Potential Term
		term2 = TERM(V(X,Y,V0,R0),'Position','Binding Potential')
		eq.add_term(term2)

		#Initializing solutions
		eq.solution = psi_0(X,Y)

	for i in range(N_eqs):
		for j in range(N_eqs):
			#Gross-Pitaevskii and Cross-Potential Terms
			if i == j:
				term_name = 'Gross-Pitaevskii Term'
				g = 1.
			else:
				term_name = 'Cross-Potential from Equation %s' %j
				g = 1e3

			eq = eqs[j]

			f = lambda eq,g: g*np.abs(eq.solution)**2

			kwargs = {'Function': f,'Variables':{'eq':eq,'g':g}} #(eq.solution)² for eq1

			term3 = TERM(eq.solution,'Position',term_name,True,**kwargs)
			eqs[i].add_term(term3)

	#for i in range(N_eqs):
	#	print(system.equations[i])

	#Solving system
	t0 = time.clock()
	system.solve(X,Y,dt,N_stride,N_steps)
	times.append(time.clock()-t0)

pl.figure('times')
pl.plot(N_eqs_vals[:N_eqs],times,'-o')
pl.xlabel('Number of Equations')
pl.ylabel('System Solving Time (s)')
pl.show()