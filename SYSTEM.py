'''

Definies the class SYSTEM composed of NLSE of the class EQUATION.
This SYSTEM class has the solvers necessary to solve the system.

'''
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt

class SYSTEM:
	def __init__(self,name=None):
		self.name = name
		self.n_eqs = 0
		self.equations = {}

	def add_equation(self,equation):
		self.equations[self.n_eqs] = equation
		self.n_eqs += 1

	def solve(self,X,Y,dt,N_stride,N_steps):
		t = 0
		for n in range(N_stride):
			print('Stride nº: %s' %n)
			for m in range(N_steps):
				t += dt
				for i in self.equations:
					self.equations[i].parts()
					self.equations[i].step(dt)

			for i in self.equations:
				prob = np.abs(self.equations[i].solution)**2

				#Make the graph and save it
				pl.figure()
				pl.contour(X,Y,self.equations[i].N.real,10)
				pl.contourf(X,Y,prob)
				pl.colorbar()
				plt.savefig('/home/joaorodrigues/Desktop/MCE/Final/figs/%s.png' %(self.equations[i].name+'_'+str(n)),dpi=200)
				pl.close()