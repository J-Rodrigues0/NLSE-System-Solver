'''

Definies the class SYSTEM composed of NLSE of the class EQUATION.
This SYSTEM class has the solvers necessary to solve the system.

'''
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt

class SYSTEM:
	def __init__(self,name=None):
		self.name = name #Name of the system (optional)
		self.n_eqs = 0 #Number of equations in the system
		self.equations = {} #Dictionary of equations in the system

	def add_equation(self,equation): #Adds an equation to the system
		self.equations[self.n_eqs] = equation
		self.n_eqs += 1

	def solve(self,X,Y,dt,N_stride,N_steps): #Solves the system using the SSSFM
		t = 0
		sum_squared = {}
		for i in self.equations:
			sum_squared[i] = np.sum(np.abs(self.equations[i].solution)**2)

		for n in range(N_stride): #Iterate through the number of strides in the run
			print('Stride nº: %s' %n)
			for m in range(N_steps): #Perform N_steps of the run
				t += dt
				for i in self.equations: #Perform one step of each equation
					self.equations[i].parts() #Update terms (L and N - and, eventually, V)
					self.equations[i].step(dt) #Do one step

					total_prob = np.sum(np.abs(self.equations[i].solution)**2)/sum_squared[i]
					if not 1. - total_prob < 1e-5:
						raise ValueError('Total Probability != 1 : %s' %float(total_prob))

			for i in self.equations: #Build figures
				prob = np.abs(self.equations[i].solution)**2 #Probability density of the solution

				pl.figure()
				if self.equations[i].V != []:
					pl.contour(X,Y,self.equations[i].V.real,10)
				pl.contourf(X,Y,prob)
				pl.colorbar()
				plt.savefig('/home/joaorodrigues/Desktop/MCE/Final/figs/%s.png' %(self.equations[i].name+'_'+str(n)),dpi=200)
				pl.close()