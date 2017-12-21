'''

Definies the class SYSTEM composed of NLSE of the class EQUATION.
This SYSTEM class has the solvers necessary to solve the system.

'''

class SYSTEM:
	def __init__(self,name=None):
		self.name = name
		self.n_eqs = 0
		self.equations = {}

	def add_equation(self,equation):
		self.equations[self.n_eqs] = equation