'''

Defines the class EQUATION which will be used as an object in the system of NLSE we want to solve.

'''

class EQUATION:
	def __init__(self):
		self.n_terms = 0
		self.terms = {}

	def add_term(self,term):
		self.terms[self.n_terms] = term  
		self.n_terms += 1

	def solve(self,other):
		
