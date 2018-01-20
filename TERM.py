'''

Defines the class TERM which is a term in an EQUATION object in the system of NLSE equations we want to solve.

The TERM has an optional name, a defining matrix and the representation in which it should be solved 
(either positions or momentums).

'''

class TERM:
	def __init__(self,matrix,representation,name = '',time_variant = False,time_derivative=False,auxiliary = False,**kwargs):
		self.matrix = matrix #Matrix representation of the term
		self.representation = representation #Representation of the term. Either 'Momentum' or 'Position'
		self.name = name #Name of the term (optional)
		self.time_variant = time_variant #If the term if time variant
		self.time_derivative = time_derivative #If the term is to be integrated in time
		self.auxiliary = auxiliary #If the term is auxiliary it will still be integrated in time but it wont 
								   #considered in the equation

		if time_variant or time_derivative:			
			self.function = kwargs['Function'] #Function to call to update the term if time variant
			self.variables = kwargs['Variables'] #kwargs to pass to self.function to update the term


	def __str__(self): #String conversion
		return '--TERM-- \nName: %s \nRepresentation: %s\n' %(self.name,self.representation)
