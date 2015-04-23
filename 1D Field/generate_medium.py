import numpy as np
from scipy import constants as spc

# Constants
epsilon_0 = spc.epsilon_0

# Simulation parameters
L = 5e-6 # approx. length of medium
dx = 1e-9 # space step

# Medium 
epsilon = [epsilon_0]*50 + [4*epsilon_0]*60 # vector to contain permittivity at each node
W = 1.4 # strength of randomness
b = 60
number_of_cells = 0
while len(epsilon) < int(L/dx):
	a = int(100*(1+W*(np.random.rand()-0.5)))
	new_region = [epsilon_0]*a + [4*epsilon_0]*b
	epsilon.extend(new_region)
	number_of_cells += 1
epsilon.extend([4*epsilon_0]*60 + [epsilon_0]*50)
epsilon = np.array(epsilon)

np.save("epsilon", epsilon)