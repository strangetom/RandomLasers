import numpy as np
from numpy import ma
from scipy import constants as spc

hbar = spc.hbar
c = spc.c
e = spc.e
epsilon_0 = spc.epsilon_0
epsilon = epsilon_0
mu_0 = spc.mu_0
m_e = spc.m_e

L = 50e-6 # approx. length of medium
dx = 1e-9 # space step
T = 7e-13 # time for simulation
dt = dx/c # time step
N = int(T/dt)

# Medium 
epsilon = [epsilon_0]*50 + [4*epsilon_0]*60 # vector to contain permittivity at each node
W = 1.4 # strength of randomness
b = 60
number_of_cells = 0
while len(epsilon) < 1000:
	a = int(100*(1+W*(np.random.rand()-0.5)))
	new_region = [epsilon_0]*a + [4*epsilon_0]*b
	epsilon.extend(new_region)
	number_of_cells += 1
epsilon.extend([4*epsilon_0]*60 + [epsilon_0]*50)
epsilon = np.array(epsilon)
# gain medium == 1, scattering medium == 0
medium_mask = np.int32(ma.masked_greater(epsilon, epsilon_0).filled(0)/(epsilon_0))
medium = (1-medium_mask)**2

medium_length = epsilon.shape[0]


# Storage arrays
E_storage = np.zeros(medium_length)
H_storage = np.zeros(medium_length)
#sig = np.hstack((1e15*np.ones(20), 1e14*np.ones(20), 1e13*np.ones(20), 1e12*np.ones(20), 1e11*np.ones(20), 1e10*np.ones(20), 1e9*np.ones(20), 1e8*np.ones(20), 1e7*np.ones(20), 1e6*np.ones(20), 1e5*np.ones(20), np.zeros(300), 1e5*np.ones(20), 1e6*np.ones(20), 1e7*np.ones(20), 1e8*np.ones(20), 1e9*np.ones(20), 1e10*np.ones(20), 1e11*np.ones(20), 1e12*np.ones(20), 1e13*np.ones(20), 1e14*np.ones(20), 1e15*np.ones(20)))
sigma = 100*np.arange(0,300,1)**3
sig = np.hstack((sigma, np.zeros(medium_length-2*300), sigma[::-1]))

# Initial conditions
E = np.zeros(medium_length)
H = np.zeros(medium_length)

for timestep in range(2000):
	
	for position in range(0,H.shape[0]-1,1):
		H[position] = H[position] + dt/(mu_0*dx)*(E[position+1] - E[position]) - dt*sig[position]*H[position]
	E[0] = E[1]

	for position in range(1,E.shape[0],1):
		E[position] = E[position] + dt/(epsilon[position]*dx)*(H[position] - H[position-1]) - dt*sig[position]*E[position]
	H[-1] = H[-2]

	E[350] += np.exp(-(timestep - 30.)**2 /100.)

	if timestep % 1 == 0:
		# store data in storage list
		E_storage = np.vstack((E_storage,E))
		H_storage = np.vstack((H_storage,H))

	print(timestep, end='\r')

np.savetxt('E_storage.txt',E_storage, delimiter=',')
np.savetxt('H_storage.txt',H_storage, delimiter=',')
np.savetxt('Medium.txt', medium, delimiter=',')




	