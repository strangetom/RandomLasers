import numpy as np
import numba as nb
from numpy import ma
from scipy import constants as spc
from numpy import fft

hbar = spc.hbar
c = spc.c
e = spc.e
epsilon_0 = spc.epsilon_0
epsilon = epsilon_0
mu_0 = spc.mu_0
m_e = spc.m_e

L = 5e-6 # approx. length of medium
dx = 1e-9 # space step
T = 7e-13 # time for simulation
dt = dx/c # time step
N = int(T/dt)

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
# gain medium == 1, scattering medium == 0
medium_mask = np.int32(ma.masked_greater(epsilon, epsilon_0).filled(0)/(epsilon_0))
medium = (1-medium_mask)**2 #swaps values so scattering medium (high refractive index) = 1

medium_length = epsilon.shape[0]


# Storage arrays
E_storage = [np.zeros(medium_length)]
H_storage = [np.zeros(medium_length)]
sigma = 100*np.arange(0,300,1)**3
sig = np.hstack((sigma, np.zeros(medium_length-2*300), sigma[::-1]))

# Initial conditions
E = np.zeros(medium_length)
H = np.zeros(medium_length)

# Set emission location
location = 350
while medium[location] != 0:
	location += 20

@nb.jit(nopython=True)
def update_H(H, E):
	for position in range(0,H.shape[0]-1,1):
		H[position] = H[position] + (E[position+1] - E[position]) - dt*sig[position]*H[position]
	return H

@nb.jit(nopython=True)	
def update_E(H, E):
	for position in range(1,E.shape[0],1):
		E[position] = E[position] + epsilon_0/epsilon[position]*(H[position] - H[position-1]) - dt*sig[position]*E[position]
	return E

for timestep in range(250000):
	
	H = update_H(H, E)
	E[0] = E[1]

	E = update_E(H, E)
	H[-1] = H[-2]

	E[location] += 1000*np.exp(-(timestep-100)**2/100.)

	if timestep % 50 == 0 and timestep > 125000:
		# store data in storage list
		E_storage.append(E.copy())
		H_storage.append(H.copy())

	print(timestep, end='\r')

E_storage = np.array(E_storage)
H_storage = np.array(H_storage)

np.save('E_storage',E_storage)
np.save('H_storage',H_storage)
np.save('Medium', medium)


# Spectrum calculations

# Sum electric field in space to get time evolution
Z = np.sum(E_storage, axis=1)

# Take fourier tranform
ftransform = np.absolute(fft.rfft(Z))

# Calculate frequency limits of transform
mintime = 125000*dt
maxtime = 250000*dt
freqs = np.linspace(1/maxtime, 1/mintime, ftransform.shape[0])

np.save('ftransform', ftransform)
np.save('ftransformfreqs', freqs)




	