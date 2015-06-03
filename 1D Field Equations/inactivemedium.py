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

epsilon = np.load('epsilon.npy')
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

	E[location] += np.exp(-(timestep-100)**2/100.)

	if timestep % 1 == 5 and timestep > 200000:
		# store data in storage list
		E_storage.append(E.copy())

	print(timestep, end='\r')

E_storage = np.array(E_storage)

# Spectrum calculations

# Sum electric field in space to get time evolution
Z = np.sum(E_storage, axis=1)

# Take fourier tranform
ftransform = np.absolute(fft.rfft(Z))

# Save data
np.save('Z', Z)
np.save('ftransform', ftransform)




	