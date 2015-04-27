import numpy as np
import numba as nb
from numpy import ma
from scipy import constants as spc
from numpy import fft

# Constants
hbar = spc.hbar
c = spc.c
e = spc.e
epsilon_0 = spc.epsilon_0
mu_0 = spc.mu_0
m_e = spc.m_e
pi = np.pi

w_a = 2*pi*6e14 # centre frequency
tau_32 = 1e-13 # lifetime of level 3
tau_21 = 1e-9 # lifetime of level 2
tau_10 = 1e-11 # lifetime of level 1
T_2 = 2.18e-14 # mean time between dephasing events, collision time
N_t = 5.5*spc.N_A # total number of excitable atoms.
E_max = 1 # Maximum E field. Guess. 

# Derived constants
gamma_r = 1/tau_21
gamma_c = spc.e**2*w_a**2/(spc.m_e*6*pi*spc.epsilon_0*spc.c**3)
dw_a = 1/tau_21 + 2/T_2

# Simulation parameters
L = 5e-6 # approx. length of medium
dx = 1e-9 # space step
T = 7e-13 # time for simulation
dt = dx/c # time step
N = int(T/dt)

# Load pregenerated medium
epsilon = np.load('epsilon.npy')
# gain medium == 1, scattering medium == 0
medium_mask = np.int32(ma.masked_greater(epsilon, epsilon_0).filled(0)/(epsilon_0))
medium = (1-medium_mask)**2 #swaps values so scattering medium (high refractive index) = 1
medium_length = epsilon.shape[0]

# PML sigma array
sigma = 100*np.arange(0,300,1)**3
sig = np.hstack((sigma, np.zeros(medium_length-2*300), sigma[::-1]))


# Storage arrays
E_storage = [np.zeros(medium_length)]
H_storage = [np.zeros(medium_length)]
P_storage = [np.zeros(medium_length)]
N0_storage = [np.ones(medium_length)]
N1_storage = [np.zeros(medium_length)]
N2_storage = [np.zeros(medium_length)]
N3_storage = [np.zeros(medium_length)]

# Initial conditions
E = np.zeros(medium_length)
H = np.zeros(medium_length)
P = np.zeros(medium_length)
P_prev = np.zeros(medium_length)
N0 = np.ones(medium_length)
N1 = np.zeros(medium_length)
N2 = np.zeros(medium_length)
N3 = np.zeros(medium_length)

# Set emission location
location = 350
while medium[location] != 0:
	location += 20

@nb.jit(nopython=True)
def update_H(H, E):
	for position in range(0,H.shape[0]-1, 1):
		H[position] = H[position] + (E[position+1] - E[position]) - dt*sig[position]*H[position]
	E[0] = E[1]
	return H, E

@nb.jit(nopython=True)	
def update_E(H, E, P, P_prev):
	for position in range(1,E.shape[0], 1):
		E[position] = E[position] + epsilon_0/epsilon[position]*(H[position] - H[position-1]) - e/(epsilon[position]*E_max)*(P[position] - P_prev[position]) - dt*sig[position]*E[position]
	H[-1] = H[-2]
	return E, H

@nb.jit()
def update_P(P, P_prev, E, N1, N2):
	temp = P.copy()
	for position in range(0, P.shape[0], 1):
		P[position] = dt**2*(gamma_r*e*E_max*N_t)/(gamma_c*m_e*(1 + dw_a*dt))*(N1[position] - N2[position])*E[position] + (2 - w_a**2*dt**2 + dw_a*dt)/(1 + dw_a*dt)*P[position] - 1/(1 + dw_a*dt)*P_prev[position]
	P_prev = temp
	return P, P_prev

@nb.jit()
def update_N(N0 ,N1, N2, N3, E, P, P_prev, P_r):
	N0_next = np.ones(N0.shape[0])
	N1_next = np.zeros(N0.shape[0])
	N2_next = np.zeros(N0.shape[0])
	N3_next = np.zeros(N0.shape[0])
	for position in range(0, N0.shape[0], 1):
		if medium_mask[position] == 1:
			N3_next[position] = N3[position] + dt*(P_r*N0[position] - N3[position]/tau_32)
			N2_next[position] = N2[position] + dt*(N3[position]/tau_32 + E_max*e/(dt*N_t*hbar*w_a)*E[position]*(P[position] - P_prev[position]) - N2[position]/tau_21)
			N1_next[position] = N1[position] + dt*(N2[position]/tau_21 - E_max*e/(dt*N_t*hbar*w_a)*E[position]*(P[position] - P_prev[position]) - N1[position]/tau_10)
			N0_next[position] = N0[position] + dt*(N1[position]/tau_10 - P_r*N0[position])
	return N0_next, N1_next, N2_next, N3_next


for timestep in range(250000):
	
	if timestep < 400:
		P_r = 1e15
	else:
		P_r = 0
	
	P, P_prev = update_P(P, P_prev, E, N1, N2)
	
	H, E = update_H(H, E)

	E, H = update_E(H, E, P, P_prev)

	N0, N1, N2, N3 = update_N(N0 ,N1, N2, N3, E, P, P_prev, P_r)

	E[location] += np.exp(-(timestep-100)**2/100.)

	if timestep % 1 == 0 and timestep > 125000:
		# store data in storage list
		E_storage.append(E.copy())
		#H_storage.append(H.copy())
		#P_storage.append(P.copy())
		#N0_storage.append(N0.copy())
		#N1_storage.append(N1.copy())
		#N2_storage.append(N2.copy())
		#N3_storage.append(N3.copy())

	print(timestep, end='\r')

#E_storage = np.array(E_storage)
H_storage = np.array(H_storage)
#P_storage = np.array(P_storage)
#N0_storage = np.array(N0_storage)
#N1_storage = np.array(N1_storage)
#N2_storage = np.array(N2_storage)
#N3_storage = np.array(N3_storage)

#np.save('E_storage',E_storage)
#np.save('H_storage',H_storage)

# Spectrum calculations

# Sum electric field in space to get time evolution
Z = np.sum(E_storage, axis=1)

# Take fourier tranform
ftransform = np.absolute(fft.rfft(Z))

# Save data
np.save('Z', Z)
np.save('ftransform', ftransform)