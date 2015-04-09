import numpy as np
from numpy import ma
import numbapro as nb
from scipy import constants as spc

# Define constants
w_a = 2*np.pi*6e14 # centre frequency
tau_32 = 1e-13 # lifetime of level 3
tau_21 = 1e-9 # lifetime of level 2
tau_10 = 1e-11 # lifetime of level 1
T_2 = 2.18e-14 # mean time between dephasing events, collision time
P_r = 1e7 # pumping rate

# Fundamental constants
hbar = spc.hbar
c = spc.c
e = spc.e
epsilon_0 = spc.epsilon_0
mu_0 = spc.mu_0
m_e = spc.m_e

# Derived constants
gamma_r = 1/tau_21
gamma_c = spc.e**2*w_a**2/(spc.m_e*6*np.pi*spc.epsilon_0*spc.c**3)
dw_a = 1/tau_21 + 2/T_2

# Simulation parameters
L = 20e-6 # approx. length of medium
dx = 1e-9 # space step
T = 7e-12 # time for simulation
dt = dx/c # time step
N = int(T/dt)

# Generate medium: 4*epsilon_0 is scattering, epsilon_0 is gain
epsilon = [epsilon_0]*int(300e-9/dx) # vector to contain permittivity at each node
W = 1.4 # strength of randomness
b = int(180e-9/dx)
number_of_cells = 0
while len(epsilon) < L/dx:
	a = int(300e-9*(1+W*(np.random.rand()-0.5))/dx)
	new_region = [epsilon_0]*a + [4*epsilon_0]*b
	epsilon.extend(new_region)
	number_of_cells += 1
epsilon.extend([epsilon_0]*int(300e-9/dx))
epsilon = np.array(epsilon)
# gain medium == 1, scattering medium == 0
medium_mask = np.int32(ma.masked_greater(epsilon, epsilon_0).filled(0)/(epsilon_0))

# Set actual length of medium
I = epsilon.shape[0]
actual_L = I*dx


# Initial conditions
P = np.zeros(I)
P_prev = P.copy()
E = np.zeros(I)
H = np.zeros(I)-1
N0 = 5.5*spc.N_A*np.ones(I)
N1 = np.zeros(I)
N2 = np.zeros(I)
N3 = np.zeros(I)

# Storage arrays
P_storage = []
E_storage = []
H_storage = []
N0_storage = []
N1_storage = []
N2_storage = []
N3_storage = []

# Function definitions
@nb.jit(nopython=True)
def update_P(P, P_prev, E, N1, N2):
	for position in range(1, P.shape[0], 1):
		P[position] = dt**2*gamma_r*e**2/(gamma_c*m_e*(1+dw_a*dt/2))*(N1[position]-N2[position])*E[position] + (2-w_a**2*dt**2)/(1+dw_a*dt/2)*P[position] - (dw_a*dt/2-1)/(1+dw_a*dt/2)*P_prev[position]
	P[0] = P[1]
	return P

@nb.jit(nopython=True)
def update_H(E, H):
	for position in range(1, H.shape[0]-1, 1):
		H[position] = H[position] + dt/(mu_0*dx)*(E[position+1]-E[position])
	H[-1] = H[-2]
	return H

@nb.jit(nopython=True)
def update_E(E, H, P_next, P):
	for position in range(1, E.shape[0], 1):
		E[position] = E[position] + dt/(epsilon[position]*dx)*(H[position]-H[position-1]) + (P[position] - P_next[position])/epsilon[position]
	E[0] = E[1]
	return E

@nb.jit(nopython=True)
def update_N(N0, N1, N2, N3, E, P_next, P, medium_mask):
	for position in range(1, N0.shape[0], 1):
		if medium_mask[position] > 0:
			N3[position] = N3[position] + dt*(P_r*N0[position]-N3[position]/tau_32)
			N2[position] = N2[position] + dt*(N3[position]/tau_32 + E[position]/(hbar*w_a)*( (P_next[position]-P[position])/dt - N2[position]/tau_21))
			N1[position] = N1[position] + dt*(N2[position]/tau_21 - E[position]/(hbar*w_a)*( (P_next[position]-P[position])/dt - N1[position]/tau_10))
			N0[position] = N0[position] + dt*(N1[position]/tau_10 - P_r*N0[position])
	return N0, N1, N2, N3


for timestep in range(300000):

	P_next = update_P(P, P_prev, E, N1, N2)

	H_next = update_H(E, H)
	E_next = update_E(E, H, P_next, P)

	N0_next, N1_next, N2_next, N3_next = update_N(N0, N1, N2, N3, E, P_next, P, medium_mask,)

	P_prev = P
	P = P_next
	E = E_next
	H = H_next
	N0 = N0_next
	N1 = N1_next
	N2 = N2_next
	N3 = N3_next

	if timestep % 100 == 0:
		#store data in storage list
		P_storage.append(P)
		E_storage.append(E)
		H_storage.append(H)
		N0_storage.append(N0)
		N1_storage.append(N1)
		N2_storage.append(N2)
		N3_storage.append(N3)

	print(timestep, end='\r')


P_storage = np.array(P_storage)
E_storage = np.array(E_storage)
H_storage = np.array(H_storage)
N0_storage = np.array(N0_storage)
N1_storage = np.array(N1_storage)
N2_storage = np.array(N2_storage)
N3_storage = np.array(N3_storage)