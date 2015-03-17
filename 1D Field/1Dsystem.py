import numpy as np
from numpy import ma
from scipy import constants as spc

# Define constants
w_a = 2*np.pi*6e14 # centre frequency
tau_32 = 1e-13 # lifetime of level 3
tau_21 = 1e-9 # lifetime of level 2
tau_10 = 1e-11 # lifetime of level 1
T_2 = 2.18e-14 # mean time between dephasing events
P_r = 1e7 # pumping rate

hbar = spc.hbar
c = spc.c
e = spc.e
epsilon_0 = spc.epsilon_0
mu_0 = spc.mu_0
m_e = spc.m_e

gamma_r = 1/tau_21
gamma_c = spc.e**2*w_a**2/(spc.m_e*6*np.pi*spc.epsilon_0*spc.c**3)
dw_a = 1/tau_21 + 2/T_2

L = 50e-6 # approx. length of medium
dx = 1e-9 # space step
#T = 7e-10 # time for simulation
T = 7e-13 # time for simulation
dt = 1e-17 # time step
N = int(T/dt)

# Generate medium
epsilon = [] # vector to contain permittivity at each node
W = 1.4 # strength of randomness
b = int(180e-9/dx)
number_of_cells = 0
while len(epsilon) < L/dx:
	a = int(300e-9*(1+W*np.random.rand()-0.5)/dx)
	new_region = [epsilon_0]*a + [4*epsilon_0]*b
	epsilon.extend(new_region)
	number_of_cells += 1
epsilon = np.array(epsilon)
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

for timestep in range(N):
	# Polarisation
	P_next = dt**2*gamma_r*e**2/(gamma_c*m_e*(1+dw_a*dt))*(N1 - N2)*E + (1+dt*dw_a-w_a**2*dt**2)/(1+dw_a*dt)*P - P_prev/(1+dw_a*dt)

	# Maxwell
	E_next = E + P - P_next + dt/(epsilon*dx)*(H-np.roll(H,-1))
	H_next = H + dt/(mu_0*dx)*(np.roll(E,1)-E)

	# Population
	N3_next = (N3 + dt*(P_r*N0 - N3/tau_32) )#*medium_mask
	N2_next = (N2 + dt*(N3/tau_32 + E/(hbar*w_a)*(P_next-P)/dt - N2/tau_21) )#*medium_mask
	N1_next = (N1 + dt*(N2/tau_21 - E/(hbar*w_a)*(P_next-P)/dt - N1/tau_10) )#*medium_mask
	N0_next = (N0 + dt*(N1/tau_10 - P_r*N0) )#*medium_mask

	P_prev = P
	P = P_next
	E = E_next
	H = H_next
	N0 = N0_next
	N1 = N1_next
	N2 = N2_next
	N3 = N3_next

	if timestep % 100 == 0:
		# store data in storage list
		P_storage.append(P)
		E_storage.append(E)
		H_storage.append(H)
		N0_storage.append(N0)
		N1_storage.append(N1)
		N2_storage.append(N2)
		N3_storage.append(N3)

	print(timestep, end='\r')

"""
P_storage = np.array(P_storage)
E_storage = np.array(E_storage)
H_storage = np.array(H_storage)
N0_storage = np.array(N0_storage)
N1_storage = np.array(N1_storage)
N2_storage = np.array(N2_storage)
N3_storage = np.array(N3_storage)"""