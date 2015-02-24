import numpy as np
from scipy import constants as spc

#define constants
sig_abs = 3.0e-24 # absorption cross section
sig_em = 3.0e-23 # emission cross section
tau_e = 3.2e-6 # lifetime of excited state
N_t = 1.6e25 # total concentration of laser particles
l = 100.e-6 # transport mean free path
kappa_e = 5.e-4 # extinction coefficient (see spec sheet in shared folder)
tau_G = 14.e-9 # pump pulse FWHM
tau_R = 14.e-9 # probe pulse FWHM
t_G = 15.e-9 # time of maxima of pump pulse 
t_R = 15.e-9 # time of maxima of probe pulse 
illum_area = np.pi*(2.e-3)**2 # illumation area
n = 1.35 # average refractive index of medium
c = spc.c/n # speed of light in medium
v = spc.c/n # transport velocity
D = v*l/3. # diffusion coeffecient

E_G = 6.63e-34*c/532e-9 # energy of pump photons
E_A = 6.63e-34*c/700e-9 # energy of emitted photons
I_G0 = 4.e10 # average pump intensity

#define space parameters
L = 0.001 # length of medium
dz = l/2 # space increment
z = np.arange(-l/2, L+l/2+dz, dz) # vector in space
z[0] = z[1] # modifiy space vector so pulses don't decay before entering medium
z[-1] = z[-2] # modifiy space vector so pulses don't decay after exiting medium
J = z.shape[0] # number of space steps

#define time parameters
N = 100000 # number of time steps
T = 50.e-9 # length of time
dt = T/N # time increment

# initial conditions
W_G = np.zeros(z.shape[0])
W_A = np.zeros(z.shape[0])
N_pop = np.zeros(z.shape[0])

# lists for storage
W_G_storage = []
W_A_storage = []
N_pop_storage = []
I_G_storage = []

# function definitions
def ASE_RHS(W_A, N_mid_step):
	"""Calculates the right hand side of the PDE for amplified spontaneous emission"""
	M = np.diagflat([0]+[ dt*D/(2*dz**2) for i in range(J-2)], -1) + np.diagflat([0]+[ (1. - 2*dt*D/(2*dz**2) + dt*sig_em*v*N_mid_step[i]*N_t) for i in range(J-2)]+[0]) + np.diagflat([ dt*D/(2*dz**2) for i in range(J-2)]+[0], 1)
	return M.dot(W_A) + dt*c*N_t*E_A/tau_e/I_G0 * N_pop

def ASE_LHS(N_mid_step):
	"""Calculates the left hand side of the PDE for amplified spontaneous emission"""
	M = np.diagflat([0]+[ -dt*D/(2*dz**2) for i in range(J-2)], -1) + np.diagflat([1]+[ (1. + 2*dt*D/(2*dz**2) - dt*sig_em*v*N_mid_step[i]*N_t) for i in range(J-2)]+[1]) + np.diagflat([ -dt*D/(2*dz**2) for i in range(J-2)]+[0], 1)
	return M

def PUMP_RHS(W_G, N_pop, I_G):
	"""Calculates the right hand side of the PDE for the pump"""
	M = np.diagflat([0]+[ dt*D/(2*dz**2) for i in range(J-2)], -1) + np.diagflat([0]+[(1. - 2*dt*D/(2*dz**2) - dt*sig_abs*v*N_t*(1-N_mid_step[i]) ) for i in range(J-2) ]+[0]) + np.diagflat([ dt*D/(2*dz**2) for i in range(J-2)]+[0], 1)
	return M.dot(W_G) + dt/tau_e *I_G

def PUMP_LHS(N_mid_step):
	"""Calculates the left hand side of the PDE for the pump"""
	M = np.diagflat([0]+[ -dt*D/(2*dz**2) for i in range(J-2)], -1) + np.diagflat([1]+[ (1. + 2*dt*D/(2*dz**2) + dt*sig_abs*v*N_t*(1-N_mid_step[i])) for i in range(J-2)]+[1]) + np.diagflat([ -dt*D/(2*dz**2) for i in range(J-2)]+[0], 1)
	return M

def PROBE_RHS(W_R, N_pop, I_R):
	"""Calculates the right hand side of the PDE for the pump"""
	M = np.diagflat([0]+[ dt*D/(2*dz**2) for i in range(J-2)], -1) + np.diagflat([0]+[(1. - 2*dt*D/(2*dz**2) + dt*sig_em*v*N_t*N_mid_step[i]) for i in range(J-2) ]+[0]) + np.diagflat([ dt*D/(2*dz**2) for i in range(J-2)]+[0], 1)
	return M.dot(W_R) + dt/tau_e *I_R

def PROBE_LHS(N_mid_step):
	"""Calculates the left hand side of the PDE for the pump"""
	M = np.diagflat([0]+[ -dt*D/(2*dz**2) for i in range(J-2)], -1) + np.diagflat([1]+[ (1. + 2*dt*D/(2*dz**2) - dt*sig_em*v*N_t*N_mid_step[i]) for i in range(J-2)]+[1]) + np.diagflat([ -dt*D/(2*dz**2) for i in range(J-2)]+[0], 1)
	return M

def POP_RHS(N_pop, W_G, W_A):
	"""Calculates the right hand side of the PDE for the excited population"""
	term_1 = dt*sig_abs*v*(1-N_pop)*I_G0/(c*E_G) * W_G
	term_2 = dt*sig_em*v*N_pop*I_G0/(c*E_A) * W_A
	term_3 = dt/tau_e * N_pop
	return N_pop + (term_1 - term_2 - term_3)

def I_G(t):
	"""Calculates the intensity through space at a given time"""
	I_vec = c*tau_e/l * np.sqrt(4*np.log(2)/np.pi)*np.exp(-kappa_e*z)*np.exp( -4*np.log(2)*(t-t_G-z/c)**2/tau_G**2)
	I_vec[0]=I_vec[-1] = 0
	return I_vec

def I_R(t):
	"""Calculates the intensity through space at a given time"""
	I_vec = c*tau_e/l * np.sqrt(4*np.log(2)/np.pi)*np.exp(-kappa_e*z)*np.exp( -4*np.log(2)*(t-t_R-z/c)**2/tau_R**2)
	I_vec[0]=I_vec[-1] = 0
	return I_vec

for timestep in range(N):

	# get intensity spatial profile at current time
	I_G_mid_step = I_G(timestep*dt + dt)

	# calculate next time steps for variables
	N_pop_next = POP_RHS(N_pop, W_G, W_A)
	N_mid_step = (N_pop_next + N_pop)/2

	# Calculate premultiplying matrices for the pump and ASE
	A_pump = PUMP_LHS(N_mid_step)
	A_ASE = ASE_LHS(N_mid_step)

	W_G_next = np.linalg.solve(A_pump, PUMP_RHS(W_G, N_mid_step, I_G_mid_step) )
	W_A_next = np.linalg.solve(A_ASE, ASE_RHS(W_A, N_mid_step) )

	# set newly calculated values to current time for next loop
	N_pop = N_pop_next
	W_G = W_G_next
	W_A = W_A_next

	if timestep % 50 == 0:
		# store data in storage list
		I_G_storage.append(I_G_mid_step)
		N_pop_storage.append(N_pop)
		W_G_storage.append(W_G)
		W_A_storage.append(W_A)

	print(timestep, end='\r')


# convert storage lists to numpy arrays
I_G_storage = np.array(I_G_storage)
N_pop_storage = np.array(N_pop_storage)
W_G_storage = np.array(W_G_storage)
W_A_storage = np.array(W_A_storage)
Flux = (W_A_storage[:,0])
"""
np.savetxt('./Data//L=1/W_A.I=4e10.L=1.txt',W_A_storage, delimiter=',',newline='\n')
np.savetxt('./Data/L=1/N_pop.I=4e10.L=1.txt',N_pop_storage, delimiter=',',newline='\n')
np.savetxt('./Data/L=1/Flux.I=4e10.L=1.txt',Flux, delimiter=',',newline='\n')
"""