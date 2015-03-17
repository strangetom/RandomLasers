import numpy as np
from scipy import constants as spc

#define constants
sig_abs = 3.0e-24 # absorption cross section
sig_em = 3.0e-23 # emission cross section
tau_e = 3.2e-6 # lifetime of excited state
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
I_G0 = 2.e11 # average pump intensity

#define space parameters
L = 0.003 # length of medium
dz = l/2 # space increment
z = np.arange(-l/2, L+l/2+dz, dz) # vector in space
z[0] = z[1] # modifiy space vector so pulses don't decay before entering medium
z[-1] = z[-2] # modifiy space vector so pulses don't decay after exiting medium
J = z.shape[0] # number of space steps


# q is the number of space steps in from the sides that we're removing
q = 15
# change the shape of N_t
N_t = np.vstack( (np.zeros((q,J)), 1.6e25*np.ones((J-2*q,J)), np.zeros((q,J)) ) ) # total concentration of laser particles

F = np.vstack( (np.zeros((q,J)), np.ones((J-2*q,J)), np.zeros((q,J)) ) )


#define time parameters
N = 100000 # number of time steps
T = 50.e-9 # length of time
dt = T/N # time increment

# initial conditions
W_G = np.zeros((J,J))
W_A = np.zeros((J,J))
N_pop = np.zeros((J,J))

# lists for storage
W_G_storage = []
W_A_storage = []
N_pop_storage = []
I_G_storage = []

# function definitions
def PUMP_RHS(W_G, N_pop, I_G, y_deriv_constant):
	"""Calculates the right hand side of the PDE for the pump"""
	if y_deriv_constant:
		M = np.diagflat([0]+[ dt*D/(2*dz**2) for i in range(J-2)], -1) + np.diagflat([0]+[(1. - 2*dt*D/(2*dz**2)) for i in range(J-2) ]+[0]) + np.diagflat([ dt*D/(2*dz**2) for i in range(J-2)]+[0], 1)
		return M.dot(W_G.T).T*F - dt/4*sig_abs*v*N_t*(1-N_pop)*W_G*F + dt/tau_e *I_G*F
	elif ~y_deriv_constant:
		M = np.diagflat([0]*(q+1)+[ dt*D/(2*dz**2) for i in range(J-1-2*q)]+(q-1)*[0], -1) + np.diagflat([0]*q+[(1. - 2*dt*D/(2*dz**2)) for i in range(J-2*q) ]+q*[0]) + np.diagflat([0]*(q-1)+[ dt*D/(2*dz**2) for i in range(J-1-2*q)]+(q+1)*[0], 1)
		return M.dot(W_G)*F - dt/4*sig_abs*v*N_t*(1-N_pop)*W_G*F + dt/tau_e *I_G*F

def PUMP_LHS(N_mid_step, y_deriv_constant):
	"""Calculates the left hand side of the PDE for the pump"""
	if y_deriv_constant:
		M = np.diagflat([0]+[ -dt*D/(2*dz**2) for i in range(J-2)], -1) + np.diagflat([1]+[ (1. + 2*dt*D/(2*dz**2) + dt/4*sig_abs*v*N_t[i][i]*(1-N_mid_step[i][i])) for i in range(J-2)]+[1]) + np.diagflat([ -dt*D/(2*dz**2) for i in range(J-2)]+[0], 1)
		return M
	elif ~y_deriv_constant:
		M = np.diagflat([0]*(q+1)+[ -dt*D/(2*dz**2) for i in range(J-1-2*q)]+(q-1)*[0], -1) + np.diagflat([1]*q+[ (1. + 2*dt*D/(2*dz**2) + dt/4*sig_abs*v*N_t[i][i]*(1-N_mid_step[i][i])) for i in range(J-2*q)]+q*[1]) + np.diagflat([0]*(q-1)+[ -dt*D/(2*dz**2) for i in range(J-1-2*q)]+(q+1)*[0], 1)
		return M 

def ASE_RHS(W_A, N_mid_step, y_deriv_constant):
	"""Calculates the right hand side of the PDE for amplified spontaneous emission"""
	if y_deriv_constant:
		M = np.diagflat([0]+[ dt*D/(2*dz**2) for i in range(J-2)], -1) + np.diagflat([0]+[ (1. - 2*dt*D/(2*dz**2)) for i in range(J-2)]+[0]) + np.diagflat([ dt*D/(2*dz**2) for i in range(J-2)]+[0], 1)
		return M.dot(W_A.T).T*F + dt/4*sig_em*v*N_mid_step*N_t*W_A*F + dt*c*N_t*E_A/tau_e/I_G0 * N_mid_step*F
	elif ~y_deriv_constant:
		M = np.diagflat([0]*(q+1)+[ dt*D/(2*dz**2) for i in range(J-1-2*q)]+(q-1)*[0], -1) + np.diagflat([0]*q+[ (1. - 2*dt*D/(2*dz**2)) for i in range(J-2*q)]+q*[0]) + np.diagflat([0]*(q-1)+[ dt*D/(2*dz**2) for i in range(J-1-2*q)]+(q+1)*[0], 1)
		return M.dot(W_A)*F + dt/4*sig_em*v*N_mid_step*N_t*W_A*F + dt*c*N_t*E_A/tau_e/I_G0 * N_mid_step*F

def ASE_LHS(N_mid_step, y_deriv_constant):
	"""Calculates the left hand side of the PDE for amplified spontaneous emission"""
	if y_deriv_constant:
		M = np.diagflat([0]+[ -dt*D/(2*dz**2) for i in range(J-2)], -1) + np.diagflat([1]+[ (1. + 2*dt*D/(2*dz**2) - dt/4*sig_em*v*N_mid_step[i][i]*N_t[i][i]) for i in range(J-2)]+[1]) + np.diagflat([ -dt*D/(2*dz**2) for i in range(J-2)]+[0], 1)
		return M
	elif ~y_deriv_constant:
		M = np.diagflat([0]*(q+1)+[ -dt*D/(2*dz**2) for i in range(J-1-2*q)]+(q-1)*[0], -1) + np.diagflat([1]*q+[ (1. + 2*dt*D/(2*dz**2) - dt/4*sig_em*v*N_mid_step[i][i]*N_t[i][i]) for i in range(J-2*q)]+q*[1]) + np.diagflat([0]*(q-1)+[ -dt*D/(2*dz**2) for i in range(J-1-2*q)]+(q+1)*[0], 1)
		return M

def POP_RHS(N_pop, W_G, W_A):
	"""Calculates the right hand side of the PDE for the excited population"""
	term_1 = dt/2*sig_abs*v*(1-N_pop)*I_G0/(c*E_G) * W_G
	term_2 = dt/2*sig_em*v*N_pop*I_G0/c * W_A/E_A
	term_3 = dt/2/tau_e * N_pop
	return N_pop + (term_1 - term_2 - term_3)

def I_G(t):
	"""Calculates the intensity through space at a given time"""
	I_vec = c*tau_e/l * np.sqrt(4*np.log(2)/np.pi)*np.exp(-kappa_e*z)*np.exp( -4*np.log(2)*(t-t_G-z/c)**2/tau_G**2)
	I_vec[0]=I_vec[-1] = 0
	return np.zeros((J,J)) + I_vec.reshape(J,1)

for timestep in range(N):

	"""
	First we calculate the n+0.5 step.
	The y derivative is assumed constant for this step.
	"""

	# get intensity spatial profile at n+0.25
	I_G_mid_step = I_G(timestep*dt + dt/4)

	# calculate intermediate time steps for variables
	N_pop_next_half = POP_RHS(N_pop, W_G, W_A)
	N_mid_step = (N_pop_next_half + N_pop)/2

	# Calculate premultiplying matrices for the pump and ASE
	A_pump = PUMP_LHS(N_mid_step, True)
	A_ASE = ASE_LHS(N_mid_step, True)

	W_G_mid_step = np.linalg.solve(A_pump, PUMP_RHS(W_G, N_mid_step, I_G_mid_step, True) ) 
	W_A_mid_step = np.linalg.solve(A_ASE, ASE_RHS(W_A, N_mid_step, True) ) 

	"""
	Then we complete the calculation to get the n+1 step.
	The x derivative is assumed constant for this step.
	"""

	# get intensity spatial profile at n+0.75
	I_G_mid_step = I_G(timestep*dt + 3*dt/4)

	# calculate next time steps for variables
	N_pop_next = POP_RHS(N_pop_next_half, W_G_mid_step, W_A_mid_step)
	N_mid_step = (N_pop_next + N_pop_next_half)/2

	# Calculate premultiplying matrices for the pump and ASE
	A_pump = PUMP_LHS(N_mid_step, False)
	A_ASE = ASE_LHS(N_mid_step, False)

	W_G_next = (np.linalg.solve(A_pump, PUMP_RHS(W_G_mid_step, N_mid_step, I_G_mid_step, False) )).T
	W_A_next = (np.linalg.solve(A_ASE, ASE_RHS(W_A_mid_step, N_mid_step, False) )).T

	"""
	Set newly calculated data ready for next loop and save data
	"""
	# set newly calculated values to current time for next loop
	N_pop = N_pop_next
	W_G = W_G_next
	W_A = W_A_next

	if timestep % 100 == 0:
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
np.savetxt('./Data/L=2/W_A.I=4e10.L=2.txt',W_A_storage[:,int(J/2)], delimiter=',',newline='\n')
np.savetxt('./Data/L=2/N_pop.I=4e10.L=2.txt',N_pop_storage[:,int(J/2)], delimiter=',',newline='\n')
np.savetxt('./Data/L=2/Flux.I=4e10.L=2.txt',np.sum(Flux,axis=1), delimiter=',',newline='\n')
"""
np.save('N_pop', N_pop_storage)
np.save('W_A', W_A_storage)
np.save('W_G', W_G_storage)