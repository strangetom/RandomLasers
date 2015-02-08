import numpy as np
from scipy import constants as spc
from scipy import sparse as sprs

#define constants
sig_abs = 3.0e-24 #absorption cross section
sig_em = 3.0e-23 #emission cross section
tau_e = 3.2e-6 #lifetime of excited state
N_t = 1.6e25 #total concentration of laser particles
l = 40.e-6 #transport mean free path
kappa_e = 3.e3 #extinction coefficient (see spec sheet in shared folder)
tau_G = 14.e-9 #pump pulse FWHM
tau_R = 14.e-9 #proble pulse FWHM 
t_G = 5.e-9 #time of maxima of pump pulse 
t_R = 15.e-9 #time of maxima of probe pulse 
illum_area = np.pi*(2.e-3)**2 #illumation area
I_G0 = 200e-3/illum_area/1 #average pump intensity. !!USING E/tau_G is the peak power, not average. We need the frequency of the pulses.!!
I_R0 = 200e-3/illum_area/1 #average probe intensity
n = 1.35 #average refractive index of medium
c = spc.c/n #speed of light in medium
v = spc.c/n #transport velocity
D = v*l/3. #diffusion coeffecient


L = 0.0008 #length of slap
dx = l/2. #space steps across slap
x = np.arange(-l,L+l+dx,dx) #array of space steps
x[0] = x[1] = 0#modify so intensity doesn't decay before medium
x[-2] = x[-1] = x[-3] #modify so intensity doesn't decay after medium
J = x.shape[0]

N = 1000000 #number of time steps
T = 5.e-8 #length of time
dt = T/N #time steps

#Define beta value
beta = D*dt/(dx**2)


def f(N_1, W_G, I_G_vals):
	"""None partial derivate function for W_G"""
	return -1*sig_abs*v*(N_t-N_1)*W_G + I_G_vals/l

def g(N_1, W_R, I_R_vals):
	"""None partial derivate function for W_R"""
	return sig_em*v*N_1*W_R + I_R_vals/l

def h(N_1, W_A):
	"""None partial derivate function for W_A"""
	return sig_em*v*N_1*W_A + N_1/tau_e

def q(N_1, W_G, W_R, W_A):
	"""None partial derivate function for N_!"""
	return sig_abs*v*(N_t-N_1)*W_G - sig_em*v*N_1*(W_R+W_A) - N_1/tau_e

def I_G(t):
	"""Calculate pump intensity vector at time t"""
	return I_G0*np.sqrt(4*np.log(2)/np.pi)*np.exp(-kappa_e*x)*np.exp(-4*np.log(2)*(t-t_G-x/c)**2/(tau_G**2))

def I_R(t):
	"""Calculate probe intensity vector at time t"""
	return I_R0*np.sqrt(4*np.log(2)/np.pi)*np.exp(-kappa_e*x)*np.exp(-4*np.log(2)*(t-t_R-x/c)**2/(tau_R**2))

def create_B_matrix(beta):
	"""Defines the matrix multipling the laplacian """
	return np.diagflat([beta for i in range(J-1)], -1) + np.diagflat([0]+[1.-2.*beta for i in range(J-2)]+[0]) + np.diagflat([beta for i in range(J-1)], 1)

if dt > dx**2/(2*D):
	#If stability criterion not met, abort.
	print("Unstable conditions.")
	print(str(dt)+" !< "+str(dx**2/(2*D)) )
else:
	#Define intial conditions
	W_G = np.zeros(x.shape[0])
	W_R = np.zeros(x.shape[0])
	W_A = np.zeros(x.shape[0])
	N_1 = np.zeros(x.shape[0])

	B = create_B_matrix(beta)

	#Storage lists (containing intial values)
	W_G_store = []
	W_R_store = []
	W_A_store = []
	N_1_store = []
	I_G_store = []
	I_R_store = []
	Outgoing_flux = [(W_A[2]-W_A[1])/dx]

	#Run numerical calculation
	for timestep in range(N):
		I_G_vals = I_G(timestep*dt)
		I_R_vals = I_R(timestep*dt)

		W_G_new = B.dot(W_G) + dt*f(N_1,W_G,I_G_vals)
		W_R_new = B.dot(W_R) + dt*g(N_1,W_R,I_R_vals)
		W_A_new = B.dot(W_A) + dt*h(N_1,W_A)
		N_1_new = N_1 + dt*q(N_1, W_G, W_R, W_A)

		W_G = W_G_new
		W_R = W_R_new
		W_A = W_A_new
		N_1 = N_1_new

		if timestep % 500 == 0:	
			if timestep*dt < 30e-9:
				W_G_store.append(W_G)
				W_R_store.append(W_R)
				W_A_store.append(W_A)
				N_1_store.append(N_1)
				I_G_store.append(I_G_vals)
				I_R_store.append(I_R_vals)
				Outgoing_flux.append((W_A[2]-W_A[1])/dx)
			else:
				W_A_store.append(W_A)
				N_1_store.append(N_1)
				Outgoing_flux.append((W_A[2]-W_A[1])/dx)

		print(timestep,end='\r')

	W_G_store = np.array(W_G_store)
	W_R_store = np.array(W_R_store)
	W_A_store = np.array(W_A_store)
	N_1_store = np.array(N_1_store)
	I_G_store = np.array(I_G_store)
	I_R_store = np.array(I_R_store)
	Outgoing_flux = np.array(Outgoing_flux)

"""
np.savetxt('./Data/W_G_store.txt',W_G_store, delimiter=',',newline='\n')
np.savetxt('./Data/W_R_store.txt',W_R_store, delimiter=',',newline='\n')
np.savetxt('./Data/W_A_store.txt',W_A_store, delimiter=',',newline='\n')
np.savetxt('./Data/N_1_store.txt',N_1_store, delimiter=',',newline='\n')
np.savetxt('./Data/I_G_store.txt',I_G_store, delimiter=',',newline='\n')
np.savetxt('./Data/I_R_store.txt',I_R_store, delimiter=',',newline='\n')
np.savetxt('./Data/Outgoing_flux.txt',Outgoing_flux, delimiter=',',newline='\n')
"""