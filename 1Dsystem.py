import numpy as np
import scipy.constants as spc

#define constants
sig_abs = 3.0e-24 #absorption cross section
sig_em = 3.0e-23 #emission cross section
tau_e = 3.2e-6 #lifetime of excited state
N_t = 1.6e25 #total concentration of laser particles
l = 40.e-6 #transport mean free path
kappa_e = 0.1e-9 #extinction coefficient
tau_G = 14.e-9 #pump pulse FWHM
tau_R = 14.e-9 #proble pulse FWHM 
t_G = 0. #time of maxima of pump pulse 
t_R = 0. #time of maxima of probe pulse 
illum_area = np.pi*(2.e-3)**2 #illumation area
I_G0 = 200.e-3/illum_area/tau_G #average pump intensity
I_R0 = 200.e-3/illum_area/tau_R #average probe intensity
n = 1.35 #average refractive index of medium
c = spc.c/n #speed of light in medium
v = spc.c/n #transport velocity
D = v*l/3. #diffusion coeffecient


J = 100 #number of space steps
L = 8.e-4 #length of slap
dx = L/J #space steps across slap
z = np.arange(0,L-dx,dx) # array of space steps

N = 100 #number of time steps
T = 1e-6 #length of time
dt = T/N #time steps



def f(N_1, W_G, I_G_vals):
	return -1*sig_abs*v*(N_t-N_1)*W_G + I_G_vals/l

def g(N_1, W_R, I_R_vals):
	return sig_em*v*N_1*W_R + I_R_vals/l

def h(N_1, W_A):
	return sig_em*v*N_1*W_A + N_1/tau_e

def q(N_1, W_G, W_R, W_A):
	return sig_abs*v*(N_t-N_1)*W_G - sig_em*v*N_1*(W_R+W_A) - N_1/tau_e

def I_G(t):
	return I_R0*np.sqrt(4*np.log(2)/np.pi)*np.exp(-kappa_e*z)*np.exp(-4*np.log(2)*(t-t_G-z/c)**2/(tau_R**2))

def I_R(t):
	return I_R0*np.sqrt(4*np.log(2)/np.pi)*np.exp(-kappa_e*z)*np.exp(-4*np.log(2)*(t-t_R-z/c)**2/(tau_R**2))

def create_A_matrix(beta):
	"""Defines the matrix multipling the time derivative and inverts it"""
	A =  np.diagflat([-beta for i in range(J-1)], -1) + np.diagflat([1.+beta]+[1.+2.*beta for i in range(J-2)]+[1.+beta]) + np.diagflat([-beta for i in range(J-1)], 1)
	return np.linalg.inv(A)
	#consider sparse matrices? using scipy.sparse.diags http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.sparse.diags.html

def create_B_matrix(beta):
	"""Defines the matrix multipling the laplacian """
	return np.diagflat([beta for i in range(J-1)], -1) + np.diagflat([1.-beta]+[1.-2.*beta for i in range(J-2)]+[1.-beta]) + np.diagflat([beta for i in range(J-1)], 1)

#Define intial conditions
W_G = np.zeros(z.shape[0])
W_R = np.zeros(z.shape[0])
W_A = np.zeros(z.shape[0])
N_1 = np.zeros(z.shape[0])

#Define beta values
beta = D*dt/(2.*dx**2)

#Declare matrices for Crank-Nicolson method
A_G = create_A_matrix(beta)
A_R = create_A_matrix(beta)
A_A = create_A_matrix(beta)

B_G = create_B_matrix(beta)
B_R = create_B_matrix(beta)
B_A = create_B_matrix(beta)

#Run numerical calculation
for timestep in range(N):
	I_G_vals = I_G(timestep*dt)
	I_R_vals = I_R(timestep*dt)


	W_G_new = np.dot(A_G, (B_G.dot(W_G) + f(N_1,W_G,I_G_vals)) )
	W_R_new = np.dot(A_R, (B_R.dot(W_R) + g(N_1,W_R,I_R_vals)) )
	W_A_new = np.dot(A_A,  (B_A.dot(W_A) + h(N_1,W_A)) )
	N_1_new = N_1 + dt*q(N_1, W_G, W_R, W_A)

	W_G = W_G_new
	W_R = W_R_new
	W_A = W_A_new
	N_1 = N_1_new

	if timestep == 0:
		W_G_store = W_G
		W_R_store = W_R
		W_A_store = W_A
		N_1_store = N_1
	else:
		W_G_store = np.vstack((W_G_store, W_G))
		W_R_store = np.vstack((W_R_store, W_R))
		W_A_store = np.vstack((W_A_store, W_A))
		N_1_store = np.vstack((N_1_store, N_1))

	print(timestep,end='\r')
