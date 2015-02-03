import numpy as np
from scipy import constants as spc

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
t_R = 5.e-9 #time of maxima of probe pulse 
illum_area = np.pi*(2.e-3)**2 #illumation area
I_G0 = 200.e-3/illum_area/1 #average pump intensity. !!USING E/tau_G is the peak power, not average. We need the frequency of the pulses.!!
I_R0 = 200.e-3/illum_area/1 #average probe intensity
n = 1.35 #average refractive index of medium
c = spc.c/n #speed of light in medium
v = spc.c/n #transport velocity
D = v*l/3. #diffusion coeffecient


J = 4 #number of space steps
L = 8.e-4 #length of slap
dx = L/J #space steps across slap
z = np.arange(0,L,dx) # array of space steps

N = 100000 #number of time steps
T = 1.e-7 #length of time
dt = T/N #time steps

#Define beta value
beta = D*dt/(dx**2)


def f(N_1, W_G, I_G_vals):
	return -1*sig_abs*v*(N_t-N_1)*W_G + I_G_vals/l


def g(N_1, W_R, I_R_vals):
	return sig_em*v*N_1*W_R + I_R_vals/l


def h(N_1, W_A):
	return sig_em*v*N_1*W_A + N_1/tau_e


def q(N_1, W_G, W_R, W_A):
	return sig_abs*v*(N_t-N_1)*W_G - sig_em*v*N_1*(W_R+W_A) - N_1/tau_e


def I_G(t):
	return I_G0*np.sqrt(4*np.log(2)/np.pi)*np.exp(-kappa_e*z)*np.exp(-4*np.log(2)*(t-t_G-z/c)**2/(tau_G**2))


def I_R(t):
	return I_R0*np.sqrt(4*np.log(2)/np.pi)*np.exp(-kappa_e*z)*np.exp(-4*np.log(2)*(t-t_R-z/c)**2/(tau_R**2))


def create_B_matrix(beta):
	"""Defines the matrix multipling the laplacian """
	return np.diagflat([beta for i in range(J-1)], -1) + np.diagflat([0]+[1.-2.*beta for i in range(J-2)]+[0]) + np.diagflat([beta for i in range(J-1)], 1)

#Define intial conditions
W_G = np.zeros(z.shape[0])
W_R = np.zeros(z.shape[0])
W_A = np.zeros(z.shape[0])
N_1 = np.zeros(z.shape[0])

B_G = create_B_matrix(beta)
B_R = create_B_matrix(beta)
B_A = create_B_matrix(beta)


#Run numerical calculation
for timestep in range(N):
	I_G_vals = I_G(timestep*dt)
	I_R_vals = I_R(timestep*dt)

	W_G_new = B_G.dot(W_G) + dt*f(N_1,W_G,I_G_vals)
	W_R_new = B_R.dot(W_R) + dt*g(N_1,W_R,I_R_vals)
	W_A_new = B_A.dot(W_A) + dt*h(N_1,W_A)
	N_1_new = N_1 + dt*q(N_1, W_G, W_R, W_A)
	

	W_G = W_G_new
	W_R = W_R_new
	W_A = W_A_new
	N_1 = N_1_new

	if timestep*dx == 0:
		W_G_store = W_G
		W_R_store = W_R
		W_A_store = W_A
		N_1_store = N_1
		I_G_store = I_G_vals
		I_R_store = I_R_vals
	elif timestep*dx < 5e-8:
		W_G_store = np.vstack((W_G_store, W_G))
		W_R_store = np.vstack((W_R_store, W_R))
		W_A_store = np.vstack((W_A_store, W_A))
		N_1_store = np.vstack((N_1_store, N_1))
		I_G_store = np.vstack((I_G_store, I_G_vals))
		I_R_store = np.vstack((I_R_store, I_R_vals))
	else:
		W_A_store = np.vstack((W_A_store, W_A))
		N_1_store = np.vstack((N_1_store, N_1))

	print(timestep,end='\r')
