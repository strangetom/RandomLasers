import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plt.hold(True)	

W_G_store = np.loadtxt("./Data/W_G_store.txt", delimiter=',')

l = 40.e-6
L = 0.8e-3 #length of slap
dx = l/2. #space steps across slap
x = np.arange(0,L+2*l+dx,dx) # array of space steps

N = 1000000 #number of time steps
T = 5.e-8 #length of time
dt = T/N #time steps
time = ( dt*500*np.arange(0, W_G_store.shape[0], 1) ) *1e9 #time array in nanoseconds

"""
#If you fancy it multicoloured
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),  
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),  
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),  
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),  
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)] 

for i in range(len(tableau20)):  
    r, g, b = tableau20[i]  
    tableau20[i] = (r / 255., g / 255., b / 255.) 
"""

for i in range(x.shape[0]):
	position = ( i*dx*np.ones(W_G_store.shape[0]) ) *1e3 # position in microns
	ax.plot(xs=time, ys=W_G_store[:,i] ,zs=position, zdir='y', color='r'])

ax.set_xlabel('Time (ns)')
ax.set_xlim(0., 30)
ax.set_ylabel('Position (mm)')
ax.set_ylim(0, 0.88)
ax.set_zlabel('Energy density')
ax.set_zlim(0, 0.004)
ax.set_title('Pump pulse energy density')


plt.show()
