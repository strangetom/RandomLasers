import numpy as np
from numpy import unravel_index
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

wa = np.load('W_A_storage.npy')
npop = np.load('N_pop_storage.npy')

Z1 = npop[unravel_index(wa.argmax(), npop.shape)[0]] #position of first peak
Z2 = npop[unravel_index(wa[360:440].argmax(), npop.shape)[0]+360]

N = 100000 # number of time steps
T = 50.e-9 # length of time
dt = T/N # time increment

l = 100.e-6
L = 0.003  # length of slap
dx = l / 2.  # space steps across slap
X = np.arange(-l/2, L + l/2 + dx, dx)/1e-3  # array of space step
Y = X.copy()
X, Y = np.meshgrid(X, Y)

fig = plt.figure(figsize=(18,8))
ax = fig.add_subplot(1, 2, 1, projection='3d')
ax.plot_surface(X, Y, Z1, rstride=1, cstride=1, linewidth=0, cmap=cm.rainbow)
ax.set_zlim(0, 1.1*Z1.max())
ax.set_xlabel('X distance (mm)')
ax.set_ylabel('Y distance (mm)')
ax.set_zlabel('Excitation level')
ax.set_xlim(-l/2/1e-3, L/1e-3+l/2/1e-3)
ax.set_ylim(-l/2/1e-3, L/1e-3+l/2/1e-3)
ax.set_title('Time: 12.9 ns')
ax.w_xaxis.gridlines.set_lw(1)
ax.w_yaxis.gridlines.set_lw(1)
ax.w_zaxis.gridlines.set_lw(1)
ax.elev = 30
ax = fig.add_subplot(1, 2, 2, projection='3d')
ax.plot_surface(X, Y, Z2, rstride=1, cstride=1, linewidth=0, cmap=cm.rainbow)
ax.set_zlim(0, 1.1*Z2.max())
ax.set_xlabel('X distance (mm)')
ax.set_ylabel('Y distance (mm)')
ax.set_zlabel('Excitation level')
ax.set_xlim(-l/2/1e-3, L/1e-3+l/2/1e-3)
ax.set_ylim(-l/2/1e-3, L/1e-3+l/2/1e-3)
ax.set_title('Time: 19.5 ns')
ax.w_xaxis.gridlines.set_lw(1)
ax.w_yaxis.gridlines.set_lw(1)
ax.w_zaxis.gridlines.set_lw(1)
ax.elev = 30
fig.subplots_adjust(hspace = .05, wspace = .01)
plt.savefig('peaks1and3.png')