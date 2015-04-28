import numpy as np
from scipy import constants as spc
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

data = np.load("N_pop_storage.npy")

l = 100.e-6
L = 0.003  # length of slap
dx = l / 2.  # space steps across slap
X = np.arange(-l/2, L + l/2 + dx, dx)/1e-3  # array of space step
Y = X.copy()
X, Y = np.meshgrid(X, Y)

fig = plt.figure()
ax = fig.gca(projection='3d')

def update(i, ax, fig):
	ax.cla()
	Z = data[i*10]
	surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.jet, linewidth=0)
	ax.set_zlim(0, 1.1*data.max())
	ax.set_xlabel('X distance (mm)')
	ax.set_ylabel('Y distance (mm)')
	ax.set_zlabel('Excitation level')
	ax.set_xlim(-l/2/1e-3, L/1e-3+l/2/1e-3)
	ax.set_ylim(-l/2/1e-3, L/1e-3+l/2/1e-3)
	ax.w_xaxis.gridlines.set_lw(0.0)
	ax.w_yaxis.gridlines.set_lw(0.0)
	ax.w_zaxis.gridlines.set_lw(0.0)
	ax.set_title("Time: {:.3f} ns".format(i*5e-10/1e-9))
	return surf

anim = animation.FuncAnimation(fig, update, fargs=(ax, fig), frames=int(data.shape[0]/10) )

anim.save('N_pop.mp4', fps=5, bitrate=-1, codec="libx264")




