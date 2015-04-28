import numpy as np
from scipy import constants as spc
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

data_1 = np.load("N_pop_storage.npy")
data_2 = np.load("W_A_storage.npy")

l = 100.e-6
L = 0.003  # length of slap
dx = l / 2.  # space steps across slap
X = np.arange(-l/2, L + l/2 + dx, dx)/1e-3  # array of space step
Y = X.copy()
X, Y = np.meshgrid(X, Y)

fig = plt.figure(figsize=(16,9))
ax1 = fig.add_subplot(1,2,1, projection='3d')
ax2 = fig.add_subplot(1,2,2, projection='3d')
fig.subplots_adjust(wspace = .05)
figtext = fig.text(0.47, 0.93, "", fontsize=18)

def update(i, ax1, ax2, fig):
	ax1.cla()
	ax2.cla()

	figtext.set_text("Time: {:.3f} ns".format(i*5e-10/1e-9))

	Z_1 = data_1[i*10]
	surf = ax1.plot_surface(X, Y, Z_1, rstride=1, cstride=1, cmap=cm.RdBu, linewidth=0)
	ax1.set_zlim(0, 1.1*data_1.max())
	ax1.set_xlabel('X distance (mm)')
	ax1.set_ylabel('Y distance (mm)')
	ax1.set_zlabel('Excitation level')
	ax1.set_xlim(-l/2/1e-3, L/1e-3+l/2/1e-3)
	ax1.set_ylim(-l/2/1e-3, L/1e-3+l/2/1e-3)
	ax1.w_xaxis.gridlines.set_lw(0.0)
	ax1.w_yaxis.gridlines.set_lw(0.0)
	ax1.w_zaxis.gridlines.set_lw(0.0)
	ax1.set_title("Population")
	ax1.view_init(elev=15.)

	Z_2 = data_2[i*10]
	surf = ax2.plot_surface(X, Y, Z_2, rstride=1, cstride=1, cmap=cm.RdBu, linewidth=0)
	ax2.set_zlim(0, 1.1*data_2.max())
	ax2.set_xlabel('X distance (mm)')
	ax2.set_ylabel('Y distance (mm)')
	ax2.set_zlabel('Energy density')
	ax2.set_xlim(-l/2/1e-3, L/1e-3+l/2/1e-3)
	ax2.set_ylim(-l/2/1e-3, L/1e-3+l/2/1e-3)
	ax2.w_xaxis.gridlines.set_lw(0.0)
	ax2.w_yaxis.gridlines.set_lw(0.0)
	ax2.w_zaxis.gridlines.set_lw(0.0)
	ax2.set_title("ASE")
	ax2.view_init(elev=15.)

anim = animation.FuncAnimation(fig, update, fargs=(ax1, ax2, fig), frames=int(data_1.shape[0]/10) )

anim.save('W_A+N_pop.mp4', fps=5, bitrate=-1, codec="libx264")
#plt.show()



