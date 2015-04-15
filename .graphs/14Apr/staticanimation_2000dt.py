import numpy as np
from scipy import constants as spc
from matplotlib import pyplot as plt

# Tableau20 definitions
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),  
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),  
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),  
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),  
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)] 
for i in range(len(tableau20)):  
    r, g, b = tableau20[i]  
    tableau20[i] = (r / 255., g / 255., b / 255.)

L = 50e-6 # approx. length of medium
dx = 1e-9 # space step
T = 7e-13 # time for simulation
dt = dx/spc.c # time step
N = int(T/dt)

medium = np.loadtxt('Medium.txt', delimiter=',')
E_storage = np.loadtxt('E_storage.txt', delimiter=',')

space = dx*np.arange(0,E_storage.shape[1], 1)

fig, ax = plt.subplots(5, 1, figsize=(16,9))
fig.subplots_adjust(hspace = .007, wspace = .001)

ax[0].plot(space/1e-6, E_storage[30], color=tableau20[6])
ax[0].fill_between(space/1e-6, 1,-1, where=medium, color='0.87')
#ax[0].barh(-0.2, space.shape[0]*dx/1e-6, height=0.4, left=0, align='edge', color='0.85', alpha=0.4, ec='None')
ax[0].set_yticks([-0.2, 0, 0.2, 0.4])
ax[0].set_xticklabels(())
ax[0].spines["top"].set_visible(False)
ax[0].spines["right"].set_visible(False)
ax[0].yaxis.set_ticks_position('left')
ax[0].tick_params(axis='x', bottom='off', top='off')
ax[0].set_ylim((-0.3, 0.5))
ax[0].set_xlim((0, 0.9))
ax[0].set_ylabel('0.1 fs', rotation=0, fontsize=15)
ax[0].yaxis.set_label_coords(-0.065, 0.5)

ax[1].plot(space/1e-6, E_storage[100], color=tableau20[6])
ax[1].fill_between(space/1e-6, 1,-1, where=medium, color='0.87')
ax[1].set_yticks([-0.2, 0, 0.2, 0.4])
ax[1].set_xticklabels(())
ax[1].spines["top"].set_visible(False)
ax[1].spines["right"].set_visible(False)
ax[1].yaxis.set_ticks_position('left')
ax[1].tick_params(axis='x', bottom='off', top='off')
ax[1].set_ylim((-0.3, 0.5))
ax[1].set_xlim((0, 0.9))
ax[1].set_ylabel('0.3 fs', rotation=0, fontsize=15)
ax[1].yaxis.set_label_coords(-0.065, 0.5)

ax[2].plot(space/1e-6, E_storage[300], color=tableau20[6])
ax[2].fill_between(space/1e-6, 1,-1, where=medium, color='0.87')
ax[2].set_yticks([-0.2, 0, 0.2, 0.4])
ax[2].set_xticklabels(())
ax[2].spines["top"].set_visible(False)
ax[2].spines["right"].set_visible(False)
#ax[2].spines["bottom"].set_visible(False)
ax[2].yaxis.set_ticks_position('left')
ax[2].tick_params(axis='x', bottom='off', top='off')
ax[2].set_ylim((-0.3, 0.5))
ax[2].set_xlim((0, 0.9))
ax[2].set_ylabel('1.0 fs', rotation=0, fontsize=15)
ax[2].yaxis.set_label_coords(-0.065, 0.5)

ax[3].plot(space/1e-6, E_storage[1000], color=tableau20[6])
ax[3].fill_between(space/1e-6, 1,-1, where=medium, color='0.87')
ax[3].set_yticks([-0.2, 0, 0.2, 0.4])
ax[3].set_xticklabels(())
ax[3].spines["top"].set_visible(False)
ax[3].spines["right"].set_visible(False)
#ax[3].spines["bottom"].set_visible(False)
ax[3].yaxis.set_ticks_position('left')
ax[3].tick_params(axis='x', bottom='off', top='off')
ax[3].set_ylim((-0.3, 0.5))
ax[3].set_xlim((0, 0.9))
ax[3].set_ylabel('3.3 fs', rotation=0, fontsize=15)
ax[3].yaxis.set_label_coords(-0.065, 0.5)

ax[4].plot(space/1e-6, E_storage[2000], color=tableau20[6])
ax[4].fill_between(space/1e-6, 1,-1, where=medium, color='0.87')
ax[4].set_yticks([-0.2, 0, 0.2, 0.4])
ax[4].spines["top"].set_visible(False)
ax[4].spines["right"].set_visible(False)
ax[4].yaxis.set_ticks_position('left')
ax[4].tick_params(axis='x', top='off')
ax[4].set_ylim((-0.3, 0.5))
ax[4].set_xticklabels(np.arange(0, 1, 0.1), fontsize=14)
ax[4].set_xlim((0, 0.9))
ax[4].set_xlabel('Distance ($\mu$m)', fontsize=16)
ax[4].set_ylabel('6.6 fs', rotation=0, fontsize=15)
ax[4].yaxis.set_label_coords(-0.065, 0.5)

fig.text(0.04, 0.5, 'E (Vm$^{-1}$)', ha='center', va='center',rotation='vertical',fontsize=16)

fig.savefig('temp.png')
#plt.show()
