import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm

plt.close('all')

# load data
npop = np.load('N_pop_storage.npy')
wa = np.load('W_A_storage.npy')

# some colours
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),  
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),  
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),  
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),  
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)] 
for i in range(len(tableau20)):  
    r, g, b = tableau20[i]  
    tableau20[i] = (r / 255., g / 255., b / 255.)          

# spatial grid
l = 100.e-6
L = 0.003  # length of slap
dx = l / 2.  # space steps across slap
X = np.arange(-l/2, L + l/2 + dx, dx)/1e-3  # array of space step
X = X[1:-1] # restrict from 0 mm to 3 mm

# time grid (not needed)
N = 100000 # number of time steps
T = 50.e-9 # length of time
dt = T/N # time increment
time = dt*100*np.arange(0, 1000, 1)

# take snapshots at 12.9ns (max) and at 14.5ns (min) 
# 20 rows = 1 ns
npop12_2 = npop[244, 32][1:-1]
npop12_5 = npop[251, 32][1:-1]
npop12_9 = npop[258, 32][1:-1]
npop13_2 = npop[265, 32][1:-1]
npop13_5 = npop[271, 32][1:-1]

snapshots = [npop12_2, npop12_5, npop12_9, npop13_2, npop13_5]
snapshot_labels = [r"$12.2\/ns$", r"$12.5\/ns$", r"$12.9\/ns$", r"$13.2\/ns$", r"$13.5\/ns$"]
greys = ['0.3', '0.5', '0.7']

fig, ax = plt.subplots(figsize=(12,6))
for i, snapshot in enumerate(snapshots):
	if i == 0:
		ax.plot(X, snapshot, color=tableau20[0], linewidth=1.5, label=snapshot_labels[i])
	elif i == 4:
		ax.plot(X, snapshot, color=tableau20[6], linewidth=1.5, label=snapshot_labels[i])
	else:
		ax.plot(X, snapshot, color=greys[i-1], linewidth=1.5, label=snapshot_labels[i])

ax.set_xlim((0.25, 2.75))
ax.set_ylim((0.1, 0.2))
ax.set_ylabel('Excitation level', fontsize=16)
ax.set_xlabel('Position inside medium (mm)', fontsize=16)
ax.set_xticks(np.arange(0.25, 3., 0.25))
ax.set_yticks(np.arange(0.11, 0.21, 0.01))
lg = ax.legend()
lg.draw_frame(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')


fig.tight_layout()
fig.savefig('snapshots.png')
plt.show()