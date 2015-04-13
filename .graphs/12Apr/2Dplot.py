import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm

plt.close('all')

# load data
npop = np.load('N_pop_storage.npy')
wa = np.load('W_A_storage.npy')

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
npop12_6 = npop[252, 32][1:-1]
npop13_0 = npop[260, 32][1:-1]
npop13_5 = npop[271, 32][1:-1]

snapshots = [npop12_2, npop12_6, npop13_0, npop13_5]
snapshot_labels = [r"$12.2\/ns$", r"$12.6\/ns$", r"$13.0\/ns$", r"$13.5\/ns$"]

fig, ax = plt.subplots(figsize=(12,6))
for i, snapshot in enumerate(snapshots):
	if i == 0 or i == 3:
		ax.plot(X, snapshot, color=plt.cm.RdBu(20+150*i), linewidth=1.5, label=snapshot_labels[i])
	else:
		ax.plot(X, snapshot, color='0.65', linewidth=1.5, label=snapshot_labels[i])

ax.set_xlim((0.25, 2.75))
ax.set_ylim((0.1, 0.2))
ax.set_ylabel('Excitation level', fontsize=16)
ax.set_xlabel('Position inside medium (mm)', fontsize=16)
ax.set_xticks(np.arange(0.25, 3., 0.25))
ax.set_yticks(np.arange(0.11, 0.21, 0.01))
ax.legend()
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')


fig.tight_layout()
fig.savefig('snapshots.png')
plt.show()