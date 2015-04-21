import numpy as np
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

N_t = 1.6e25*(1e-2)**3
sig_em = 3e-23

n_L2_2IG0 = np.loadtxt('./Data/L=2/N_pop.I=2e10.L=2.txt', delimiter=',')
n_L2_4IG0 = np.loadtxt('./Data/L=2/N_pop.I=4e10.L=2.txt', delimiter=',')
n_L2_8IG0 = np.loadtxt('./Data/L=2/N_pop.I=8e10.L=2.txt', delimiter=',')

n_L3_2IG0 = np.loadtxt('./Data/L=3/N_pop.I=2e10.L=3.txt', delimiter=',')
n_L3_4IG0 = np.loadtxt('./Data/L=3/N_pop.I=4e10.L=3.txt', delimiter=',')
n_L3_8IG0 = np.loadtxt('./Data/L=3/N_pop.I=8e10.L=3.txt', delimiter=',')

n_L4_2IG0 = np.loadtxt('./Data/L=4/N_pop.I=2e10.L=4.txt', delimiter=',')
n_L4_4IG0 = np.loadtxt('./Data/L=4/N_pop.I=4e10.L=4.txt', delimiter=',')
n_L4_8IG0 = np.loadtxt('./Data/L=4/N_pop.I=8e10.L=4.txt', delimiter=',')

n_L2_2IG0 = N_t*sig_em*np.mean(np.mean(n_L2_2IG0, axis=1), axis=0)
n_L2_4IG0 = N_t*sig_em*np.mean(np.mean(n_L2_4IG0, axis=1), axis=0)
n_L2_8IG0 = N_t*sig_em*np.mean(np.mean(n_L2_8IG0, axis=1), axis=0)

n_L3_2IG0 = N_t*sig_em*np.mean(np.mean(n_L3_2IG0, axis=1), axis=0)
n_L3_4IG0 = N_t*sig_em*np.mean(np.mean(n_L3_4IG0, axis=1), axis=0)
n_L3_8IG0 = N_t*sig_em*np.mean(np.mean(n_L3_8IG0, axis=1), axis=0)

n_L4_2IG0 = N_t*sig_em*np.mean(np.mean(n_L4_2IG0, axis=1), axis=0)
n_L4_4IG0 = N_t*sig_em*np.mean(np.mean(n_L4_4IG0, axis=1), axis=0)
n_L4_8IG0 = N_t*sig_em*np.mean(np.mean(n_L4_8IG0, axis=1), axis=0)


fig, ax = plt.subplots(figsize=(12,8))

ax.plot([2, 3, 4], [n_L2_2IG0*1e3, n_L3_2IG0*1e3, n_L4_2IG0*1e3], color=tableau20[0], label=r'$I_{G0}=2\times 10^{10}\/Wm^{-2}$')
ax.plot([2, 3, 4], [n_L2_4IG0*1e3, n_L3_4IG0*1e3, n_L4_4IG0*1e3], color=tableau20[6], label=r'$I_{G0}=4\times 10^{10}\/Wm^{-2}$')
ax.plot([2, 3, 4], [n_L2_8IG0*1e3, n_L3_8IG0*1e3, n_L4_8IG0*1e3], color=tableau20[2], label=r'$I_{G0}=8\times 10^{10}\/Wm^{-2}$')
ax.set_xlim((1.9, 4.1))
#ax.set_ylim((0, 0.11))
ax.set_xlabel('Medium dimensions (mm)', fontsize=16)
ax.set_ylabel('Average gain coefficient ($10^{-3}cm^{-1}$)', fontsize=16)

lg = ax.legend()
lg.draw_frame(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

fig.savefig('Averagegain2D.png')
#plt.show()