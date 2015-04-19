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

# 1D data
n_L1_IG0_1D = np.loadtxt('./Data/L=1/N_pop.I=1e11.L=1.txt', delimiter=',')
n_L1_2IG0_1D = np.loadtxt('./Data/L=1/N_pop.I=2e11.L=1.txt', delimiter=',')
n_L1_4IG0_1D = np.loadtxt('./Data/L=1/N_pop.I=4e11.L=1.txt', delimiter=',')
n_L1_8IG0_1D = np.loadtxt('./Data/L=1/N_pop.I=8e11.L=1.txt', delimiter=',')

n_L2_IG0_1D = np.loadtxt('./Data/L=2/N_pop.I=1e11.L=2.txt', delimiter=',')
n_L2_2IG0_1D = np.loadtxt('./Data/L=2/N_pop.I=2e11.L=2.txt', delimiter=',')
n_L2_4IG0_1D = np.loadtxt('./Data/L=2/N_pop.I=4e11.L=2.txt', delimiter=',')
n_L2_8IG0_1D = np.loadtxt('./Data/L=2/N_pop.I=8e11.L=2.txt', delimiter=',')

n_L3_IG0_1D = np.loadtxt('./Data/L=3/N_pop.I=1e11.L=3.txt', delimiter=',')
n_L3_2IG0_1D = np.loadtxt('./Data/L=3/N_pop.I=2e11.L=3.txt', delimiter=',')
n_L3_4IG0_1D = np.loadtxt('./Data/L=3/N_pop.I=4e11.L=3.txt', delimiter=',')
n_L3_8IG0_1D = np.loadtxt('./Data/L=3/N_pop.I=8e11.L=3.txt', delimiter=',')

n_L4_IG0_1D = np.loadtxt('./Data/L=4/N_pop.I=1e11.L=4.txt', delimiter=',')
n_L4_2IG0_1D = np.loadtxt('./Data/L=4/N_pop.I=2e11.L=4.txt', delimiter=',')
n_L4_4IG0_1D = np.loadtxt('./Data/L=4/N_pop.I=4e11.L=4.txt', delimiter=',')
n_L4_8IG0_1D = np.loadtxt('./Data/L=4/N_pop.I=8e11.L=4.txt', delimiter=',')

n_L1_IG0_1D = N_t*sig_em*np.mean(np.mean(n_L1_IG0_1D, axis=1), axis=0)
n_L1_2IG0_1D = N_t*sig_em*np.mean(np.mean(n_L1_2IG0_1D, axis=1), axis=0)
n_L1_4IG0_1D = N_t*sig_em*np.mean(np.mean(n_L1_4IG0_1D, axis=1), axis=0)
n_L1_8IG0_1D = N_t*sig_em*np.mean(np.mean(n_L1_8IG0_1D, axis=1), axis=0)

n_L2_IG0_1D = N_t*sig_em*np.mean(np.mean(n_L2_IG0_1D, axis=1), axis=0)
n_L2_2IG0_1D = N_t*sig_em*np.mean(np.mean(n_L2_2IG0_1D, axis=1), axis=0)
n_L2_4IG0_1D = N_t*sig_em*np.mean(np.mean(n_L2_4IG0_1D, axis=1), axis=0)
n_L2_8IG0_1D = N_t*sig_em*np.mean(np.mean(n_L2_8IG0_1D, axis=1), axis=0)

n_L3_IG0_1D = N_t*sig_em*np.mean(np.mean(n_L3_IG0_1D, axis=1), axis=0)
n_L3_2IG0_1D = N_t*sig_em*np.mean(np.mean(n_L3_2IG0_1D, axis=1), axis=0)
n_L3_4IG0_1D = N_t*sig_em*np.mean(np.mean(n_L3_4IG0_1D, axis=1), axis=0)
n_L3_8IG0_1D = N_t*sig_em*np.mean(np.mean(n_L3_8IG0_1D, axis=1), axis=0)

n_L4_IG0_1D = N_t*sig_em*np.mean(np.mean(n_L4_IG0_1D, axis=1), axis=0)
n_L4_2IG0_1D = N_t*sig_em*np.mean(np.mean(n_L4_2IG0_1D, axis=1), axis=0)
n_L4_4IG0_1D = N_t*sig_em*np.mean(np.mean(n_L4_4IG0_1D, axis=1), axis=0)
n_L4_8IG0_1D = N_t*sig_em*np.mean(np.mean(n_L4_8IG0_1D, axis=1), axis=0)

# 2D data
n_L2_2IG0_2D = np.loadtxt('./Data/L=2/N_pop.I=2e10.L=2.txt', delimiter=',')
n_L2_4IG0_2D = np.loadtxt('./Data/L=2/N_pop.I=4e10.L=2.txt', delimiter=',')
n_L2_8IG0_2D = np.loadtxt('./Data/L=2/N_pop.I=8e10.L=2.txt', delimiter=',')

n_L3_2IG0_2D = np.loadtxt('./Data/L=3/N_pop.I=2e10.L=3.txt', delimiter=',')
n_L3_4IG0_2D = np.loadtxt('./Data/L=3/N_pop.I=4e10.L=3.txt', delimiter=',')
n_L3_8IG0_2D = np.loadtxt('./Data/L=3/N_pop.I=8e10.L=3.txt', delimiter=',')

n_L4_2IG0_2D = np.loadtxt('./Data/L=4/N_pop.I=2e10.L=4.txt', delimiter=',')
n_L4_4IG0_2D = np.loadtxt('./Data/L=4/N_pop.I=4e10.L=4.txt', delimiter=',')
n_L4_8IG0_2D = np.loadtxt('./Data/L=4/N_pop.I=8e10.L=4.txt', delimiter=',')

n_L2_2IG0_2D = N_t*sig_em*np.mean(np.mean(n_L2_2IG0_2D, axis=1), axis=0)
n_L2_4IG0_2D = N_t*sig_em*np.mean(np.mean(n_L2_4IG0_2D, axis=1), axis=0)
n_L2_8IG0_2D = N_t*sig_em*np.mean(np.mean(n_L2_8IG0_2D, axis=1), axis=0)

n_L3_2IG0_2D = N_t*sig_em*np.mean(np.mean(n_L3_2IG0_2D, axis=1), axis=0)
n_L3_4IG0_2D = N_t*sig_em*np.mean(np.mean(n_L3_4IG0_2D, axis=1), axis=0)
n_L3_8IG0_2D = N_t*sig_em*np.mean(np.mean(n_L3_8IG0_2D, axis=1), axis=0)

n_L4_2IG0_2D = N_t*sig_em*np.mean(np.mean(n_L4_2IG0_2D, axis=1), axis=0)
n_L4_4IG0_2D = N_t*sig_em*np.mean(np.mean(n_L4_4IG0_2D, axis=1), axis=0)
n_L4_8IG0_2D = N_t*sig_em*np.mean(np.mean(n_L4_8IG0_2D, axis=1), axis=0)


fig, ax = plt.subplots(1, 2, figsize=(20,10))


ax[0].plot([1, 2, 3, 4], [n_L1_IG0_1D*1e3, n_L2_IG0_1D*1e3, n_L3_IG0_1D*1e3, n_L4_IG0_1D*1e3], color=tableau20[0], label=r'$I_{G0}=1\times 10^{11}\/Wm^{-2}$')
ax[0].plot([1, 2, 3, 4], [n_L1_2IG0_1D*1e3, n_L2_2IG0_1D*1e3, n_L3_2IG0_1D*1e3, n_L4_2IG0_1D*1e3], color=tableau20[6], label=r'$I_{G0}=2\times 10^{11}\/Wm^{-2}$')
ax[0].plot([1, 2, 3, 4], [n_L1_4IG0_1D*1e3, n_L2_4IG0_1D*1e3, n_L3_4IG0_1D*1e3, n_L4_4IG0_1D*1e3], color=tableau20[2], label=r'$I_{G0}=4\times 10^{11}\/Wm^{-2}$')
ax[0].plot([1, 2, 3, 4], [n_L1_8IG0_1D*1e3, n_L2_8IG0_1D*1e3, n_L3_8IG0_1D*1e3, n_L4_4IG0_1D*1e3], color=tableau20[4], label=r'$I_{G0}=8\times 10^{11}\/Wm^{-2}$')
ax[0].set_xlim((0.9, 4.1))
ax[0].set_ylim((0, 0.24))
#ax[0].set_xlabel('Medium length (mm)', fontsize=16)
ax[0].set_ylabel('Average gain coefficient ($10^{-3}cm^{-1}$)', fontsize=16)
ax[0].set_yticks(np.arange(0.0, 0.25, 0.02))
lg = ax[0].legend()
lg.draw_frame(False)
ax[0].xaxis.set_ticks_position('bottom')
ax[0].yaxis.set_ticks_position('left')

ax[1].plot([2, 3, 4], [n_L2_2IG0_2D*1e3, n_L3_2IG0_2D*1e3, n_L4_2IG0_2D*1e3], color=tableau20[0], label=r'$I_{G0}=2\times 10^{10}\/Wm^{-2}$')
ax[1].plot([2, 3, 4], [n_L2_4IG0_2D*1e3, n_L3_4IG0_2D*1e3, n_L4_4IG0_2D*1e3], color=tableau20[6], label=r'$I_{G0}=4\times 10^{10}\/Wm^{-2}$')
ax[1].plot([2, 3, 4], [n_L2_8IG0_2D*1e3, n_L3_8IG0_2D*1e3, n_L4_8IG0_2D*1e3], color=tableau20[2], label=r'$I_{G0}=8\times 10^{10}\/Wm^{-2}$')
ax[1].set_xlim((1.9, 4.1))
ax[1].set_ylim((0, 0.12))
ax[1].set_yticks(np.arange(0.0, 0.13, 0.01))
#ax[1].set_xlabel('Medium length (mm)', fontsize=16)
#ax[1].set_ylabel('Average gain coefficient ($10^{-3}cm^{-1}$)', fontsize=16)
lg = ax[1].legend()
lg.draw_frame(False)
ax[1].xaxis.set_ticks_position('bottom')
ax[1].yaxis.set_ticks_position('left')

fig.text(0.5, 0.04, 'Medium dimensions (mm)', ha='center', va='center',fontsize=16)

fig.savefig('Averagegain.png')
#plt.show()