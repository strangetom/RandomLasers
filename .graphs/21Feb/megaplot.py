import numpy as np
from matplotlib import pyplot as plt
from scipy import constants as spc

c = spc.c/1.35 # speed of light in medium
E_A = 6.63e-34*c/700e-9

#Load data
flux1x1 = 1.e11/(c*E_A*1e22)*np.loadtxt('./Data/L=1/Flux.I=1e11.L=1.txt',delimiter=',')
flux1x2 = 2.e11/(c*E_A*1e22)*np.loadtxt('./Data/L=1/Flux.I=2e11.L=1.txt',delimiter=',')
flux1x4 = 4.e11/(c*E_A*1e22)*np.loadtxt('./Data/L=1/Flux.I=4e11.L=1.txt',delimiter=',')

N1x1 = np.loadtxt('./Data/L=1/N_pop.I=1e11.L=1.txt',delimiter=',')
N1x2 = np.loadtxt('./Data/L=1/N_pop.I=2e11.L=1.txt',delimiter=',')
N1x4 = np.loadtxt('./Data/L=1/N_pop.I=4e11.L=1.txt',delimiter=',')

flux2x1 = 1.e11/(c*E_A*1e22)*np.loadtxt('./Data/L=2/Flux.I=1e11.L=2.txt',delimiter=',')
flux2x2 = 2.e11/(c*E_A*1e22)*np.loadtxt('./Data/L=2/Flux.I=2e11.L=2.txt',delimiter=',')
flux2x4 = 4.e11/(c*E_A*1e22)*np.loadtxt('./Data/L=2/Flux.I=4e11.L=2.txt',delimiter=',')

N2x1 = np.loadtxt('./Data/L=2/N_pop.I=1e11.L=2.txt',delimiter=',')
N2x2 = np.loadtxt('./Data/L=2/N_pop.I=2e11.L=2.txt',delimiter=',')
N2x4 = np.loadtxt('./Data/L=2/N_pop.I=4e11.L=2.txt',delimiter=',')

flux3x1 = 1.e11/(c*E_A*1e22)*np.loadtxt('./Data/L=3/Flux.I=1e11.L=3.txt',delimiter=',')
flux3x2 = 2.e11/(c*E_A*1e22)*np.loadtxt('./Data/L=3/Flux.I=2e11.L=3.txt',delimiter=',')
flux3x4 = 4.e11/(c*E_A*1e22)*np.loadtxt('./Data/L=3/Flux.I=4e11.L=3.txt',delimiter=',')

N3x1 = np.loadtxt('./Data/L=3/N_pop.I=1e11.L=3.txt',delimiter=',')
N3x2 = np.loadtxt('./Data/L=3/N_pop.I=2e11.L=3.txt',delimiter=',')
N3x4 = np.loadtxt('./Data/L=3/N_pop.I=4e11.L=3.txt',delimiter=',')

time = (5.e-13 * 50 * np.arange(0, flux1x4.shape[0], 1)) * 1e9

#Setup axes
fig, ax = plt.subplots(3,3,figsize=(12,12))
fig.subplots_adjust(hspace = .001, wspace = .001)

#Plot everything
ax[0,0].plot(time, flux1x4,'r')
ax[0,0].set_ylim((0,12))
ax[0,0].set_xticklabels(())
ax[0,0].set_yticklabels(np.arange(0,12+1,2))
ax2 = ax[0,0].twinx()
ax2.plot(time, N1x4[:,12],'b')
ax2.set_yticklabels(())
ax2.set_ylim(0,0.8)

ax[0,1].plot(time, flux2x4,'r')
ax[0,1].set_ylim((0,12))
ax[0,1].set_xticklabels(())
ax[0,1].set_yticklabels(np.arange(0,12+1,1))
ax[0,1].set_yticklabels(())
ax2 = ax[0,1].twinx()
ax2.plot(time, N2x4[:,22],'b')
ax2.set_yticklabels(())
ax2.set_ylim(0,0.8)

ax[0,2].plot(time, flux3x4,'r')
ax[0,2].set_ylim((0,12))
ax[0,2].set_xticklabels(())
ax[0,2].set_yticklabels(np.arange(0,12+1,1))
ax[0,2].set_yticklabels(())
ax2 = ax[0,2].twinx()
ax2.plot(time, N3x4[:,32],'b')
ax2.set_yticks(np.arange(0, 0.8+0.1, 0.1))
ax2.set_ylim(0,0.8)

ax[1,0].plot(time, flux1x2,'r')
ax[1,0].set_ylim((0,7))
ax[1,0].set_xticklabels(())
ax[1,0].set_yticklabels(np.arange(0,6+1,1))
ax2 = ax[1,0].twinx()
ax2.plot(time, N1x2[:,12],'b')
ax2.set_yticklabels(())
ax2.set_ylim(0,0.8)

ax[1,1].plot(time, flux2x2,'r')
ax[1,1].set_ylim((0,7))
ax[1,1].set_xticklabels(())
ax[1,1].set_yticklabels(np.arange(0,6+1,1))
ax[1,1].set_yticklabels(())
ax2 = ax[1,1].twinx()
ax2.plot(time, N2x2[:,10],'b')
ax2.set_yticklabels(())
ax2.set_ylim(0,0.8)

ax[1,2].plot(time, flux3x2,'r')
ax[1,2].set_ylim((0,7))
ax[1,2].set_xticklabels(())
ax[1,2].set_yticklabels(np.arange(0,6+1,1))
ax[1,2].set_yticklabels(())
ax2 = ax[1,2].twinx()
ax2.plot(time, N3x2[:,12],'b')
ax2.set_yticks(np.arange(0, 0.7+0.1, 0.1))
ax2.set_ylim(0,0.8)

ax[2,0].plot(time, flux1x1,'r')
ax[2,0].set_ylim((0,4))
ax[2,0].set_xticklabels([0,10,20,30,40])
ax[2,0].set_yticklabels(np.arange(0,3+1,0.5))
ax2 = ax[2,0].twinx()
ax2.plot(time, N1x1[:,12],'b')
ax2.set_yticklabels(())
ax2.set_ylim(0,0.6)

ax[2,1].plot(time, flux2x1,'r')
ax[2,1].set_ylim((0,4))
ax[2,1].set_xticklabels([0,10,20,30,40])
ax[2,1].set_yticklabels(np.arange(0,3+1,0.5))
ax[2,1].set_yticklabels(())
ax2 = ax[2,1].twinx()
ax2.plot(time, N2x1[:,12],'b')
ax2.set_yticklabels(())
ax2.set_ylim(0,0.6)

ax[2,2].plot(time, flux3x1,'r')
ax[2,2].set_ylim((0,4))
ax[2,2].set_xticklabels([0,10,20,30,40,50])
ax[2,2].set_yticklabels(np.arange(0,3+1,0.5))
ax[2,2].set_yticklabels(())
ax2 = ax[2,2].twinx()
ax2.plot(time, N3x1[:,12],'b')
ax2.set_yticks(np.arange(0, 0.5+0.1, 0.1))
ax2.set_ylim(0,0.6)

#Axis labels
fig.text(0.5, 0.05, 'Time (ns)', ha='center', va='center',fontsize=18)
fig.text(0.08, 0.5, 'Flux ($10^{22}s^{-1}$)', ha='center', va='center',rotation='vertical',fontsize=18)
fig.text(0.95, 0.5, 'Excitation level', ha='center', va='center',rotation=270,fontsize=18)

#Label each figure
fig.text(0.26, 0.325, r'$L = 1 mm$''\n'r'$I_{G0} = 1\times10^{11}Wm^{-2}$', bbox=dict(facecolor='none',edgecolor='black'))
fig.text(0.26, 0.59, r'$L = 1 mm$''\n'r'$I_{G0} = 2\times10^{11}Wm^{-2}$', bbox=dict(facecolor='none',edgecolor='black'))
fig.text(0.253, 0.855, r'$L = 1 mm$''\n'r'$I_{G0} = 4\times10^{11}Wm^{-2}$', bbox=dict(facecolor='none',edgecolor='black'))

fig.text(0.517, 0.325, r'$L = 2 mm$''\n'r'$I_{G0} = 1\times10^{11}Wm^{-2}$', bbox=dict(facecolor='none',edgecolor='black'))
fig.text(0.517, 0.59, r'$L = 2 mm$''\n'r'$I_{G0} = 2\times10^{11}Wm^{-2}$', bbox=dict(facecolor='none',edgecolor='black'))
fig.text(0.51, 0.855, r'$L = 2 mm$''\n'r'$I_{G0} = 4\times10^{11}Wm^{-2}$', bbox=dict(facecolor='none',edgecolor='black'))

fig.text(0.777, 0.325, r'$L = 3 mm$''\n'r'$I_{G0} = 1\times10^{11}Wm^{-2}$', bbox=dict(facecolor='none',edgecolor='black'))
fig.text(0.777, 0.59, r'$L = 3 mm$''\n'r'$I_{G0} = 2\times10^{11}Wm^{-2}$', bbox=dict(facecolor='none',edgecolor='black'))
fig.text(0.77, 0.855, r'$L = 3 mm$''\n'r'$I_{G0} = 4\times10^{11}Wm^{-2}$', bbox=dict(facecolor='none',edgecolor='black'))

fig.savefig('1Dsystem length and pump variation.png')