import numpy as np
from scipy import constants as spc
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize=(20, 10))
ax = fig.add_subplot(111, projection='3d')
plt.hold(True)

store = np.loadtxt("./Data/L=1/W_A.I=1e12.L=1.txt", delimiter=',')

l = 100.e-6
L = 0.001  # length of slap
dx = l / 2.  # space steps across slap
x = np.arange(-l, L + l + dx, dx)  # array of space steps

N = 1000000  # number of time steps
T = 5.e-8  # length of time
dt = T / N  # time steps
# time array in nanoseconds
time = (dt * 500 * np.arange(0, store.shape[0], 1)) * 1e9

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
c = spc.c/1.35
E_A = 6.63e-34*c/700e-9
I_G0 = 1.e12


for i in range(x.shape[0]):
    position = (i * dx * np.ones(store.shape[0])) * 1e3 - l*1e3  # position in microns
    ax.plot(xs=time, ys=I_G0/c/E_A*store[:, i], zs=position, zdir='y', color='r', label="1 mm sample thickness\n1x10$^{12}$ $Jm^{-3}$ pump intensity" if i ==0 else "")

ax.set_xlabel('Time (ns)')
ax.set_xlim(0, 50)
ax.set_ylabel('Position (mm)')
ax.set_ylim(0, 1)
ax.set_zlabel('Intensity')
ax.set_zlim(0, 1.2e23)
ax.set_title('Amp. Spont. Emission intensity profile')
plt.legend(loc=1)
#fig.savefig('Graph.png')

plt.show()
