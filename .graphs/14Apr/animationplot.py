import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

# Tableau20 definitions
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),  
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),  
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),  
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),  
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)] 
for i in range(len(tableau20)):  
    r, g, b = tableau20[i]  
    tableau20[i] = (r / 255., g / 255., b / 255.)

E_storage = np.loadtxt('./Data/E_storage_2000dt_1000dx.txt', delimiter=',')
medium = np.loadtxt('./Data/Medium_2000dt_1000dx.txt', delimiter=',')


fig, ax = plt.subplots()
def animate(i):
    ax.cla()
    line = ax.plot(E_storage[i*3],color=tableau20[6])
    ax.fill_between(np.arange(0,medium.shape[0],1),1,-1,where=medium,color='0.7')
    ax.plot([0, 0], [-1, 1], color='k', linewidth=3) # boundary at 0
    ax.plot([medium.shape[0], medium.shape[0]], [-1, 1], color='k', linewidth=3) # bounary at the end
    ax.set_ylim((-1.1, 1.1))
    ax.set_xlim((-10, medium.shape[0]+10))
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.set_xlabel('Position')
    ax.set_ylabel('Electric field')
    ax.set_title("Time: "+str(i*3))
    return line
anim = animation.FuncAnimation(fig, animate, interval=1, frames=int(E_storage.shape[0]/3))
anim.save('movie.mp4', fps=17, bitrate=-1, codec='libx264')
#plt.show()