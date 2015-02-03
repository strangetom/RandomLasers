# RandomLasers
Masters project modelling random lasers

###Shortcuts  
**[29 January](#29-january)**  
**[1 February](#1-february)**  
**[2 February](#2-february)**  
**[3 February](#3-february)**

___________
# 3 February  
* Crank-Nicolson implementation seems to be unstable and barfs after around 140 steps (see Fig 1).    
![alt-text](https://github.com/strangetom/RandomLasers/blob/master/.graphs/03Feb/03.02.15.instability.of.method.png "Figure 1")
* Using the method Lagendijk used in his paper (Physical Review E, 1996), we get much more sensible results.
![alt-text](https://github.com/strangetom/RandomLasers/blob/master/.graphs/03Feb/03.02.15.W_G.png "Figure 2")  
In this Figure 2, the vertical blue line is where the maxima of the pump (and probe) pulse occur.  
Figure 3 (below), shows how N_1 evolves.
![alt-text](https://github.com/strangetom/RandomLasers/blob/master/.graphs/03Feb/03.02.15.N_1.png "Figure 3")  
_(Note: all times are in nanoseconds)_  
Figures 2 and 3 follow the same shape as similar figures in Lagendijk's paper, as does the graph for the probe (not shown).    
The graph for the amplified spontaneous emmission is a bit funky, it may be that it decays over a much longer time.   

# 2 February  
* Need to indentify the dominant terms in the equations and find out why the calculation runs away.  
* Numbapro academic license obtained from continuum.io so optimisation using @jit decorator can be implemented at some point.  



# 1 February  
Update  
* Boundary conditions fixed, makes no difference at this time.  
* IPython notebook examples makes sense, suspect problem with parameters of our calcalation.  
*  Lagendijk finite difference method isn't Crank Nicolson since the 2nd order space dereivate is independent of time. It maybe be worth looking at implementing their method on Tuesday (it would highlight any problems with my code.)  
* Boost and lapack are installed (on laptop and netbook). It's going to be a faff to translate the python code into C++.  

 

# 29 January  
Thoughts for before Tuesday  
* Go through IPython notebook example fully  
* Use finite difference method as detailed in paper by Wiersma & Lagendijk  
* Look into doing this in C++: Boost library, lapack etc.  
* Fix boundary conditions

