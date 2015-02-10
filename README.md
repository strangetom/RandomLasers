# RandomLasers
Masters project modelling random lasers

###Shortcuts  
**[29 January](#29-january)**  
**[1 February](#1-february)**  
**[2 February](#2-february)**  
**[3 February](#3-february)**  
**[5 February](#5-february)**  
**[8 February](#8-february)**  
**[10 February](#10-february)**  

___________
# 10 February  
* The values of I_G0 and I_R0 in the intensity functions should be in terms of photon count, rather than energy units. This also means that the energy densities are actually photon densities (so the 4 coupled equations now have consistant units).  
* It look like Lagendijk uses kappa_e of about 1e4, judging by the results of the gain coefficient as a function of depth into the material. Fig 8 shows how the gain coefficient varies as a function of depth. This is after the pump pulse has ended (set as 1 pulsewidth after the peak of the pulse).  
![alt-text](https://github.com/strangetom/RandomLasers/blob/master/.graphs/10Feb/Gain%20comparison.png "Figure 8")  
* Fig 8 matches has the correct shape. 
  * for the single sided pump, the gain coefficient decays rapidly into the medium.  
  * for the double sided pump, the gain coefficent is roughly uniform. It does dip in the centre (slightly).  

  
# 8 February    
* Sparse matrices make a significant difference when the B matrix is large, so the create_B_matrix function has been modified to use sparse matrices (csr) if the are more the 100 position elements. When the matrix is small, they actually slow things down. See Figure 7.  
![alt-text](https://github.com/strangetom/RandomLasers/blob/master/.graphs/08Feb/Sparse.Matrix.Benchmark.png "Figure 7")  
* Figure 6 can also be found on plotly [here](https://plot.ly/~strangetom/32/trace-1/). This allows for some interactivity, but it is a major pain in the arse to set up.  
  
# 5 February  
* Major correction to the numerical integration method. The space steps need to be smaller than the transport mean free path, which means the time steps get even smaller. The code has been modified to save less data than before since we don't need (1e6 x 45) data points per variable to plot graphs.  
* The amplified spontaneous emmission response seems to be only decaying with a lifetime given by tau_e (as one would expect). We are expecting to see some peaks in the time evolution when we vary certain parameters, but we haven't observed this.  
* 3D graphs! As examples, Fig 5 shows pump energy density, Fig 6 shows pump intensity.   
![alt-text](https://github.com/strangetom/RandomLasers/blob/master/.graphs/05Feb/3D.Pump.Energy.png "Figure 5")  
![alt-text](https://github.com/strangetom/RandomLasers/blob/master/.graphs/05Feb/3D.Pump.Intensity.png "Figure 6")  

# 3 February  
* Crank-Nicolson implementation seems to be unstable and barfs after around 140 steps (see Fig 1).    
![alt-text](https://github.com/strangetom/RandomLasers/blob/master/.graphs/03Feb/03.02.15.instability.of.method.png "Figure 1")
* Using the method Lagendijk used in his paper (Physical Review E, 1996), we get much more sensible results.
![alt-text](https://github.com/strangetom/RandomLasers/blob/master/.graphs/03Feb/03.02.15.W_G.png "Figure 2")  
In this Figure 2, the vertical blue line is where the maxima of the pump (and probe) pulse occur.  
Figure 3 (below), shows how N_1 evolves.
![alt-text](https://github.com/strangetom/RandomLasers/blob/master/.graphs/03Feb/03.02.15.N_1.png "Figure 3")  
_(Note: all times are in seconds)_  
Figures 2 and 3 follow the same shape as similar figures in Lagendijk's paper, as does the graph for the probe (not shown).    
The graph for the amplified spontaneous emmission is a bit funky, it may be that it decays over a much longer time.  
* Using @jit in numba doesn't speed things up at all in it's current state. Other optimizations applied that massively increase speed (For 100,000 iterations, netbook down from 14 mins to 40 seconds).  
* 10 microsecond calculation run. Only the data for W_A and N_1 were saved (graphs only, see Figure 4, since 10 million iterations is 1gb file size). Expected exponential decay observed. Woop!.  
![alt-text](https://github.com/strangetom/RandomLasers/blob/master/.graphs/03Feb/W_A.10million.png "Figure 4")  

# 2 February  
* Need to identify the dominant terms in the equations and find out why the calculation runs away.  
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

