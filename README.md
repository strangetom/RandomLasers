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
**[14 February](#14-february)**   
**[19 February](#19-february)**  
**[21 February](#21-february)**  
**[24 February](#24-february)**  
**[26 February](#26-february)**  
**[28 February](#28-february)**  
**[3 March](#3-march)**  
**[5 March](#5-march)**  
**[8 March](#8-march)**  
**[10 March](#10-march)**  
**[17 March](#17-march)**  
**[19 March](#19-march)**  
**[7 April](#7-april)**  
**[9 April](#9-april)**  
**[14 April](#14-april)**  
**[18 April](#18-april)**  
**[21 April](#21-april)**  


___________
# 21 April  
* Continuing with the inactive scattering medium: The medium was excited by a Gaussian pulse of arbitrary amplitude (either 100 or 1000) and left to evolve over 250,000 time steps. The electric field between 125,000 and 250,000 time steps was recorded, summed over distance to get the time evolution of the field. This was Fourier Transformed and the absolute value taken to get the spectrum, which is shown in Figure 25.  
![alt-text](https://github.com/strangetom/RandomLasers/blob/master/.graphs/21Apr/100%20pulse/Spectrum.png "Figure 25")  
* For this, the medium length was 5 microns.  


# 18 April  
* Some work on the diffusion simulations: The average gain for the 1D systems was plotted for each pump intensity to show the variation and the affect the length of the medium has, see Figure 24.  
![alt-text](https://github.com/strangetom/RandomLasers/blob/master/.graphs/18Apr/Averagegain1D.png "Figure 24")  
* The 2D one will also be plotted. Since that will probably give the same results-ish, it maybe worth plotting them side by side for a) ease of comparison, and b) to save space.  

# 14 April  
* Little progress has been made with the complete system. The polarisation seems to explode when a source of any type or magnitude is added.  
* Taking a step backwards, a more simple case has been looked at. This is simply solving Maxwell's equations in 1D in a strongly scattering medium. A single Gaussian source (width 100 timesteps) is launched at timestep 30 in what would be the gain medium. No gain or polarisation is present. The time evolution is shown in Figure 22.  
![alt-text](https://github.com/strangetom/RandomLasers/blob/master/.graphs/14Apr/movie_2000dt_1000dx.gif "Figure 22")  
* This does everything one would expect it to do.
 * There a slower velocity in the high refractive index parts (grey)
 * There's reflection and transmission at each interface
 * The PML boundary conditions do a fine job of stopping any reflections from the edges, to simulate an open system.  
* To this, gain needs to be introduced  
 * I'm thinking maybe pump the atoms with P<sub>r</sub>=1e7 (ish) at the timestep corresponding to the peak of the gaussian pulse only.  
* The gain needs to be linked to the electric field via the polarisation.  
* For display in the report, the above gif has been made made into a static graph with snapshots at differing times in Figure 23.  
![alt-text](https://github.com/strangetom/RandomLasers/blob/master/.graphs/14Apr/static_2000dt_1000dx.png "Figure 23") 

# 9 April  
* Well here's something interesting: with the population equations corrected, it's super obvious that E, H and P will never deviate from 0 unless there's some kind of pump to E and sometime. This is what Jiang and Soukoulis are talking about in there paper when they mention homogenously distributed sources to stimulate spontaneous emission.  
* Code has been added that generates the positions of around 25 sources. These are roughly uniformly distributed in the medium (but they can only occur in the gain medium). 
* Jiang and Soukoulis' paper says these sources should generate waves of a Lorentzian frequency distribution centred around w<sub>a</sub> with an amplitude dependent on N<sub>2</sub>.


# 7 April  
Time to do more work on the field equation stuff.
* The polarization equation has been updated so the 1<sup>st</sup> order time derivative uses P<sup>n+1</sup> and P<sup>n-1</sup> so it's consistent with the 2<sup>nd</sup> order time derivative.  
* Variable naming conventions have been cleaned up and clarified in the main for loop for time progression.  
* It doesn't yet work... Basically nothing but a load of 0s comes out when run for 300,000 iterations (10<sup>-12</sup>s). The cause of this may be the lack of emitters, however the Sebbah and Vanneste paper does not use emitters and still manages to get results. They pump the four level system uniformly over the whole system. The numbers they use are missing (e.g, the pump rate), it would be kinda handy to have those.  
More updates to come later today. Probably. 
* The pump should not happen continuously I think. It should happen for a short period of time only at the beginning of the simulation (maybe only for the first iteration, or maybe as a gaussian pulse of short, say 10<sup>-16</sup> s, duration). Needs investigating.  


# 19 March  
* Writing the code utilising for loops is a lot slower than using broadcasting of numpy arrays, however it does work without overflowing and it becomes simpler to add boundary conditions.  
* The @jit decorator from numba brings some tremendous performance increases, so yay.  
* The code doesn't seem to produce the expected results when its pumped uniformly with the P<sub>r</sub> term in the rate equations. It looks like it may be necessary to add emitters to the fields.  
* N2 and N1 never seem to change from 0 (the inital condition) -> need to investigate why.  


# 17 March  
* Initial work on the simulations has begun. There are a number of points that need to be clarified.  
 * Absorbing boundary conditions at the edges (Liao method?)  
 * 300 nm leads at each side - are these just inacive regions for the boundary conditions?  
 * Need a more robust way to only make the population equations evaluate at the correct points in the medium. At the moment a binary mask is used, but I'm not sure forcing the population to 0 in the scattering regions is the correct way to go about it.  
* The Jiang paper mentions adding sources of spontaneous emission spread homogenously throughout the medium. A later paper by Sebbah does similar work in 2D, but does not include these sources. Whether or not these are required and how to implement them is unclear.  


# 10 March
* Diffusion equation work is now over.
 * There is some further work that can be done on the laser spiking. Eqn. 16 in the Lagendijk paper gives the frequency of the spikes, which we can use to compare with our results.  
* **Field equation work has begun.** 
* Jiang and Soukoulis published a paper in 2000 that simulated random laser odes by coupling Maxwell's equations to the rate equations of electron occupying each energy level (4 level system). They solve this in 1D, and use the time variation of the electric field to get the lasing modes (by doing a Fourier Transform).  
* The equations themselves may be reasonably simple to implement in code (for 1D), they're mostly first order time derivatives or curls (There's one 2nd order time derivate, but meh). A simple forward Euler method can be used.  
* The paper states that sources must by introduced to simulate spontaneous emission but doesn't state how. For this we could add the code for a Lorentzian emitter at each space step, then randomly allow them to emit at each time step. This idea needs work.  


# 8 March
* No progress with the non-square geometry.  
* Looked at the gain profiles through the centre of the slab, after the pump pulse has ended (at t = 30 ns). Figures 21 and 22 shows the gain profile at different lengths for an incident pulse below threshold (4x10<sup>9</sup> Wm<sup>-2</sup>) and above threshold (4x10<sup>10</sup> Wm<sup>-2</sup>) respectively.  
![alt-text](https://github.com/strangetom/RandomLasers/blob/master/.graphs/08Mar/Gainprofile.4e9Wm-2.png "Figure 21")  
![alt-text](https://github.com/strangetom/RandomLasers/blob/master/.graphs/08Mar/Gainprofile.4e10Wm-2.png "Figure 22")    
* The small extinction coefficient means that plots look more similar to the double sided pump from the Lagendijk paper, this is a result of the gain in the centre of the medium being depleted faster due to the higher excitation level. 
*It may be worth comparing one of these plots with a similar one with a higher extinction coefficient to show how it affects the gain profile*  
  

# 5 March  
* Using the Crank-Nicolson method whilst varying the geometry seems to be very difficult. The stems from the need for the A matrix to contain differnt N_pop values for each column of the W matrix, which isn't mathematically possible without defining a new A matrix for each column. This take a stupid amount of time.  
* An alternative might be to create a matrix for N<sub>t</sub> and set outer rows of this to 0 in order to change the effective geometry of the sample. Updates incoming.  
* **Update**: Progress was made by restricting N<sub>t</sub> and modifying the matrices on the RHS and LHS to include zeros in the upper and lower corners when calculating the x derivatives. The idea was that this would restrict gain in these areas. This was successful but the resulting animations show an unusual squareness to the excitation level and energy densities (see [here](https://github.com/strangetom/RandomLasers/blob/master/.graphs/05Mar/N_pop.L%3D3.q%3D15.IG0%3D2e11.mp4) and [here](https://github.com/strangetom/RandomLasers/blob/master/.graphs/05Mar/W_A.L%3D3.q%3D15.IG0%3D2e11.mp4) respectively)  that scaled with the width of the available gain.  
Why the atoms near the x = 0 and x = L boundaries are not excitated I do not know.  
* An alternative method that would deal with the issue in the first bullet point of this update would involve remapping the matrix each to a column vector.
This would require different remappings when calculating x and y derivatives but functions have been written that achieve this (M is the original matrix).  
```python
X = M.reshape(M.shape[0]*M.shape[1],1)
Y = M.T.reshape(M.shape[0]*M.shape[1],1)

def ytox(V):
	return V.reshape(M.shape[1],M.shape[0]).T.reshape(M.shape[0]*M.shape[1],1)

def xtoy(V):
	return V.reshape(M.shape[0],M.shape[1]).T.reshape(M.shape[0]*M.shape[1],1)
```  
Thus, using this method, a large matrix could be constructed with the correct elements that would allow the system to be solved. In theory.  
  
# 3 March  
* Code updated to correctly calculate 2D diffusion. The transpositions are now done in the correct places, so in the first half step the y derivative is kept constant, and in the second half step the x derivative is kept constant. **This doesn't actually change the results at all.**   
* An animation showing the time evolution for the excitation level over 500 ns can be found [here](https://github.com/strangetom/RandomLasers/blob/master/.graphs/03Mar/Super%20long%20time/N_pop.mp4). The purpose of this is to show that the population decays with the natural lifetime after lasing action has ended.  
* Things to do next:
 * Gut the probe from the code. It isn't used and probably won't be.
 * Make the 2D geometry variable (i.e. rectangular not just square)
  * Doing this with the current matrix methods may be trickier. The system won't be tridiagonal as such. It may be easier to allocate the right hand side as a vector operation. 
  
# 28 February  
* The 2D system has been correct to ensure the pump pulse is only applied along the x-axis. This was done by making the I_G and I_R functions return matrices that were then transposed when all the other matrices were. This appears to have made no change to the original Figure 17 (the one below has been updated anyway).  
* Animations showing how the 2D system evolves have been created:
 * [Pump energy density](https://github.com/strangetom/RandomLasers/blob/master/.graphs/27Feb/W_G.mp4)  
![alt-text](https://github.com/strangetom/RandomLasers/blob/master/.graphs/27Feb/W_G.gif "Figure 18")    
 * [Excitation level](https://github.com/strangetom/RandomLasers/blob/master/.graphs/27Feb/N_pop.mp4)  
![alt-text](https://github.com/strangetom/RandomLasers/blob/master/.graphs/27Feb/N_pop.gif "Figure 19")   
 * [Amplified spontaneous emission](https://github.com/strangetom/RandomLasers/blob/master/.graphs/27Feb/W_A.mp4)  
![alt-text](https://github.com/strangetom/RandomLasers/blob/master/.graphs/27Feb/W_A.gif "Figure 20") 
  
# 26 February  
* Further corrections made to the 1D CN code, mostly involved with the probe pulse so this hasn't changed any of the figures yet.
* The 2D system code has been written and appears to work successfully. Figure 16 shown below shows the back scattered flux at the x=0 interface from a pump pulse of intensity 4x10<sup>10</sup> Wm<sup>-2</sup> and a thickness of 2 mm. Note that thickness of 1 mm do not seem to allow for laser action even with extremely large pump intensities (~10<sup>16</sup> Wm<sup>-2</sup>)
![alt-text](https://github.com/strangetom/RandomLasers/blob/master/.graphs/26Feb/2D%20flux%20%2B%20excitation.png "Figure 16")  
* A mega plot similar to the 1D one has been plotted, for L = 2, 3, 4 mm. See Figure 17 below:  
![alt-text](https://github.com/strangetom/RandomLasers/blob/master/.graphs/26Feb/2D%20system%20length%20and%20pump%20variation.png "Figure 17")  
* Note how a lower pump pulse intensity results in higher output intensities than the 1D system. This is due to the plotted flux being the sum of the outgoing flux at each point along the interface.  
* In it's current state, a 2 mm thick sample takes 5.5 minutes to run, a 4 mm thick sample takes almost 13 minutes. Some optimization is clearly needed. Using multiprocessing to run the W_'X'_next calculations simultaneously seems to actual slows things down.  
* To do over the weekend: See if its possible to plot the energy density at each spatial grid point, and animiate it to show the evolution over time. Do this for 3 mm, 2x10<sup>10</sup> Wm<sup>-2</sup>, since it'll probably give the clearest peaks. Will need to regenerate the data, [use this to save the data](http://docs.scipy.org/doc/numpy/reference/generated/numpy.save.html).   

# 24 February  
* Some corrections were made to the 1D CN code, and the code to include and implement the probe pulse was added in, but hasn't been used. The graphs from [21 February](#21-february) have been updated accordingly.  
* The theory behind solving a 2 dimensional set of equations has been looked into. See section 20.3.2 of Numerical Methods (3<sup>rd</sup> edition).  
* We can solve the 2 dimensional equations in a similar way to the tridiagonal matrix equation we get in 1 dimension. The difference is that we have 2 equations and we move in half steps, alternating between solving the x-derivatives and the y-derivatives.  


# 21 February  
* Crank-Nicolson method has been implemented. Still slightly unsure of boundary conditions, but the results seem to be as expected.  
* Since CN is always stable, we can reduce the number of time steps. However since some parameters are multiplied by the time step, this changes the results. See Figure 13 below:  
![alt-text](https://github.com/strangetom/RandomLasers/blob/master/.graphs/21Feb/Step%20size%20comparison.png "Figure 13")  
* We can also show how the pump intensity affects the outgoing flux. See Figure 14 below:  
![alt-text](https://github.com/strangetom/RandomLasers/blob/master/.graphs/21Feb/Pump%20intensity%20comparison.png "Figure 14")  
* And to finish off, a big plot showing how thickness and pump intensity affect the flux and excitation level. The plot in the bottom left corner has a thickness of 1 mm and a pump intensity of 4x10<sup>10</sup> Wm<sup>-2</sup>.  
 Going right along the plots, the thickness increases in 1 mm steps.  
 Going up along the plots, the pump intensity doubles.  
![alt-text](https://github.com/strangetom/RandomLasers/blob/master/.graphs/21Feb/1Dsystem%20length%20and%20pump%20variation.png "Figure 15")  

# 19 February
* The PDEs were made dimensionless by introducing introducing dimensionless variables. Variables were made dimensionless by dividing by the largest value with the dimension.  
* These were discretized and implemented in code using the forward Euler method. Crank-Nicolson has not been attmepted yet (it will be).  
* The extinction coefficent has been changed from 2x10<sup>4</sup> to 2x10<sup>-4</sup>, which gives more reasonable results (this needs confirming/referencing).  
* Figure 11 shows the flux obtained with 100 um mean free path and 1 mm thickness. The only every 500th data point was saved.   
![alt-text](https://github.com/strangetom/RandomLasers/blob/master/.graphs/19Feb/Flux.png "Figure 11")  
* Now in 3D!  
![alt-text](https://github.com/strangetom/RandomLasers/blob/master/.graphs/19Feb/W_A.3D.png "Figure 12")


# 14 February  
* Yay for working on valentines day!  
* Code does not do what it's meant to.  Attempted to implement Crank-Nicolson again but it's just as crap as before. Cannot get the system to lase.
* We have an example fortran code written using CN method and in 3 dimensions. It mostly makes sense, however there are a lot extra multiplying factors I don't yet understand the purpose of.  
* Giving up with life.  

# 10 February  
* The values of I<sub>G0</sub> and I<sub>R0</sub> in the intensity functions should be in terms of photon count, rather than energy units. This also means that the energy densities are actually photon densities (so the 4 coupled equations now have consistant units).  
* It look like Lagendijk uses K<sub>e</sub> of 1x10<sup>4</sup>, judging by the results of the gain coefficient as a function of depth into the material. Fig 8 shows how the gain coefficient varies as a function of depth. This is after the pump pulse has ended (set as 1 pulsewidth after the peak of the pulse).  
![alt-text](https://github.com/strangetom/RandomLasers/blob/master/.graphs/10Feb/Gain%20comparison.png "Figure 8")  
* Fig 8 matches has the correct shape. 
  * for the single sided pump, the gain coefficient decays rapidly into the medium.  
  * for the double sided pump, the gain coefficent is roughly uniform. It does dip in the centre (slightly).  
* The 3D plots [below](#5-february) have been recreated with the corrected intensities. This changes the shape of them quite a bit. Shown below are Fig 9 (single sided pump intensity) and Fig 10 (double sided pump intensity).  
![alt-text](https://github.com/strangetom/RandomLasers/blob/master/.graphs/10Feb/3D.Pump.Intensity.Single.png "Figure 9")
![alt-text](https://github.com/strangetom/RandomLasers/blob/master/.graphs/10Feb/3D.Pump.Intensity.Double.png "Figure 10")

  
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

