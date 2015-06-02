# RandomLasers  

#### Abstract  
Simulations of a random laser were written in Python using two methods. Diffusion equations for
the energy densities present in the gain medium were coupled to the population of excited atoms
and solved using the Crank-Nicolson method. This allowed the magnitude of the output flux of
the laser to be investigated for different dimensions and pump intensities It was demonstrated
that this method was able to successfully model a random laser, but yielded no information
about the spectral properties of the emission. Maxwell’s time dependent field equations were
coupled to the polarization and the population equations for a four level laser system and solved
using finite-difference time domain numerical methods. This was successfully demonstrated
for an inactive medium and the spectral properties of the long lived modes of the electric field
were obtained. This method was unsuccessful in modelling a random laser with an active gain
medium, however.
____  
This simulations in this project were written using Python 3.4.3.  
The required libraries are:  
* __numpy 1.9.2__  
* __scipy 0.15.1__  

For the field equation simulations an additional library is required:  
* __numba 0.18.1__  
