! This header contains defintions of global variables

module variables

implicit none

integer, parameter :: maxtsteps = 5000
integer, parameter :: maxpsteps = 250
integer, parameter :: minpsteps = 4
integer, parameter :: instream = 30
integer, parameter :: outstream = 40

real, parameter :: c0 = 2.9979e8

integer tsteps
integer xsteps
integer ysteps
integer zsteps

real dcoeff
real, allocatable, dimension(:,:,:) :: wa, wg, wr, n1
real, allocatable, dimension(:,:,:) :: ia, iga, ira

end module variables
