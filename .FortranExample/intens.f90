! This program contains the functions for the intensities of the pump and probe 
! pulses

function intens(ig0, ke, tg, taug, c, z, t)

implicit none

integer nz

real pi
real ig, ig0, ke, tg, taug, c, z, t
real b, par4
real aloss

real intens

!_____________________________________________________________________________

! define pi
pi = 4.0 * atan(1.0)
! b= 4ln2
b = 4.0 * log(2.0)
! par4 = (square root of 4ln2/pi)
par4 = sqrt( b / pi )

! Absorption loss
aloss = exp ( -1.0 * ke * z )


! tg = (time when maximum pump pulse is incident on surface)
! tr = (time when maximum probe pulse is incident on surface)
! igo = (average intensity of pump pulse)
! iro = (average intensity of probe pulse)
! Tr = (pulse lengths full width at half max)
! Tg = (-----------------"------------------)

ig = ig0*par4*aloss
ig = ig * exp(-b*(t-tg-z/c)**2/taug**2)
!   ir = iro*par4exp(-b(t-tr-z/c)**2/Tr**2)

!write(*,*) ig0
!write(*,*) b
!write(*,*) par4
!write(*,*) tg, c, z, ke, taug, t, aloss

intens = ig

return

end function intens
