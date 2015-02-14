! This is a program to solve the Random Laser problem

program diffusion

use variables

implicit none

interface
   function matinv(inmatrix, m2, m3, p1, itfrac)
   real, dimension(:,:,:) :: inmatrix, m2, m3
   real p1
   integer itfrac
   real, dimension(size(inmatrix,1),size(inmatrix,2),size(inmatrix,3)) :: matinv
   end function matinv
end interface

interface
   function pops(inmatrix, ma, mg, mr, nx, ny, nz, par1, par2, par3, delt, cn, ig0, mua, mur, mug)
   real, dimension(:,:,:) :: inmatrix, ma, mg, mr
   real par1, par2, par3, delt, cn, ig0, mua, mug, mur
   integer nx, ny, nz
   real, dimension(size(inmatrix,1),size(inmatrix,2),size(inmatrix,3)) :: pops
   end function pops
end interface

interface
   function intens(ig0, ke, tg, taug, c, z, t)
   real ig0, ke, tg, taug, c, z, t
   real intens
   end function intens
end interface

integer i, j, k, l, tfrac, pstp

real sabs, sem, aval, mfp, ig0, ir0, ext0, timeg, timer, pulseg, pulser, cn
real refract, life, ntot, mug, mur, mua, h
real a1, a2, da
real t, ig, ir, dt, dx, dy, dz, modt
real, allocatable, dimension(:,:,:) :: ncop
real, allocatable, dimension(:) :: x, y, z

! variables involved in model

!_____________________________________________________________________________

open(unit = instream, file = 'laser.d', status = 'old' )

read(instream, *) tsteps
read(instream, *) xsteps
read(instream, *) ysteps
read(instream, *) zsteps
read(instream, *) pstp

tsteps = min(maxtsteps, tsteps)		! setting number of time steps
xsteps = min(maxpsteps, xsteps)		! size of grid in each dimension
ysteps = min(maxpsteps, ysteps)
zsteps = min(maxpsteps, zsteps)
xsteps = max(minpsteps, xsteps)
ysteps = max(minpsteps, ysteps)
zsteps = max(minpsteps, zsteps)

read(instream, *) refract
read(instream, *) sabs
read(instream, *) sem
read(instream, *) mfp
read(instream, *) life
read(instream, *) ntot
read(instream, *) ext0
read(instream, *) pulseg
read(instream, *) pulser
read(instream, *) timeg
read(instream, *) timer
read(instream, *) ig0
read(instream, *) mua
read(instream, *) mug
read(instream, *) mur
read(instream, *) h

close (instream)

allocate(wa(1:xsteps,1:ysteps,1:zsteps))	! allocating sizes to variables
allocate(wg(1:xsteps,1:ysteps,1:zsteps))	! from 1 to minimum
allocate(wr(1:xsteps,1:ysteps,1:zsteps))
allocate(n1(1:xsteps,1:ysteps,1:zsteps))
allocate(ncop(1:xsteps,1:ysteps,1:zsteps))
allocate(ia(1:xsteps,1:ysteps,1:zsteps))
allocate(iga(1:xsteps,1:ysteps,1:zsteps))
allocate(ira(1:xsteps,1:ysteps,1:zsteps))
allocate(x(1:xsteps))
allocate(y(1:ysteps))
allocate(z(1:zsteps))

mua = c0*h/mua
mug = c0*h/mug
mua = c0*h/mur
cn = c0 / refract
aval = 1.0 / life
dcoeff = cn * mfp / 3.0
dt = 5.0e-12
dz = 1.0e-4
dy = dz
dx = dz
a1 = cn * sabs
a2 = cn * sem
da = dcoeff * dt / dz / dz		! alpha in discrete version
da = da / 3.0	
ir0 = 0.0				! with modified time step
mfp = mfp

!ig0 = ig0 / 1.0e17

!write(*,*) 'Transport speed = ',cn
!write(*,*) 'Einstein-A = ',aval
!write(*,*) 'D = ',dcoeff
!write(*,*) 'Timestep = ',dt
!write(*,*) 'Distance step = ',dx
!write(*,*) 'Absorption = ',a1
!write(*,*) 'Emission = ',a2
!write(*,*) 'Alpha = ',da
!write(*,*) 'Intensity ig0 = ',ig0

! Assign initial values

t = 0.0
x(1) = 0.0
do j = 2, xsteps
   x(j) = x(j-1) + dx
enddo
y(1) = 0.0
do k = 2, ysteps
   y(k) = y(k-1) + dy
enddo
z(1) = 0.0
do l = 2, zsteps
   z(l) = z(l-1) + dz
enddo
do j = 1, xsteps
   do k = 1, ysteps
      do l = 1, zsteps
         n1(j,k,l) = 0.0
         wr(j,k,l) = 0.0
         wg(j,k,l) = 0.0
         wa(j,k,l) = 0.0
      enddo
   enddo
enddo

!write(*,*) 'Dimensions of box are ',x(xsteps), y(ysteps), z(zsteps)

open(unit = outstream, file = 'laser.dat', status = 'unknown')

do i = 1, tsteps
   write(*,*) 'Time is ',i
   do tfrac = 1, 3
      modt = t + real(tfrac-1) * dt / 3.0
      ncop = n1
      n1 = pops(n1, wa, wg, wr, xsteps, ysteps, zsteps, a1, a2, aval, dt/3.0, cn, ig0, mua, mur, mug)
!wa = wa/(ig0/cn)
!wg = wg/(ig0/cn) ! addition to C
!wr = wr/(ig0/cn) ! addition to C
      ia = ncop * aval
      wa = matinv(wa, ncop*a2*dt/3.0*ntot, ia*cn/ig0*mua*ntot*dt/3.0, da, tfrac)	! matrix inversion
      do l = 1, zsteps
         ig = intens(ig0, ext0, timeg, pulseg, cn, z(l), modt)
         iga(:,:,l) = ig
      enddo
      wg = matinv(wg, (ncop-1.0)*a1*ntot*dt/3.0, iga*cn/ig0*dt/3.0/mfp, da, tfrac)	! da is diffusion coeff
!if (i == 1 .and. tfrac == 3) then
!do j = 1, xsteps
!do k = 1, ysteps
!do l = 1, zsteps
!write(*,*) wg(j,k,l)
!enddo
!enddo
!enddo
!stop
!endif
      do l = 1, zsteps
         ir = intens(ir0, ext0, timer, pulser, cn, z(l), modt)
         ira(:,:,l) = ir
      enddo
      wr = matinv(wr, ncop*a2*ntot*dt/3.0, ira*cn*dt/3.0/ig0/mfp, da, tfrac)

   enddo

      ! wa is the amplified laser light.


   ! Write out the whole data block for the current timestep
if ( mod(i,pstp) == 0 ) then
     write(outstream, *) t
       do j = 1, xsteps
         do k = 1, ysteps
            do l = 1, zsteps
            write(outstream,*) x(j), y(k), z(l), wa(j,k,l), wg(j,k,l), n1(j,k,l)
            enddo
         enddo
      enddo
   endif

!write out data block for one point, all tsteps
!if ( mod(i,pstp) == 0 ) then
!write(outstream, *) t
!       j = 7
!       k = 3
!       l = 2
!            write(outstream, *) x(j), y(k), z(l), wa(j,k,l)
!write(outstream, *) t
!       j = 20
!       k = 20
!       l = 20
!           write(outstream, *) x(j), y(k), z(l), wa(j,k,l)
!write(outstream, *) t
!       j = 3
!       k = 3
!       l = 2
!            write(outstream, *) x(j), y(k), z(l), wa(j,k,l), wg(j,k,l), !n1(j,k,l)
!endif

   t = t + dt

enddo

deallocate(wa)	! allocating sizes to variables
deallocate(wg)	! from 1 to minimum
deallocate(wr)
deallocate(n1)
deallocate(ncop)
deallocate(ia)
deallocate(iga)
deallocate(ira)
deallocate(x)
deallocate(y)
deallocate(z)

close(outstream)

stop

end program diffusion
