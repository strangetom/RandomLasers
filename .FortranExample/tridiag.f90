! Solves a tridiagonal system of linear equations by LU decomposition
! translated to FORTRAN90 from Numerical Recipes by Press et al.
!_______________________________________________________________________________

subroutine tridiag(a, b, c, r, u)

implicit none

real, dimension(:) :: a	! The vector of non-zero coefficients below diagonal
real, dimension(:) :: b ! -----------------"-----------------  on     "
real, dimension(:) :: c ! -----------------"----------------- above   "
real, dimension(:) :: r ! The right-hand side vector
real, dimension(:) :: u ! The solution space

integer n		! The number of equations to be solved
integer j

real bet
real, allocatable, dimension(:) :: gam

logical trifail

!_______________________________________________________________________________

! Determine n, assuming that the r.h.s. is the master size.
n = size(r)
trifail = .FALSE.
if (size(a) .ne. n) then
   write(*,*) 'tridiag failed: below diagonal vector not of required size.'
   trifail = .TRUE.
endif
if (size(b) .ne. n) then
   write(*,*) 'tridiag failed: diagonal vector not of required size.'
   trifail = .TRUE.
endif
if (size(c) .ne. n) then
   write(*,*) 'tridiag failed: above diagonal vector not of required size.'
   trifail = .TRUE.
endif
if (size(u) .ne. n) then
   write(*,*) 'tridiag failed: solution space not of required size.'
   trifail = .TRUE.
endif
if (trifail) stop

!__________________________________________________________________

! On reaching here, we know the input vectors are OK for size.
if (b(1) .eq. 0.0) then
   write(*,*) 'tridiag requests a rewrite of equations: b(1)=0.0'
   stop
endif
allocate(gam(1:n))  		! allocate workspace array

!_________________________________________________________________

! It's now OK to have a go at the actual algorithm.
bet = b(1)
u(1) = r(1) / bet
do j = 2, n				! Decomp. and forward subs.
   gam(j) = c(j-1) / bet
   bet = b(j) - a(j)*gam(j)
   if (bet .eq. 0.0) then
      write(*,*) 'tridiag failed'	! Algorithm fails: zero pivot
      stop
   endif
   u(j) = (r(j)-a(j)*u(j-1))/bet
enddo 

do j = n-1, 1, -1			! Backsubstitution
   u(j) = u(j) - gam(j+1)*u(j+1)
enddo

! Solution space complete; release work array
deallocate(gam)

return

end subroutine tridiag
