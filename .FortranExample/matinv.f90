! This function inverts the sparse input matix, 'inmatrix', and returns the
! inverse as a matrix of the same size and shape.

function matinv(inmatrix, m2, m3, p1, itfrac)

implicit none

interface
   subroutine tridiag(a, b, c, r, u)
   real, dimension(:) :: a, b, c, r, u
   end subroutine tridiag
end interface

integer itfrac

real, dimension(:,:,:) :: inmatrix, m2
real, dimension(:,:,:) :: m3
real p1
real, dimension(size(inmatrix,1),size(inmatrix,2),size(inmatrix,3)) :: matinv

integer nx, ny, nz, ncubed
integer nsq, ndash
integer j, k, l, m, mm1
integer out, mid, inn, count

real, allocatable, dimension(:) :: rhs, diag, below, above, soln

!__________________________________________________________________________

! enter real code below

! Set up the rmight-hand side vector
nx = size(inmatrix,1)
ny = size(inmatrix,2)
nz = size(inmatrix,3)
ncubed = nx * ny * nz
allocate(rhs(1:ncubed))
allocate(diag(1:ncubed))
allocate(below(1:ncubed))
allocate(above(1:ncubed))
allocate(soln(1:ncubed))
count = 0

select case (itfrac)
case(1)

do l = 1, nz
   do k = 1, ny
      do j = 1, nx
         count = count + 1
         rhs(count) = m3(j,k,l)			! Adds the variable-free r.h.s
         rhs(count) = rhs(count) + inmatrix(j,k,l)*m2(j,k,l)
         if (k < ny) rhs(count) = rhs(count) + p1*inmatrix(j,k+1,l)/2.0
         if (k > 1) rhs(count) = rhs(count) + p1*inmatrix(j,k-1,l)/2.0
         if (l < nz) rhs(count) = rhs(count) + p1*inmatrix(j,k,l+1)/2.0
         if (l > 1) rhs(count) = rhs(count) + p1*inmatrix(j,k,l-1)/2.0
         rhs(count) = rhs(count) + inmatrix(j,k,l)*(1.0-2.0*p1)
         diag(count) = (1.0 + p1)		! defines the diagonal coeff.
         below(count) = 0.0
         above(count) = 0.0
         if (j > 1) then        
            below(count) = -p1/2.0
         endif
         if (j < nx ) then           
            above(count) = -p1/2.0
         endif
         soln(count) = 0.0
      enddo
   enddo
enddo

case(2)

do j = 1, nx
   do l = 1, nz
      do k = 1, ny
         count = count + 1
         rhs(count) = m3(j,k,l)			! Adds the variable-free r.h.s
         rhs(count) = rhs(count) + inmatrix(j,k,l)*m2(j,k,l)
         if (l < nz) rhs(count) = rhs(count) + p1*inmatrix(j,k,l+1)/2.0
         if (l > 1) rhs(count) = rhs(count) + p1*inmatrix(j,k,l-1)/2.0
         if (j < nx) rhs(count) = rhs(count) + p1*inmatrix(j+1,k,l)/2.0
         if (j > 1) rhs(count) = rhs(count) + p1*inmatrix(j-1,k,l)/2.0
         rhs(count) = rhs(count) + inmatrix(j,k,l)*(1.0-2.0*p1)
         diag(count) = (1.0 + p1)		! defines the diagonal coeff.
         below(count) = 0.0
         above(count) = 0.0
         if (k > 1) then
            below(count) = -p1/2.0
         endif
         if (k < ny ) then
            above(count) = -p1/2.0
         endif
         soln(count) = 0.0
!     write(*,*) count, rhs(count), below(count), diag(count), above(count)
!write(*,*) j,k,l, 'converts to ',count
      enddo
   enddo
enddo

case(3)

do k = 1, ny
   do j = 1, nx
      do l = 1, nz
         count = count + 1
         rhs(count) = m3(j,k,l)			! Adds the variable-free r.h.s
         rhs(count) = rhs(count) + inmatrix(j,k,l)*m2(j,k,l)
         if (j < nx) rhs(count) = rhs(count) + p1*inmatrix(j+1,k,l)/2.0
         if (j > 1) rhs(count) = rhs(count) + p1*inmatrix(j-1,k,l)/2.0
         if (k < ny) rhs(count) = rhs(count) + p1*inmatrix(j,k+1,l)/2.0
         if (k > 1) rhs(count) = rhs(count) + p1*inmatrix(j,k-1,l)/2.0
         rhs(count) = rhs(count) + inmatrix(j,k,l)*(1.0-2.0*p1)
         diag(count) = (1.0 + p1)		! defines the diagonal coeff.
         below(count) = 0.0
         above(count) = 0.0
         if (l > 1) then
            below(count) = -p1/2.0
         endif
         if (l < nz ) then
            above(count) = -p1/2.0
         endif
         soln(count) = 0.0
      enddo
   enddo
enddo

end select

!do out = 1, n
!   do mid = 1, n
!      do inn = 1, n
!         select case (itfrac)
!            case(1)
!               j = inn
!               k = mid
!               l = out
!            case(2)
!               k = inn
!               l = mid
!               j = out
!            case(3)
!               l = inn
!               j = mid
!               k = out
!         end select
!         count = count + 1
!         rhs(count) = m3(j,k,l)		! Adds the variable-free r.h.s
!         rhs(count) = rhs(count) + inmatrix(j,k,l)*m2(j,k,l)
!         if (k < n) rhs(count) = rhs(count) + p1*inmatrix(j,k+1,l)/2.0
!         if (k > 1) rhs(count) = rhs(count) + p1*inmatrix(j,k-1,l)/2.0
!         if (l < n) rhs(count) = rhs(count) + p1*inmatrix(j,k,l+1)/2.0
!         if (l > 1) rhs(count) = rhs(count) + p1*inmatrix(j,k,l-1)/2.0
!         rhs(count) = rhs(count) + p1*inmatrix(j,k,l)*(1.0-2.0*p1)
!         diag(count) = (1.0 + p1)		! defines the diagonal coeff.
!         if (count > 1) then
!            below(count) = 0.0
!            if (inn .gt. 1) below(count) = -p1/2.0
!         endif
!         if (count < ncubed ) then
!            above(count) = 0.0
!            if (inn .lt. n) above(count) = -p1/2.0
!         endif
!         soln(count) = 0.0
!      enddo
!   enddo
!enddo

call tridiag(below, diag, above, rhs, soln)

!do j=1, count
!write(*,*) soln(j)
!enddo

!count = 0
select case(itfrac)
   case(1)
      nsq = ny * nx
      ndash = nx
   case(2)
      nsq = nz * ny
      ndash = ny
   case(3)
      nsq = nx * nz
      ndash = nz
end select
do m = 1, count
   mm1 = m - 1
   out = 1 + mm1 / nsq
   mid = m - (out-1)*nsq
   mm1 = mid - 1
   mid = 1 + mm1 / ndash
   inn = m - (out-1)*nsq - (mid-1)*ndash
   select case(itfrac)
      case(1)
         j = inn
         k = mid
         l = out
      case(2)
         k = inn
         l = mid
         j = out
      case(3)
         l = inn
         j = mid
         k = out
   end select
   matinv(j,k,l) = soln(m)
!if (itfrac == 2) write(*,*) 'Eqn. ',m,' --> ',j,k,l
!   write(*,*) j,k,l,matinv(j,k,l)
enddo

deallocate(soln)
deallocate(rhs)
deallocate(diag)
deallocate(below)
deallocate(above)

!do out = 1, n
!   do mid = 1, n
!      do inn = 1, n
!         select case (itfrac)
!            case(1)
!               j = inn
!               k = mid
!               l = out
!            case(2)
!               k = inn
!               l = mid
!               j = out
!            case(3)
!               l = inn
!               j = mid
!               k = out
!         end select
!         count = count + 1
!         matinv(j,k,l) = soln(count)
!      enddo
!   enddo
! enddo

return 

end function matinv
