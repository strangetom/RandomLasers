! This program does a simple forward time step of the population equation

function pops(inmatrix, ma, mg, mr, nx, ny, nz, par1, par2, par3, delt, cn, ig0, mua, mur, mug)

implicit none

integer nx, ny, nz

real par1, par2, par3, delt, cn, ig0, mua, mur, mug
real, dimension(:,:,:) :: inmatrix, ma, mg, mr
real, dimension(size(inmatrix,1),size(inmatrix,2),size(inmatrix,3)) :: pops

integer j, k, l

real term1, term2, term3

!_________________________________________________________________________

! par1 = (absorption cross-section) x (local speed of light)
! par2 = (emission   ------"------) x (--------"-----------)
! par3 = 1.0 / (upper state lifetime)

do j = 1, nx
   do k = 1, ny
      do l = 1, nz
         term1 = par1*ig0/mug/cn *  ( 1.0 - inmatrix(j,k,l) ) * mg(j,k,l)
         term2 = par2*ig0/mua/cn * inmatrix(j,k,l) * (mua*mr(j,k,l)/mur + ma(j,k,l))
         term3 = par3 * inmatrix(j,k,l)
!write(*,*) term1, term2, term3, delt
         pops(j,k,l) = inmatrix(j,k,l) + delt * (term1-term2-term3)
!         write(*,*) 'At grid point ',j,k,l
      enddo
   enddo
enddo

return

end function pops
