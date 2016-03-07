program check3

implicit none
real(8),allocatable :: rrad(:), hu(:), su(:), alaki(:)   !we don't need alaki vector
real(8),allocatable :: urad(:),g(:)
real(8) :: znuc,ener,sumn,gsum, gnorm
integer :: lang,nrad,i

nrad = 100

allocate(rrad(0:nrad),urad(0:nrad),hu(0:nrad),su(0:nrad),g(0:nrad),alaki(0:nrad))

call radgrid(nrad,rrad)
ener = 0.d0
znuc = 1.d0
lang = 0
hu   = 0.d0
su = 0.d0
g  = 0.d0
do i = 0, nrad
 urad(i) = 2*znuc**(3/2)*exp(-znuc*rrad(i))
enddo
! To normalize the urad once call the overlap
call overlap(nrad,rrad,urad,sumn,su)
! So the normalized urad is 
urad =  urad /sqrt(sumn)

call energr(nrad,lang,znuc,rrad,urad,ener,hu)

g(0:nrad) = hu(0:nrad) - ener * su(0:nrad)

write(*,*) ener + 5.d-1 

! gnorm is not simply g'g
!do i = 0,nrad
!    gnorm = gnorm + g(i)*g(i)
!enddo

! To calculate the gnorm we should write
call overlap(nrad, rrad, g, gnorm, alaki)
write(*,*) "gnorm is ", gnorm

end program check3
