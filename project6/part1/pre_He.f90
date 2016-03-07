program pre_He
implicit none
real(8),allocatable :: rrad(:), hu(:), su(:)
real(8),allocatable :: urad(:),g(:),hhp(:,:),ssp(:,:),p(:)
real(8) :: znuc,ener,sumn,gsum,alpha,gmax,tol,gamm,gnorm
complex(8),allocatable::x(:)
integer :: lang,nrad,i
integer(kind=8)::it
integer :: clck_counts_beg, clck_counts_end, clck_rate
!*****************************************************************************

call system_clock ( clck_counts_beg, clck_rate )

nrad = 99000       !nwrk in ctridag is 100000 and nrad sould be < nwrk 
allocate(rrad(0:nrad),urad(0:nrad),hu(0:nrad),su(0:nrad),g(0:nrad))
allocate(ssp(2,0:nrad),hhp(2,0:nrad),x(0:nrad),p(0:nrad))
call radgrid(nrad,rrad)
ener = 0.d0
znuc = 1.d0
lang = 0
! gamma should not be zero to determinant not be singular(probably)
! gnorm should be small enough to show this < 1e-14
gamm = 0.d0
hu = 0.d0
su = 0.d0
g  = 0.d0
it = 0
tol= 1.d-7

do i = 0, nrad
  urad(i) = exp(-0.5d0 * rrad(i) * rrad(i))
enddo

do while (1==1)
   it = it + 1
   call overlap(nrad,rrad,urad,sumn,su)
   urad = urad / sqrt(sumn)
   call overlap(nrad,rrad,urad,sumn,su)
   call energr(nrad,lang,znuc,rrad,urad,ener,hu)
   g(0:nrad) = hu(0:nrad) - ener * su(0:nrad)
   call overlap(nrad,rrad,g,gnorm,su)
   if (sqrt(gnorm) < tol) exit
   call crtssp(nrad,rrad,ssp)
   call crthhp(nrad,lang,znuc,rrad,hhp)
   x = cmplx(g,0.d0,8)
   call ctridag(nrad + 1, hhp, ssp, ener, gamm, x) ! arrays are (0:nrad) which equals nrad+1
   p = real(x,8)
   urad = urad - p
   write(*,'(1x,i5,i7,3es24.15)') nrad,it,ener, -5.d-1 - ener , sqrt(gnorm)
enddo

call system_clock ( clck_counts_end, clck_rate )
write(*,'(1x,a4,2x,a7,9x,a10,10x,a8,13x,a8)')'nrad','iter','energy','dE','gmax'
print *, 'Elapsed time in second'
write (*, *)  (clck_counts_end - clck_counts_beg) / real (clck_rate)

end program pre_He

