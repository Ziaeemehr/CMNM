program part2
use He_ionization_finite_element
! part 2-1
! check normalization of the charge density
implicit none
integer,parameter :: nrad=10000
real(8) :: rrad(0:nrad), hu(0:nrad), su(0:nrad),w(0:nrad)
real(8) :: urad(0:nrad),g(0:nrad),hhp(2,nrad),ssp(2,nrad),p(0,nrad)
real(8) :: znuc,ener,sumn,gsum,alpha,gmax,tol,gamm,gnorm,integ,pi
complex(8) ::x(0,nrad)
integer :: lang,i
integer(kind=8)::it
integer :: clck_counts_beg, clck_counts_end, clck_rate

pi = 4*atan(1.d0)
call radgrid(nrad,rrad)
ener = 0.d0
znuc = 1.d0
lang = 0
hu = 0.d0
su = 0.d0
g  = 0.d0
it = 0
tol=1.d-7
integ = 0.d0
gamm  = 1.d-1
do i = 0, nrad
   urad(i) = exp(-0.5d0*rrad(i)*rrad(i))
enddo
   call overlap(nrad,rrad,urad,sumn,su)
   urad = urad/sqrt(sumn)
do i = 0,nrad
   if (i==0) then 
       w(i) = rrad(i+1)/2.d0
   elseif (i==nrad) then
       w(i) = ( rrad(nrad)-rrad(nrad-1) )/2.d0
   else
       w(i) = (rrad(i+1)-rrad(i-1))/2.d0
   endif
   integ = integ + urad(i)*urad(i) * rrad(i)*rrad(i) * w(i)
enddo

write(*,*) 'number of grid = ', nrad
write(*,*) 'Total charge   = ', integ

end program part2
