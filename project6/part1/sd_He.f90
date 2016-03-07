program sd_He
implicit none
real(8),allocatable :: rrad(:), hu(:), su(:)
real(8),allocatable :: urad(:),g(:)
real(8) :: znuc,ener,sumn,gsum,alpha,gnorm,tol
integer :: lang,nrad,i
integer(kind=8)::it
integer :: clck_counts_beg, clck_counts_end, clck_rate
!******************************************************************************
call system_clock ( clck_counts_beg, clck_rate )

nrad = 100
allocate(rrad(0:nrad),urad(0:nrad),hu(0:nrad),su(0:nrad),g(0:nrad))
call radgrid(nrad,rrad)

ener = 0.d0
znuc = 1.d0
lang = 0
hu = 0.d0
su = 0.d0
g = 0.d0
it=0
alpha=2.d-5
tol=1.d-3

do i=0,nrad
  urad(i) = exp(-rrad(i) * rrad(i))
enddo

do while (1==1)
   it = it+1
   call overlap(nrad,rrad,urad,sumn,su)
   urad = urad/sqrt(sumn)
   call overlap(nrad,rrad,urad,sumn,su)
   call energr(nrad,lang,znuc,rrad,urad,ener,hu)
   g(0:nrad) = hu(0:nrad) - ener * su(0:nrad)
   call overlap(nrad,rrad,g,gnorm,su)
   !gnorm=sqrt(sum(g**2))
   if (sqrt(gnorm) < tol) exit
   urad = urad - alpha * g
   if (mod(it,100000) == 0) then
        write(*,'(1x,i3,4x,i7,3es24.15)') nrad,it,ener,-5.d-1-ener ,sqrt(gnorm)
   endif
enddo
call system_clock ( clck_counts_end, clck_rate )
write(*,'(1x,a4,2x,a7,9x,a10,10x,a8,13x,a8)')'nrad','iter','energy','dE','gmax'
!write(*,'(1x,i3,4x,i7,3es24.15)') nrad,it,ener,-5.d-1-ener ,sqrt(gnorm)
print *, 'Elapsed time in second'
write (*, *)  (clck_counts_end - clck_counts_beg) / real (clck_rate)
end program sd_He
