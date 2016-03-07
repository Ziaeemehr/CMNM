program dft_He
use He_ionization_finite_element
implicit none
real(8),allocatable :: rrad(:), hu(:), su(:),rho(:),a(:),b(:),vca(:),vcb(:),da(:)
real(8),allocatable :: urad(:),g(:),hhp(:,:),ssp(:,:),p(:),v1(:),v2(:),VH(:),w(:)
real(8) :: znuc,ener,sumn,gsum,alpha,gmax,tol,gamm,gnorm,pi,EH,EC,EXC,EXCtrap,eta
complex(8),allocatable::x(:)
integer :: lang,nrad,i
integer(kind=8)::it
!*****************************************************************************
nrad=10000
allocate(rrad(0:nrad),urad(0:nrad),hu(0:nrad),su(0:nrad),g(0:nrad),vca(0:nrad),vcb(0:nrad))
allocate(ssp(2,0:nrad),hhp(2,0:nrad),x(0:nrad),p(0:nrad),w(0:nrad),da(0:nrad))
allocate(v1(0:nrad),v2(0:nrad),rho(0:nrad),a(0:nrad),b(0:nrad),VH(0:nrad))
!*********************************
ener=0.d0; lang=0; hu=0.d0; su=0.d0; g=0.d0; it=0
gamm=1.d-1
znuc=2.d0
tol=1.d-7
!*********************************
call radgrid(nrad,rrad)
do i=0,nrad
      if (i==0) then 
         w(i)=rrad(i+1)/2.d0
      elseif (i==nrad) then
         w(i)=(rrad(i)-rrad(i-1))/2.d0
      else
         w(i)=(rrad(i+1)-rrad(i-1))/2.d0
      endif
      urad(i)=exp(-rrad(i))
enddo
!*************************************************************************
do while (1==1)
   it=it+1
   call overlap(nrad,rrad,urad,sumn,su)
   urad=urad/sqrt(sumn)
   call overlap(nrad,rrad,urad,sumn,su)
   call energr(nrad,lang,znuc,rrad,urad,ener,hu)
   call exhar(nrad,rrad,urad,VH,vca,vcb,w)
   da(0:nrad)=(VH(0:nrad)+vca(0:nrad))*rrad(0:nrad)**2*w(0:nrad)*urad(0:nrad)
   g(0:nrad)=hu(0:nrad)+da(0:nrad)-(ener+dot_product(da,urad))*su(0:nrad)
!   g(0:nrad)=hu(0:nrad)+da(0:nrad)-(ener*su(0:nrad))
   call overlap(nrad,rrad,g,gnorm,su)
   if (sqrt(gnorm) < tol) exit
   call crtssp(nrad,rrad,ssp)
   call crthhp(nrad,lang,znuc,rrad,hhp)
   x=cmplx(g,0.d0,8)
   hhp(1,1:nrad)=hhp(1,1:nrad)+(VH(1:nrad)+vca(1:nrad))*rrad(1:nrad)**2*w(1:nrad)
   call ctridag(nrad+1,hhp,ssp,ener,gamm,x)
   p=real(x,8)
   urad=urad-p
   write(*,'(1x,i5,4x,i7,es24.15,es25.15)') nrad,it,ener,gnorm
enddo
!write(*,'(1x,i5,1x,i7,es24.15,es25.15)') nrad,it,ener,gnorm
write(*,'(2x,a4,2x,a7,9x,a10,13x,a8)')'nrad','iter','energy','gmax'

contains
!*************************************************************************
subroutine exhar(nrad,rrad,urad,VH,vca,vcb,w)
implicit none
integer :: nrad,i
real(8) :: rrad(0:nrad), rho(0:nrad),a(0:nrad),b(0:nrad),vca(0:nrad),vcb(0:nrad)
real(8) :: urad(0:nrad),v1(0:nrad),v2(0:nrad),VH(0:nrad),w(0:nrad)
real(8) :: pi,EH,EC,EXC,EXCtrap,eta,up,down
!*************************************************************************
eta=0.d0
pi = 4*atan(1.d0)
a=0.d0; b=0.d0
VH=0.d0; rho=0.d0
EC=0.d0; EXC=0.d0; EXCtrap=0.d0;
do i=0,nrad
   !rho(i)=1.d0/(4*pi)*urad(i)**2
   rho(i)=2.d0/(4*pi)*urad(i)**2
enddo
do i=0,nrad
      a(i) =4*pi*rrad(i)**2*rho(i)*w(i)
      v1(i)=sum(a(0:i))
      b(i) =4*pi*rrad(i)*rho(i)*w(i)
      v2(i)=sum(b(0:i))
enddo
VH(0)=v2(nrad)
EH=0.d0
do i=1,nrad
VH(i)=1.d0/rrad(i)*v1(i)+(v2(nrad)-v2(i))
enddo
do i=1,nrad
   EH=EH+VH(i)*rho(i)*4*pi*rrad(i)**2*w(i)
   !if (rho(i)>1.d-20) then
      call LSD_PADE(rho(i),eta,EC,up,down)
      vca(i)=up; vcb(i)=down
   !else
   !   vca(i)=0.d0; vcb(i)=0.d0; EC=0.d0
   !endif
   EXC=EXC+EC*rho(i)*w(i)*4*pi*rrad(i)**2
enddo
EH=EH*5.d-1
!write(*,*) 'nrad    = ', nrad
!write(*,*) 'EH      = ', EH
!write(*,*) 'EXC     = ', EXC
!write(*,*) 'EH+EXC  =', EH+EXC

end subroutine exhar

!******************************************************************************
end program dft_He

