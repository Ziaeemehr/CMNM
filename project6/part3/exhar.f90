subroutine exhar(nrad,rrad,urad,VH,vca,vcb,w)
implicit none
integer :: nrad,i
real(8) :: rrad(0:nrad), rho(0:nrad),a(0:nrad),b(0:nrad),vca(0:nrad),vcb(0:nrad)
real(8) :: urad(0:nrad),v1(0:nrad),v2(0:nrad),VH(0:nrad),w(0:nrad)
real(8) :: pi,EH,EC,EXC,EXCtrap,eta,up,down
!*************************************************************************
pi = 4*atan(1.d0)
a=0.d0; b=0.d0
VH=0.d0; rho=0.d0
EC=0.d0; EXC=0.d0; EXCtrap=0.d0;
do i=0,nrad
   rho(i)=1.d0/(4*pi)*urad(i)**2
enddo
do i=0,nrad
      if (i==0) then 
         w(i)=rrad(i+1)/2.d0
      elseif (i==nrad) then
         w(i)=(rrad(i)-rrad(i-1))/2.d0
      else
         w(i)=(rrad(i+1)-rrad(i-1))/2.d0
      endif
      a(i) =4*pi*rrad(i)**2*rho(i)*w(i)
      v1(i)=sum(a(0:i))
      b(i) =4*pi*rrad(i)*rho(i)*w(i)
      v2(i)=sum(b(0:i))
enddo
VH(0)=v2(nrad)
EH=0.d0
eta=1.d0
do i=1,nrad
VH(i)=1.d0/rrad(i)*v1(i)+(v2(nrad)-v2(i))
enddo
do i=1,nrad
   EH=EH+VH(i)*rho(i)*4*pi*rrad(i)**2*w(i)
   if (rho(i)>1.d-20) then
      call LSD_PADE(rho(i),eta,EC,up,down)
      vca(i)=up; vcb(i)=down
   else
      vca=0.d0; vcb=0.d0; EC=0.d0
   endif
   EXC=EXC+EC*rho(i)*w(i)*4*pi*rrad(i)**2
enddo
EH=EH*5.d-1
write(*,*) 'nrad    = ', nrad
write(*,*) 'EH      = ', EH
write(*,*) 'EXC     = ', EXC
write(*,*) 'EH+EXC  =', EH+EXC

end subroutine exhar
