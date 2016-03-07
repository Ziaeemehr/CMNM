program dft_He
! Natral Helium
use He_ionization_finite_element
implicit none
real(8),allocatable :: rrad(:), hu(:), su(:),rho(:),a(:),b(:),vca(:),vcb(:),da(:)
real(8),allocatable :: urad(:),g(:),hhp(:,:),ssp(:,:),p(:),v1(:),v2(:),VH(:),w(:)
real(8) :: znuc,ener,sumn,gsum,alpha,gmax,tol,gamm,gnorm,pi,EH,EC,EXC,EXCtrap,eta
complex(8),allocatable::x(:)
integer :: lang,nrad,i
integer(kind=8)::it
!------------------------------------------------------------------
WRITE(*,'(A)',ADVANCE='NO')'Enter nrad:'
READ(*,*)nrad
allocate(rrad(0:nrad),urad(0:nrad),hu(0:nrad),su(0:nrad),g(0:nrad),vca(0:nrad),vcb(0:nrad))
allocate(ssp(2,0:nrad),hhp(2,0:nrad),x(0:nrad),p(0:nrad),w(0:nrad),da(0:nrad))
allocate(v1(0:nrad),v2(0:nrad),rho(0:nrad),a(0:nrad),b(0:nrad),VH(0:nrad))
!------------------------------------------------------------------
     CALL radgrid(nrad,rrad)
     lang = 0.d0
!------------------------------------------------------------------
     znuc = 1.d0
     eta  = 1.d0
     CALL hyd_sch_psd(nrad,lang,znuc,rrad,ener,urad)
     !WRITE(*,*)'  ener= ',ener
     rho = urad * urad
     CALL hrtre_ener(nrad,rrad,rho,ener_hrtre,vh)
     WRITE(*,*)'   Eh=  ',ener_hrtre  
     CALL exchngcorr_ener(nrad,eta,rrad,rho,ener_xc)
     WRITE(*,*)'   Exc= ',ener_xc
     WRITE(*,*)' Exc+Eh=',ener_xc+ener_hrtre
!------------------------------------------------------------------
     znuc=1.d0
     eta=1.d0 
     c=1.d0
     CALL ks_orbit_ener(nrad,lang,znuc,c,eta,rrad,urad,i,enerhyd,ener_hrtrehyd,ener_xchyd)
     edfthyd= enerhyd + ener_hrtrehyd + ener_xchyd  
     WRITE(*,*)'i=', i
!------------------------------------------------------------------
     znuc=2.d0
     eta=1.d0 
     c=1.d0
     CALL ks_orbit_ener(nrad,lang,znuc,c,eta,rrad,urad,j,enerhelp,ener_hrtrehelp,ener_xchelp)
     edfthelp= enerhelp + ener_hrtrehelp + ener_xchelp
     WRITE(*,*)'i=',j
!------------------------------------------------------------------
     znuc=2
     eta=0.d0 
     c=2.d0
     CALL ks_orbit_ener(nrad,lang,znuc,c,eta,rrad,urad,k,enerhel,ener_hrtrehel,ener_xchel)
     edfthel = c * enerhel + ener_hrtrehel + ener_xchel
     WRITE(*,*)'i=',k
!------------------------------------------------------------------
SUBROUTINE hyd_sch_psd(nrad,lang,znuc,rrad,ener,urad)
USE He_ionization_finite_element
      IMPLICIT NONE
      INTEGER,INTENT(IN)::nrad
      REAL(8),INTENT(IN)::lang,znuc,rrad(0:nrad)
      REAL(8),INTENT(OUT)::ener
      REAL(8),INTENT(OUT)::urad(0:nrad)
      INTEGER::i,j
      REAL(8)::gama,t,sumn,gmax,normg,alpha
      REAL(8)::g(0:nrad),hgrad(0:nrad),sgrad(0:nrad),ssp(2,0:nrad),hhp(2,0:nrad),p(0:nrad)
      COMPLEX(8)::x(0:nrad)
      gama=0.4d0; t=1.d0
      DO i=0,nrad
         urad(i) = EXP(-rrad(i)*rrad(i))
      END DO
      CALL overlap(nrad,rrad,urad,sumn,sgrad)
      urad=urad/SQRT(sumn)
      CALL crtssp(nrad,rrad,ssp)
      CALL crthhp(nrad,lang,znuc,rrad,hhp)
      j=0
      DO 
         CALL overlap(nrad,rrad,urad,sumn,sgrad)
         CALL energr(nrad,lang,znuc,rrad,urad,ener,hgrad)
         g(0:nrad)=hgrad(0:nrad)-ener*sgrad(0:nrad)  
         x=CMPLX(g,0.d0,8); ;
         CALL ctridag(nrad+1,hhp,ssp,ener,gama,x)
         p(0:nrad)=REAL(x(0:nrad),8)
         urad(0:nrad)=(urad(0:nrad)-t*p(0:nrad))
         CALL overlap(nrad,rrad,urad,sumn,sgrad)
         urad=urad/SQRT(sumn)
         j=j+1
         CALL overlap(nrad,rrad,g,normg,su)
         su = 0.
         WRITE(*,*)'          normg=',normg, ener
         IF (sqrt(normg) < tol) EXIT        
      END DO
      !WRITE(*,*)'***************************PSD*****************************************'
      !WRITE(*,*)'normg=',normg
      !WRITE(*,*)'ener= ',ener
      WRITE(*,*)' j=',j
END SUBROUTINE hyd_sch_psd
!------------------------------------------------------------------
SUBROUTINE hrtre_ener(nrad,rrad,rho,ener_hrtre,vh)
USE He_ionization_finite_element
      IMPLICIT NONE
      INTEGER,INTENT(IN)::nrad
      REAL(8),INTENT(IN)::rrad(0:nrad),rho(0:nrad)
      REAL(8),INTENT(OUT)::ener_hrtre,vh(0:nrad)
      INTEGER::i,j,open_status
      REAL(8)::intgrl,pi,int_rho_dr,sumn,rhor1,rhor2,intgrl1
      REAL(8)::sgrad(0:nrad),vprm(0:nrad),vprm1(0:nrad),urad(0:nrad)
      pi=4.d0*ATAN(1.d0)
!--------------------------------------------------------------------------------------------------------
      !vprm(0)=0.d0
      !DO i=1,nrad
      !   CALL overlap(i,rrad,urad,sumn,sgrad)
      !   vprm(i)=sumn/rrad(i) 
      !END Do
      !WRITE(*,*)vprm(nrad)
!--------------------------------------------------------------------------------------------------------
      vh(0)=0.d0
      intgrl=(rrad(1)**2)*(rho(1))*((rrad(1)-rrad(0))/2.d0)
      vh(1)=intgrl/rrad(1)
      DO i=2,nrad
         rhor1=(rho(i-1))*(rrad(i-1)**2)*((rrad(i)-rrad(i-2))/2.d0)
         rhor2=(rho(i-1))*(rrad(i-1)**2)*((rrad(i-1)-rrad(i-2))/2.d0)
         intgrl=intgrl+rhor1-rhor2+(rho(i))*(rrad(i)**2)*((rrad(i)-rrad(i-1))/2)
         vh(i)=(intgrl/rrad(i))
      END Do
!--------------------------------------------------------------------------------------------------------
      vprm1=0.d0 
      intgrl1=(rho(nrad))*rrad(nrad)*((rrad(nrad)-rrad(nrad-1))/2.d0)
      intgrl1=intgrl1+(rho(nrad-1))*rrad(nrad-1)*((rrad(nrad)-rrad(nrad-2))/2) 
      vprm1(nrad-1)=intgrl1
      DO i=nrad-2,1,-1
         rhor1=(rho(i+1))*(rrad(i+1))*((rrad(i+2)-rrad(i+1))/2.d0)
         rhor2=(rho(i+1))*(rrad(i+1))*((rrad(i+2)-rrad(i))/2.d0)
         intgrl1=intgrl1+rhor2-rhor1+(rho(i))*(rrad(i))*((rrad(i+1)-rrad(i))/2.d0)
         vprm1(i)=intgrl1
      END DO
      vprm1(0)=intgrl1
!--------------------------------------------------------------------------------------------------------
      IF(nrad<10001) THEN !du to time consuming of writing huge data 
          OPEN(UNIT=2,FILE='hartree_pot',STATUS='REPLACE',ACTION='WRITE',IOSTAT=open_status)
          WRITE(2,'(1x,ES23.16)')(vh(i),i=0,nrad)

          OPEN(UNIT=3,FILE='rrad',STATUS='REPLACE',ACTION='WRITE',IOSTAT=open_status)
          WRITE(3,'(1x,ES23.16)')(rrad(i),i=0,nrad)

          OPEN(UNIT=4,FILE='rho',STATUS='REPLACE',ACTION='WRITE',IOSTAT=open_status)
          WRITE(4,'(1x,ES23.16)')(rho(i),i=0,nrad)
      END IF
!--------------------------------------------------------------------------------------------------------
      vh=vh+vprm1
      ener_hrtre=(rho(0))*(rrad(0)**2)*vh(0)*((rrad(1)-rrad(0))/4.d0)
      DO j=1,nrad-1
         ener_hrtre=ener_hrtre+(rho(j))*(rrad(j)**2)*vh(j)*((rrad(j+1)-rrad(j-1))/4.d0)
      END DO
      ener_hrtre=ener_hrtre+(rho(nrad))*(rrad(nrad)**2)*vh(nrad)*((rrad(nrad)-rrad(nrad-1))/4.d0)
!--------------------------------------------------------------------------------------------------------
      int_rho_dr=(rho(0))*(rrad(0)**2)*((rrad(1)-rrad(0))/2.d0)
      DO j=1,nrad-1
         int_rho_dr=int_rho_dr+(rho(j))*(rrad(j)**2)*((rrad(j+1)-rrad(j-1))/2.d0)
      END DO
      int_rho_dr=int_rho_dr+(rho(nrad))*(rrad(nrad)**2)*((rrad(nrad)-rrad(nrad-1))/2.d0)
      !WRITE(*,*)'int_rho_dr=',int_rho_dr
      DO i=0,nrad
         urad(i)=SQRT(rho(i))
      END DO
      CALL overlap(nrad,rrad,urad,sumn,sgrad)
      !WRITE(*,*)'     usu=  ',DOT_PRODUCT(urad,sgrad)
END SUBROUTINE hrtre_ener 
























ener=0.d0; lang=0; hu=0.d0; su=0.d0; g=0.d0; it=0
gamm=1.d-1
znuc=2.d0
tol=1.d-7
hhp = 0.d0;
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
   it = it+1
   call overlap(nrad,rrad,urad,sumn,su)
   urad=urad/sqrt(sumn)
   call overlap(nrad,rrad,urad,sumn,su)
   call energr(nrad,lang,znuc,rrad,urad,ener,hu)
   ener=ener*2
   call exhar(nrad,rrad,urad,VH,vca,vcb,w)
   da(0:nrad) = (VH(0:nrad)+vca(0:nrad)+vcb(0:nrad))*rrad(0:nrad)**2 *w(0:nrad)*urad(0:nrad)
   g(0:nrad) = hu(0:nrad)+da(0:nrad)-(ener+dot_product(da,urad))*su(0:nrad)
   !g(0:nrad)=hu(0:nrad)+da(0:nrad)-ener*su(0:nrad)
   call overlap(nrad,rrad,g,gnorm,su)
   su = 0.
   if (sqrt(gnorm) < tol) exit
   call crtssp(nrad,rrad,ssp)
   call crthhp(nrad,lang,znuc,rrad,hhp)
   x=cmplx(g,0.d0,8)
   hhp(1,1:nrad)=hhp(1,1:nrad)+(VH(1:nrad)+vca(1:nrad)+vcb(0:nrad))*rrad(1:nrad)**2*w(1:nrad)
   call ctridag(nrad+1,hhp,ssp,ener,gamm,x)
   p=real(x,8)
   urad = urad-p
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
eta = 0.d0
pi = 4 * atan(1.d0)
a  = 0.d0
b  = 0.d0
VH = 0.d0
rho = 0.d0
EC  = 0.d0
EXC = 0.d0 
EXCtrap=0.d0

do i = 0,nrad
   rho(i) = 2.d0/(4.d0 * pi) * urad(i) * urad(i)
enddo

do i = 0,nrad
      a(i)  = 4.d0 * pi * rrad(i) * rrad(i) * rho(i) * w(i)
      v1(i) = sum(a(0:i))
      b(i)  = 4.d0 * pi * rrad(i) * rho(i) * w(i)
      v2(i) = sum(b(0:i))
enddo

VH(0)=v2(nrad)
EH=0.d0

do i = 1,nrad
   VH(i)= 1.d0/rrad(i) * v1(i) + ( v2(nrad) - v2(i) )
enddo

do i = 1, nrad
   EH = EH + VH(i) * rho(i) * 4.d0 * pi * rrad(i)*rrad(i) * w(i)
   call LSD_PADE(rho(i),eta,EC,vca,vcb)
   vca(i) = up; vcb(i) = down
   EXC = EXC + EC * rho(i) * w(i) * 4.d0 * pi * rrad(i)*rrad(i)
enddo
EH = EH * 5.d-1
end subroutine exhar
end program dft_He

