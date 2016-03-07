PROGRAM finite_element 
USE He_ionization_finite_element
     IMPLICIT NONE
     INTEGER::allocstat,nrad,i,j,k
     REAL(8)::ener,sumn,gmax,normg,alpha,ener_hrtre,ener_xc,pi,eta,c,znuc,lang,edfthyd,edfthelp,edfthel
     REAL(8)::enerhyd,ener_hrtrehyd,ener_xchyd,enerhelp,ener_hrtrehelp,ener_xchelp,enerhel,ener_hrtrehel,ener_xchel
     REAL(8),ALLOCATABLE::rrad(:),urad(:),vh(:),rho(:)
     WRITE(*,'(A)',ADVANCE='NO')'Enter nrad:'
     READ(*,*)nrad
     ALLOCATE(rrad(0:nrad))
     ALLOCATE(urad(0:nrad))
     ALLOCATE(rho(0:nrad))
     ALLOCATE(vh(0:nrad))
!------------------------------------------------------------------
     CALL radgrid(nrad,rrad)
     lang=0.d0
!------------------------------------------------------------------
     znuc=1.d0
     eta=1.d0 
     !CALL hyd_sch(nrad,lang,znuc,rrad,ener,urad)
     !WRITE(*,*)'  ener= ',ener
     !CALL hyd_sch_sd(nrad,lang,znuc,rrad,ener)
     CALL hyd_sch_psd(nrad,lang,znuc,rrad,ener,urad)
     !WRITE(*,*)'  ener= ',ener
     rho=urad**2
     CALL hrtre_ener(nrad,rrad,rho,ener_hrtre,vh)
     WRITE(*,*)'   Eh=  ',ener_hrtre  
     CALL exchngcorr_ener(nrad,eta,rrad,rho,ener_xc)
     WRITE(*,*)'   Exc= ',ener_xc
     WRITE(*,*)' Exc+Eh=',ener_xc+ener_hrtre
!----------------------------------------------------------------------------------------------------------------------
     znuc=1.d0
     eta=1.d0 
     c=1.d0
     CALL ks_orbit_ener(nrad,lang,znuc,c,eta,rrad,urad,i,enerhyd,ener_hrtrehyd,ener_xchyd)
     edfthyd= enerhyd+ener_hrtrehyd+ener_xchyd  
     WRITE(*,*)'i=',i
!----------------------------------------------------------------------------------------------------------------------
     znuc=2.d0
     eta=1.d0 
     c=1.d0
     CALL ks_orbit_ener(nrad,lang,znuc,c,eta,rrad,urad,j,enerhelp,ener_hrtrehelp,ener_xchelp)
     edfthelp= enerhelp+ener_hrtrehelp+ener_xchelp
     WRITE(*,*)'i=',j
!----------------------------------------------------------------------------------------------------------------------
     znuc=2
     eta=0.d0 
     c=2.d0
     CALL ks_orbit_ener(nrad,lang,znuc,c,eta,rrad,urad,k,enerhel,ener_hrtrehel,ener_xchel)
     edfthel=c*enerhel+ener_hrtrehel+ener_xchel
     WRITE(*,*)'i=',k
!----!------------------------------------------------------------------------------------------------------------------
     WRITE(*,*)
     WRITE(*,*)'************************** DFT Hydrogen **********************************'
     WRITE(*,'(1x,4(2x,A,ES23.16))'),'ener=',enerhyd,'Eh=',ener_hrtrehyd,'Exc=',ener_xchyd,'Etot=',edfthyd
     WRITE(*,*)
     WRITE(*,*)'*************************** DFT Helium1+ *********************************'
     WRITE(*,'(1x,4(2x,A,ES23.16))'),'ener=',enerhelp,'Eh=',ener_hrtrehelp,'Exc=',ener_xchelp,'Etot=',edfthelp
     WRITE(*,*)
     WRITE(*,*)'*************************** DFT Helium *********************************'
     WRITE(*,'(1x,4(2x,A,ES23.16))'),'ener=',enerhel,'Eh=',ener_hrtrehel,'Exc=',ener_xchel,'Etot=',edfthel
     WRITE(*,*)
     WRITE(*,*)'E ionization=',edfthelp-edfthel
END PROGRAM finite_element 
!***********************************************************************************************************************
!***********************************************************************************************************************
SUBROUTINE hyd_sch(nrad,lang,znuc,rrad,ener,urad)
USE He_ionization_finite_element
      IMPLICIT NONE
      INTEGER,INTENT(IN)::nrad
      REAL(8),INTENT(IN)::rrad(0:nrad),znuc,lang
      REAL(8),INTENT(OUT)::ener,urad(0:nrad)
      INTEGER::i
      REAL(8)::sumn,gmax,normg,pi
      REAL(8)::g(0:nrad),hgrad(0:nrad),sgrad(0:nrad)
      DO i=0,nrad
         urad(i)=2*(znuc**(3.d0/2.d0))*EXP(-znuc*rrad(i))
      END DO
      CALL overlap(nrad,rrad,urad,sumn,sgrad)
      urad=urad/SQRT(sumn)
      CALL overlap(nrad,rrad,urad,sumn,sgrad)
      CALL energr(nrad,lang,znuc,rrad,urad,ener,hgrad)
      g(0:nrad)=hgrad(0:nrad)-ener*sgrad(0:nrad)  
      gmax=MAXVAL(ABS(g))
      normg=SQRT(SUM(g**2))
      !WRITE(*,*)'**************************R wave****************************************'
      !WRITE(*,*)'normg=',normg
      !WRITE(*,*)'ener= ',ener
      !WRITE(*,*)'DE=',(-0.5d0-ener)
END SUBROUTINE hyd_sch
!***********************************************************************************************************************
SUBROUTINE hyd_sch_sd(nrad,lang,znuc,rrad,ener,urad)
USE He_ionization_finite_element
      IMPLICIT NONE
      INTEGER,INTENT(IN)::nrad
      REAL(8),INTENT(IN)::lang,znuc,rrad(0:nrad)
      REAL(8),INTENT(OUT)::ener,urad(0:nrad)
      INTEGER::i,j
      REAL(8)::alpha,sumn,normg
      REAL(8)::g(0:nrad),hgrad(0:nrad),sgrad(0:nrad),gu(0:nrad)
      DO i=0,nrad
         urad(i)=EXP(-(rrad(i)**2))
      END DO
      CALL overlap(nrad,rrad,urad,sumn,sgrad)
      urad=urad/SQRT(sumn)
      alpha=2.d-5
      j=0
      DO !i=1,100
         CALL energr(nrad,lang,znuc,rrad,urad,ener,hgrad)
         CALL overlap(nrad,rrad,urad,sumn,sgrad)
         g(0:nrad)=hgrad(0:nrad)-ener*sgrad(0:nrad)  
         CALL overlap(nrad,rrad,g,normg,gu)
         WRITE(*,*)j,normg, ener
         IF(normg<3.d-3) EXIT
         urad(0:nrad)=(urad(0:nrad)-alpha*g(0:nrad))
         CALL overlap(nrad,rrad,urad,sumn,sgrad)
         urad=urad/SQRT(sumn)
         j=j+1
      END DO
      WRITE(*,*)'**************************SD****************************************'
      WRITE(*,*)'normg=',normg
      WRITE(*,*)'ener= ',ener
      WRITE(*,*)' j=',j
END SUBROUTINE hyd_sch_sd
!************************************************************************************************************************
SUBROUTINE hyd_sch_psd(nrad,lang,znuc,rrad,ener,urad)
USE He_ionization_finite_element
      IMPLICIT NONE
      INTEGER,INTENT(IN)::nrad
      REAL(8),INTENT(IN)::lang,znuc,rrad(0:nrad)
      REAL(8),INTENT(OUT)::ener
      REAL(8),INTENT(OUT)::urad(0:nrad)
      INTEGER::i,j
      REAL(8)::shtr,shti,gama,t,sumn,gmax,normg,alpha
      REAL(8)::g(0:nrad),hgrad(0:nrad),sgrad(0:nrad),ssp(2,0:nrad),hhp(2,0:nrad),p(0:nrad)
      COMPLEX(8)::x(0:nrad)
      gama=0.4d0; t=1.d0
      DO i=0,nrad
         urad(i)=EXP(-(rrad(i)**2))
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
         x=CMPLX(g,0.d0,8); shtr=ener; shti=gama
         CALL ctridag(nrad+1,hhp,ssp,shtr,shti,x)
         p(0:nrad)=REAL(x(0:nrad),8)
         urad(0:nrad)=(urad(0:nrad)-t*p(0:nrad))
         CALL overlap(nrad,rrad,urad,sumn,sgrad)
         urad=urad/SQRT(sumn)
         j=j+1
         normg=SQRT(SUM(g**2))
         WRITE(*,*)'          normg=',normg, ener
         IF(normg<1.d-11) EXIT
      END DO
      !WRITE(*,*)'***************************PSD*****************************************'
      !WRITE(*,*)'normg=',normg
      !WRITE(*,*)'ener= ',ener
      WRITE(*,*)' j=',j
END SUBROUTINE hyd_sch_psd
!************************************************************************************************************************
SUBROUTINE hrtre_ener(nrad,rrad,rho,ener_hrtre,vh)
USE He_ionization_finite_element
      IMPLICIT NONE
      INTEGER,INTENT(IN)::nrad
      REAL(8),INTENT(IN)::rrad(0:nrad),rho(0:nrad)
      REAL(8),INTENT(OUT)::ener_hrtre,vh(0:nrad)
      INTEGER::i,j,open_status
      REAL(8)::intgrl,pi,int_rho_dr,sumn,ener_hrtreprm,rhor1,rhor2,intgrl1,intgrl2
      REAL(8)::sgrad(0:nrad),vprm(0:nrad),vh1(0:nrad),vprm1(0:nrad),urad(0:nrad)
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
      write(*,*) "EH"
      write(*,*) ener_hrtre
      
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
!************************************************************************************************************************
SUBROUTINE exchngcorr_ener(nrad,eta,rrad,rho,ener_xc)
      IMPLICIT NONE
      INTEGER,INTENT(IN)::nrad
      REAL(8),INTENT(IN)::rrad(0:nrad),rho(0:nrad),eta
      REAL(8),INTENT(OUT)::ener_xc
      INTEGER::j,i
      REAL(8)::pi
      REAL(8)::ec(0:nrad),vca(0:nrad),vcb(0:nrad)
      pi=4.d0*ATAN(1.d0)
      ec=0.d0
      DO i=0,nrad
         IF(rho(i)>1.d-20) CALL lsd_pade((RHO(i)/(4.d0*pi)),ETA,ec(i),vca(i),vcb(i))
      END DO
      ener_xc=(rho(0))*(rrad(0)**2)*ec(0)*((rrad(1)-rrad(0))/2.d0)
      DO j=1,nrad-1
         ener_xc=ener_xc+(rho(j))*(rrad(j)**2)*ec(j)*((rrad(j+1)-rrad(j-1))/2.d0)
      END DO
      ener_xc=ener_xc+(rho(nrad))*(rrad(nrad)**2)*ec(nrad)*((rrad(nrad)-rrad(nrad-1))/2.d0)
END SUBROUTINE exchngcorr_ener
!************************************************************************************************************************
SUBROUTINE ks_orbit_ener(nrad,lang,znuc,c,eta,rrad,urad1,j,ener,ener_hrtre,ener_xc)
      IMPLICIT NONE
      INTEGER,INTENT(IN)::nrad
      REAL(8),INTENT(IN)::lang,znuc,c,eta,rrad(0:nrad),urad1(0:nrad)
      REAL(8),INTENT(OUT)::ener_xc,ener,ener_hrtre
      INTEGER,INTENT(OUT)::j
      INTEGER::i,k
      REAL(8)::shtr,shti,gama,t,sumn,normg,evks,enerhyd,pi
      REAL(8)::g(0:nrad),hu(0:nrad),su(0:nrad),ssp(2,0:nrad),hhp(2,0:nrad),hhphyd(2,0:nrad),p(0:nrad)
      REAL(8)::rho(0:nrad),ec(0:nrad),vca(0:nrad),vcb(0:nrad),vh(0:nrad),da(0:nrad),urad(0:nrad),gu(0:nrad)
      COMPLEX(8)::x(0:nrad)
      pi=4.d0*ATAN(1.d0); t=1.d0; gama=0.5d0 
      DO i=0,nrad
         urad(i)=EXP(-(rrad(i)**2))
      END DO
      !urad=urad1
      CALL overlap(nrad,rrad,urad,sumn,su)
      urad=urad/SQRT(sumn)
      CALL crtssp(nrad,rrad,ssp)
      CALL crthhp(nrad,lang,znuc,rrad,hhphyd)
      j=0
      DO !k=1,3 
!--------------------------------------------------------------------------------------------------------
         rho(0:nrad)=c*urad(0:nrad)**2; vca=0.d0; ec=0.d0
         DO i=0,nrad
            IF(rho(i)>1.d-20) CALL lsd_pade((RHO(i)/(4.d0*pi)),ETA,ec(i),vca(i),vcb(i))
         END DO
         CALL hrtre_ener(nrad,rrad,rho,ener_hrtre,vh)
         da(0)=(vh(0)+vca(0))*(rrad(0)**2)*(rrad(1)-rrad(0))/2.d0
         DO i=1,nrad-1
             da(i)=(vh(i)+vca(i))*(rrad(i)**2)*(rrad(i+1)-rrad(i-1))/2.d0
         END DO 
         da(nrad)=(vh(nrad)+vca(nrad))*(rrad(nrad)**2)*(rrad(nrad)-rrad(nrad-1))/2.d0
!--------------------------------------------------------------------------------------------------------
         CALL overlap(nrad,rrad,urad,sumn,su)
         CALL energr(nrad,lang,znuc,rrad,urad,ener,hu)
         evks=c*ener+DOT_PRODUCT(da(0:nrad),urad(0:nrad)**2)  
         g(0:nrad)=c*hu(0:nrad)+da(0:nrad)*urad(0:nrad)-evks*su(0:nrad)  
         
         x=CMPLX(g,0.d0,8); shtr=evks; shti=gama
         hhp(1,:)=c*hhphyd(1,:)+da(0:nrad) 
         hhp(2,:)=c*hhphyd(2,:)
         CALL ctridag(nrad+1,hhp,ssp,shtr,shti,x)
         p(0:nrad)=REAL(x(0:nrad),8)

         urad(0:nrad)=urad(0:nrad)-t*p(0:nrad)
         CALL overlap(nrad,rrad,urad,sumn,su)
         urad=urad/SQRT(sumn)

         j=j+1
         CALL overlap(nrad,rrad,g,normg,gu)
         CALL exchngcorr_ener(nrad,eta,rrad,rho,ener_xc)
         WRITE(*,'(1x,I6,3(2x,ES23.16),2(x,ES9.2))')j,normg,(c*ener+ener_hrtre+ener_xc),evks,znuc,c
         IF(normg<1.d-12) EXIT
      END DO
END SUBROUTINE ks_orbit_ener


