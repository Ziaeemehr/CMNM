!******************************************************************************************
subroutine radgrid(nrad,rrad)
    !returns a logarithmic radial grid
    implicit real(8) (a-h,o-z)
    dimension rrad(0:nrad)
    rad0=1.d-12
    radn=1.d+2
    alpha=log((radn+rad0)/rad0)/nrad
    do i=0,nrad
        rrad(i) = rad0*(exp(alpha*i)-1.d0)
    enddo
    return
end
!******************************************************************************************
subroutine energr(nrad,lang,znuc,rrad,urad,ener,ugrad)
    implicit real(8) (a-h,o-z)
    dimension rrad(0:nrad),urad(0:nrad),ugrad(0:nrad)
!    returns the energy ener and the ugrad  =  urad * Hamiltonian  for the radial 
!    wavefunction urad in the potential znuc/r
    z=-znuc
    hllp1=.5d0*lang*(lang+1)
    ener=0.d0
    ugrad(0)=0.d0
!    urad(nrad)=0.d0
    do 100,i=0,nrad-1
        ul=urad(i)
        ulp1=urad(i+1)
        t2 = ulp1-ul
        t3 = t2**2
        rl=rrad(i)
        rlp1=rrad(i+1)
        t6 = rlp1-rl
        t7 = t6**2
        t8 = 1.d0/t7
        t9 = rlp1**2
        t10 = t9**2
        t11 = rl**2
        t12 = t11**2
        t15 = t8*(t10-t12)
        t23 = 1.d0/t6
        t26 = ul-rl*t2*t23
        t27 = z*t26
        t28 = t2*t23
        t35 = t9*rlp1-t11*rl
        t38 = t26**2
        t40 = hllp1*t26
        t45 = t9-t11
        ener = ener +                            &
             0.25D0*z*t3*t15+0.33333333333333D0*         &
             (hllp1*t3*t8+0.5D0*t3*t8+2.D0*t27*t28)*t35+     &
             0.5D0*(z*t38+2.D0*t40*t28)*t45+hllp1*t38*t6 
        t52 = 0.5D0*z*t2*t15
        t55 = 2.D0*hllp1*t2*t8
        t56 = t2*t8
        t58 = rl*t23
        t59 = 1.D0+t58
        t64 = 2.D0*t27*t23
        t74 = 2.D0*t40*t23
        ugrad(i) = ugrad(i) + .5d0* (            &
            -t52+0.33333333333333D0*                 &
            (-t55-t56+2.D0*z*t59*t28-t64)*t35+           &
            0.5D0*(2.D0*t27*t59+2.D0*hllp1*t59*t28-t74)*t45+   &
            2.D0*t40*t6*t59 )
        ugrad(i+1) = .5d0* (                                &
        t52+0.33333333333333D0*(t55+t56-2.D0*z*rl*t56+t64)*t35+0.5d0*    &
        (-2.D0*t27*t58-2.D0*hllp1*rl*t56+t74)*t45-2.D0*t40*rl )
100    continue
    return
end
!******************************************************************************************
subroutine overlap(nrad,rrad,urad,sumn,ugrad)
    implicit real(8) (a-h,o-z)
    dimension rrad(0:nrad),urad(0:nrad),ugrad(0:nrad)
!    returns the overlap sumn=urad^T * S * urad and ugrad = S* urad for the radial 
!    wavefunction urad 
    sumn=0.d0
    ugrad(0)=0.d0
!    urad(nrad)=0.d0
    do 100,i=0,nrad-1
        ul=urad(i)
        ulp1=urad(i+1)
        t2 = ulp1-ul
        t3 = t2**2
        rl=rrad(i)
        rlp1=rrad(i+1)
        t5 = rlp1-rl
        t6 = t5**2
        t7 = 1.d0/t6
        t9 = rlp1**2
        t10 = t9**2
        t12 = rl**2
        t13 = t12**2
        t16 = t10*rlp1-t13*rl
        t20 = 1.d0/t5
        t23 = ul-rl*t2*t20
        t26 = t10-t13
        t27 = t20*t26
        t30 = t23**2
        t34 = t9*rlp1-t12*rl
        sumn = sumn + 0.2D0*t3*t7*t16+0.5D0*t23*t2*t27+ 0.33333333333333D0*t30*t34
        t39 = 0.2D0*t2*t7*t16
        t40 = rl*t20
        t41 = 1.D0+t40
        t47 = 0.25D0*t23*t20*t26
        t48 = t23*t34
        ugrad(i) = ugrad(i) -t39+0.25D0*t41*t2*t27-t47+ 0.33333333333333D0*t48*t41
        ugrad(i+1) = t39-0.25D0*rl*t7*t2*t26+t47- 0.33333333333333D0*t48*t40
100    continue
    return
end
!******************************************************************************************
subroutine crtssp(nrad,rrad,ssp)
    implicit real(8) (a-h,o-z)
    dimension rrad(0:nrad),ssp(2,0:nrad)
    !creates tridiagonal overlap matrix
    ssp(1,0) = 0.d0
    do 100,i=0,nrad-1
    rl=rrad(i)
    rlp1=rrad(i+1)
    t2 = rlp1-rl
    t3 = t2**2
    t4 = 1.d0/t3
    t5 = rlp1**2
    t6 = t5**2
    t8 = rl**2
    t9 = t8**2
    t14 = 0.2D0*t4*(t6*rlp1-t9*rl)
    t15 = 1.d0/t2
    t16 = rl*t15
    t17 = 1.D0+t16
    t20 = t6-t9
    t21 = t17*t15*t20
    t23 = t17**2
    t27 = t5*rlp1-t8*rl
    ssp(1,i) = ssp(1,i) + t14-0.5D0*t21+0.33333333333333D0*t23*t27
    t32 = rl*t4*t20
    ssp(2,i) = -t14+0.25D0*t21+0.25D0*t32-0.33333333333333D0*t16*t27*t17
    ssp(1,i+1) = t14-0.5D0*t32+0.33333333333333D0*t8*t4*t27
100     continue
    return
end
!******************************************************************************************
subroutine crthhp(nrad,lang,znuc,rrad,hhp)
    implicit real(8) (a-h,o-z)
    dimension rrad(0:nrad),hhp(2,0:nrad)
!    creates tridiagonal hamiltonian matrix
    z=znuc
    hllp1=.5d0*lang*(lang+1)
    hhp(1,0) = 0.d0
    do 100,i=0,nrad-1
        rl=rrad(i)
        rlp1=rrad(i+1)
        t2 = rlp1-rl
        t3 = t2**2
        t4 = 1.d0/t3
        t6 = rlp1**2
        t7 = t6**2
        t8 = rl**2
        t9 = t8**2
        t13 = 0.25D0*Z*t4*(t7-t9)
        t15 = 2.D0*hllp1*t4
        t16 = 1.d0/t2
        t18 = 1.D0+rl*t16
        t20 = Z*t18*t16
        t26 = t6*rlp1-t8*rl
        t29 = t18**2
        t33 = hllp1*t18*t16
        t37 = t6-t8
        hhp(1,i) = hhp(1,i)                         &
             -t13+0.16666666666667D0*(t15+t4+4.D0*t20)*t26+   &
                0.25D0*(-2.D0*Z*t29-4.D0*t33)*t37+hllp1*t29*t2
        t44 = Z*rl
        t45 = t44*t4
        t54 = hllp1*rl
        t55 = t54*t4
        hhp(2,i) = t13+0.16666666666667D0*(-t15-t4-2.D0*t20-2.D0*t45)    &
            *t26+0.25D0*(2.D0*t44*t16*t18+2.D0*t33+2.D0*t55)*t37-t54*t18
        hhp(1,i+1) = -t13+0.16666666666667D0*(t15+t4+4.D0*t45)*t26+    &
                    0.25D0*(-2.D0*Z*t8*t4-4.D0*t55)*t37+hllp1*t8*t16

100    continue
    return
end
!******************************************************************************************    
subroutine ctridag(n,hhp,ssp,shtr,shti,x)
!    Solves the linear system of equations 
!    [hhp - (shtr + I*shti)*ssp ] x = y
!    hhp and ssp are both tridiagonal symmetric matrices.
!    Their diagonals elements are stored in hhp(1,i) and ssp(1,i) for i=1, ... n, 
!    their off-diagonal elements in hhp(2,i) and ssp(2,i) for i=1, ... n-1.
!    On input the subroutine argument x contains the right hand side y, 
!    on output it contains the solution x. wrk is a work array. 
    implicit real(8) (a-h,o-z)
    parameter(nwrk=100000)
    complex(8) x,wrk,tt,bet,var
    dimension hhp(2,n),ssp(2,n),x(n)
    dimension wrk(nwrk)
    if (nwrk.lt.n) stop 'enlarge nwrk in ctridag'

    BET=dcmplx(hhp(1,1)-shtr*ssp(1,1),-shti*ssp(1,1))
    var=x(1)/BET
    x(1)=var
    DO 11 J=2,N
    tt=dcmplx(hhp(2,J-1)-shtr*ssp(2,J-1),-shti*ssp(2,J-1))
    wrk(J)=tt/BET
    BET=dcmplx(hhp(1,J)-shtr*ssp(1,J),-shti*ssp(1,J)) - tt*wrk(J)
    var=(x(J)-tt*x(J-1))/BET
    x(J)=var
11    CONTINUE
    DO 12 J=N-1,1,-1
    var=x(J)-wrk(J+1)*x(J+1)
    x(J)=var
12    CONTINUE
    RETURN
END
!******************************************************************************************
SUBROUTINE LSD_PADE(RHO,ETA,EC,VCA,VCB)
!     Evaluates the exchange correlation energy and potential in the Pade Approximation at a single point
!     Input: rho: total charge density
!        eta: spin polarization, eta = (rho_up - rho_down)/(rho_up + rho_down)
!     Output: ec: exchange correlation energy density
!        vca: exhange correlation potential for spin up orbitals
!        vcb: exhange correlation potential for spin down orbitals
    IMPLICIT REAL(8) (A-H,O-Z)
    PARAMETER (A0=.4581652932831429d0,A1=2.217058676663745d0,A2=0.7405551735357053d0,A3=0.01968227878617998d0)
    PARAMETER (B1=1.0D0,B2=4.504130959426697d0,B3=1.110667363742916d0,B4=0.02359291751427506d0)
    PARAMETER (DA0=.119086804055547D0,DA1=.6157402568883345d0,DA2=.1574201515892867d0,DA3=.003532336663397157d0)
    PARAMETER (DB1=0.0d0,DB2=.2673612973836267d0,DB3=.2052004607777787d0,DB4=.004200005045691381d0)
    PARAMETER (RSFAC=.6203504908994000d0)
    PARAMETER (FSFAC=1.92366105093153617d0)

    RS=RSFAC*RHO**(-1.d0/3.d0)
    FS=FSFAC*((1.d0+ETA)**(4.d0/3.d0)+(1.d0-ETA)**(4.d0/3.d0)-2.d0)
    DFS=FSFAC*4.d0/3.d0*((1.d0+ETA)**(1.d0/3.d0)-(1.d0-ETA)**(1.d0/3.d0))
    DFSA=DFS*(1.d0-ETA)
    DFSB=DFS*(-1.d0-ETA)
    A0P=A0+FS*DA0
    A1P=A1+FS*DA1
    A2P=A2+FS*DA2
    A3P=A3+FS*DA3
    B1P=B1+FS*DB1
    B2P=B2+FS*DB2
    B3P=B3+FS*DB3
    B4P=B4+FS*DB4
    TOP=A0P+RS*(A1P+RS*(A2P+RS*A3P))
    DTOP=A1P+RS*(2.d0*A2P+RS*3.d0*A3P)
    TOPX=DA0+RS*(DA1+RS*(DA2+RS*DA3))
    BOT=RS*(B1P+RS*(B2P+RS*(B3P+RS*B4P)))
    DBOT=B1P+RS*(2.d0*B2P+RS*(3.d0*B3P+RS*4.d0*B4P))
    BOTX=RS*(DB1+RS*(DB2+RS*(DB3+RS*DB4)))
    EC=-TOP/BOT
    VC=EC+RS*(DTOP/BOT-TOP*DBOT/(BOT*BOT))/3.d0
    DX=-(TOPX/BOT-TOP*BOTX/(BOT*BOT))
    VCA=VC+DX*DFSA
    VCB=VC+DX*DFSB
    RETURN
END
!    ******************************************************************************************
subroutine tmatmul(nrad,hhp,xrad,yrad)
    implicit real(8) (a-h,o-z)
    dimension xrad(0:nrad),yrad(0:nrad),hhp(2,0:nrad)
!    multiplies a vector with a tridiagonal matrix
    i=0
    yrad(i) = hhp(1,i)*xrad(i) + hhp(2,i)*xrad(i+1)
    do 100, i=1,nrad-1
        yrad(i) = hhp(2,i-1)*xrad(i-1) + hhp(1,i)*xrad(i) + hhp(2,i)*xrad(i+1)
100     continue
    i=nrad
    yrad(i) = hhp(2,i-1)*xrad(i-1) + hhp(1,i)*xrad(i) 
    return
end
!******************************************************************************************
