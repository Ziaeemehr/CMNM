!Molecular Dynamics Simulation of Lennard-Jones gas

!---------------------------------------------------------------
!define constants
MODULE constants
IMPLICIT NONE
REAL(8) :: rcut1 !first cut-off radius
PARAMETER (rcut1 = 5.5)
REAL(8) :: rcut2 !second cut-off radius
PARAMETER (rcut2 = 6.d0)
REAL(8) :: V !cut-off energy
PARAMETER (V=4.0*(1.0/rcut1**12-1.0/rcut1**6)) !(V=4.0/(rcut1**6)*(1/(rcut1**6)-1))
REAL(8) :: Vd !first derivative of potential
PARAMETER (Vd=-6.0*(8.0/rcut1**13-4.0/rcut1**7))!(Vd=-6.d0*(V+4.d0/(rcut1**12))/rcut1)
REAL(8) :: Vdd !second derivative of potential
PARAMETER (Vdd=24.d0*(26.d0/rcut1**14-7.d0/rcut1**8))!(Vdd=24.d0/rcut1*(26.d0/rcut1**12-7.d0/rcut1**6))
REAL(8) :: a0
PARAMETER (a0=(rcut2**3*(-2.d0*(10.d0*rcut1**2 - 5.d0*rcut1*rcut2 + rcut2**2)*V + 2.d0*rcut1*(4.d0*rcut1**2 - 5.d0*rcut1*rcut2 + rcut2**2)*Vd - rcut1**2*(rcut1 - rcut2)**2*Vdd))/(2.d0*(rcut1 - rcut2)**5))
REAL(8) :: a1
PARAMETER (a1=(rcut2**2*(60.d0*rcut1**2*V - 2.d0*(rcut1 - rcut2)*(6.d0*rcut1 - rcut2)*(2.d0*rcut1 + rcut2)*Vd + rcut1*(rcut1 - rcut2)**2*(3.d0*rcut1 + 2.d0*rcut2)*Vdd))/(2.d0*(rcut1 - rcut2)**5))
REAL(8) :: a2
PARAMETER (a2=(12.d0*rcut1*rcut2*(-5.d0*(rcut1 + rcut2)*V + (rcut1 - rcut2)*(2.d0*rcut1 + 3.d0*rcut2)*Vd) - (rcut1 - rcut2)**2*rcut2*(3.d0*rcut1**2 + 6.d0*rcut1*rcut2 + rcut2**2)*Vdd)/(2.d0*(rcut1 - rcut2)**5))
REAL(8) :: a3
PARAMETER (a3=(20.d0*(rcut1**2 + 4.d0*rcut1*rcut2 + rcut2**2)*V - 4.d0*(rcut1 - rcut2)*(2.d0*rcut1**2 + 10.d0*rcut1*rcut2 + 3.d0*rcut2**2)*Vd + (rcut1 - rcut2)**2*(rcut1**2 + 6.d0*rcut1*rcut2 + 3.d0*rcut2**2)*Vdd)/(2.d0*(rcut1 - rcut2)**5))
REAL(8) :: a4
PARAMETER (a4=-(30.d0*(rcut1 + rcut2)*V - 2.d0*(rcut1 - rcut2)*(7.d0*rcut1 + 8.d0*rcut2)*Vd + (rcut1 - rcut2)**2*(2.d0*rcut1 + 3.d0*rcut2)*Vdd)/(2.d0*(rcut1 - rcut2)**5))
REAL(8) :: a5
PARAMETER (a5=(12.d0*V + (rcut1 - rcut2)*(-6.d0*Vd + (rcut1 - rcut2)*Vdd))/(2.d0*(rcut1 - rcut2)**5))
END MODULE constants
