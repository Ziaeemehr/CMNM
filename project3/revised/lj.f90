SUBROUTINE force_energy(nat,epot,rat,fat)
!
! The purpose of this subrutine is to evaluate the Lennard-Jones potential energy
! and force between atoms.
!
IMPLICIT NONE
REAL(8), INTENT(IN) :: rat(3,nat) 
REAL(8), INTENT(OUT):: fat(3,nat) 
INTEGER, INTENT(IN) :: nat
CHARACTER(len=3)    :: sat
REAL(8),INTENT(out) :: epot
INTEGER :: iat, jat, k  
REAL(8) :: r2, r6, r8, r12, r14, dx, dy, dz, d
REAL(8) :: ftot

fat=0
epot = 0
ftot = 0

DO iat = 1, nat-1
    DO jat = iat+1, nat
            dx = rat(1,jat)-rat(1,iat)
            dy = rat(2,jat)-rat(2,iat)
            dz = rat(3,jat)-rat(3,iat)
            r2 = dx*dx + dy*dy + dz*dz
        r6  = r2 * r2 * r2
        r8  = r2 * r6
        r14 = r6 * r8
        r12 = r6 * r6
            d   = 4*(-12/r14 + 6/r8)       
            epot= epot+4*(1/r12-1/r6)
            fat(1,iat) = fat(1,iat) + d * dx  
            fat(2,iat) = fat(2,iat) + d * dy  
            fat(3,iat) = fat(3,iat) + d * dz  
            fat(1,jat) = fat(1,jat) - d * dx  
            fat(2,jat) = fat(2,jat) - d * dy  
            fat(3,jat) = fat(3,jat) - d * dz  
    ENDDO
ENDDO
END SUBROUTINE force_energy