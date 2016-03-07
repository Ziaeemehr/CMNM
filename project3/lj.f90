SUBROUTINE force_energy(nat,epot,rat,fat)

IMPLICIT NONE
REAL(8), INTENT(IN)  :: rat(3,nat) 
REAL(8), INTENT(OUT) :: fat(3,nat) 
INTEGER, INTENT(IN)  :: nat
CHARACTER(len=3)     :: sat
REAL(8),INTENT(out)  :: epot
INTEGER :: iat, jat, k     ! local variables
REAL(8) :: r , dx, dy, dz, d
REAL(8) :: ftot
INTEGER :: status              ! I/O status: 0 for success

fat=0
epot = 0
ftot = 0

DO iat = 1, nat-1
    DO jat = iat+1, nat
            dx = rat(1,jat)-rat(1,iat)
            dy = rat(2,jat)-rat(2,iat)
            dz = rat(3,jat)-rat(3,iat)
            r  = dx**2 + dy**2 + dz**2      ! its r^2
            d  = 4.d0*(-12.d0/r**7 + 6.d0/r**4)       
            epot = epot+4.d0*(1.d0/r**6-1.d0/r**3)
            fat(1,iat) = fat(1,iat) + d * dx  
            fat(2,iat) = fat(2,iat) + d * dy  
            fat(3,iat) = fat(3,iat) + d * dz  
            fat(1,jat) = fat(1,jat) - d * dx  
            fat(2,jat) = fat(2,jat) - d * dy  
            fat(3,jat) = fat(3,jat) - d * dz  
    ENDDO
ENDDO
END SUBROUTINE force_energy
