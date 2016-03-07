PROGRAM LJ_CG
! purpose:
! This program minimize the Lennard Jones potential of cluster of particles
! through the conjugate gradient method. since the Lennard Jones potential is 
! nonlinear, we used line search to find the scale length (alph).
! for calculating the energy and force of the system this program nees to call 
! force_energy subroutine which is beside the main program file (LJ.f90)
! to run simply type:
!
!   ifort CG.f90 LJ.f90
!
! Result: it, energy, max(force), difference of energy 
!   31  -1.739284226044216E+02   4.090742742692915E-03  -2.424023506364392E-06
IMPLICIT NONE
REAL(8), allocatable :: rat(:,:)    ! coordinates of atoms
REAL(8), allocatable :: fat(:,:)    ! force between atoms
INTEGER :: nat                      ! number of atoms
REAL(8) :: epot         ! potential energy
REAL(8) :: ftot                     ! total force on atoms
REAL(8) :: arr(3)
CHARACTER(len=20) :: filename  ! Input data file name
CHARACTER(len=3)  :: sat            ! for reading the file
INTEGER :: status                   ! I/O status: 0 for success
INTEGER :: n
INTEGER :: iat                      ! local variables
integer :: clck_counts_beg, clck_counts_end, clck_rate


filename ="s2.xyz"

! Open input data file.  
OPEN ( UNIT=21, FILE=filename, status='OLD', ACTION='READ', &
       IOSTAT=status )
READ (21,*) nat
READ (21,*)

allocate (rat(3,nat))
allocate (fat(3,nat))

fileopen: IF ( status == 0 ) THEN       ! Open successful
    DO iat = 1, nat
        READ (21,*) sat, rat(1, iat), rat(2, iat), rat(3, iat)
    END DO
    CLOSE(21)   
END if fileopen
call system_clock ( clck_counts_beg, clck_rate )
call force_energy(nat,epot,rat,fat)
call system_clock ( clck_counts_end, clck_rate )
write(*,*) epot
write (*, *) 'elapsed wall clock time in seconds'
write (*, *)  (clck_counts_end - clck_counts_beg) / real (clck_rate)

CONTAINS

!*****************************************************************
!*****************************************************************
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

END PROGRAM LJ_CG