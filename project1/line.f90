PROGRAM SDLINE
IMPLICIT NONE
REAL(8), allocatable :: rat(:,:) ! coordinates of atoms
REAL(8), allocatable :: fat(:,:) ! force between atoms
INTEGER :: nat                   ! number of atoms
REAL(8) :: epot                  ! potential energy
REAL(8) :: ftot                  ! total force on atoms
REAL(8) :: arr(3)
CHARACTER(len=20) :: filename    ! Input data file name
CHARACTER(len=3)  :: sat         ! for reading the file
INTEGER :: status                ! I/O status: 0 for success
INTEGER :: n
INTEGER :: iat, jat              ! local variables

! Get the name of the file containing the input data.
WRITE (*,1000) 
1000 FORMAT (1X,'Enter the file name with the data to be sorted: ')
WRITE (*,*) 'The file name is: 10,38,100,150,629,1000'

READ (*,'(A)') filename

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

DO iat =1,nat
    CALL random_number(arr)
    arr(1:3) = (arr(1:3)-0.5)*0.1
    rat(1:3,iat) = rat(1:3,iat)+ arr(1:3)
ENDDO
CALL LINE(nat,rat,epot,fat)

CONTAINS
!*****************************************************************

SUBROUTINE LINE(nat,rat,epot,fat)
!
! The purpose of this subrutine is to minimize the Lennard-Jones potential energy
! of atoms by use of  Steepest Descent with Line minimization
!

IMPLICIT NONE
INTEGER :: nat, it
REAL(8) :: alphmin,alpht
REAL(8) :: fat(3,nat), rat(3,nat), fat_t(3,nat),rat_t(3,nat)
REAL(8) :: epot, epot0, a, b, de
REAL(8) :: fmax
REAL(8) :: t

      print *      
      print *,' Line search Method '
      print *
alpht = 2.d-3

DO it = 0,1000
    CALL force_energy(nat,epot0,rat,fat)
    fmax = maxval(abs(fat))
    IF(it==0) epot = epot0
    write(*,'(i6,es24.15,2es24.15,2es24.15)') it,epot,alphmin/alpht,fmax,de
    IF (fmax < 5.d-3) exit
    rat_t(1:3,1:nat) = rat(1:3,1:nat)+alpht*fat(1:3,1:nat)
    CALL force_energy(nat,epot,rat_t,fat_t)
    de=epot-epot0
    b  = sum(fat(1:3,1:nat)*fat(1:3,1:nat))
    a  = sum(fat_t(1:3,1:nat)*fat(1:3,1:nat))
    a = (a-b)/(2*alpht)
    alphmin = -b/(2*a)    
    rat(1:3,1:nat) = rat(1:3,1:nat)+alphmin*fat(1:3,1:nat)
ENDDO

END SUBROUTINE LINE

!*****************************************************************

SUBROUTINE force_energy(nat,epot,rat,fat)

IMPLICIT NONE
REAL(8), INTENT(IN) :: rat(3,nat) 
REAL(8), INTENT(OUT) :: fat(3,nat) 
INTEGER, INTENT(IN) :: nat
CHARACTER(len=3)    :: sat
INTEGER     :: iat, jat, k     ! local variables
REAL(8),INTENT(out)     :: epot
REAL(8) :: r2, r6, r8, r12, r14, dx, dy, dz, d
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
END PROGRAM SDLINE
