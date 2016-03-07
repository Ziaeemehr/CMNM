PROGRAM optim_SD
!
! purpose:
!       This program minimize the potential energy of a system of
!       atoms which interact with each other via lennard-
!       jones potential U=a0(1/r^12-1/1^6)
!
IMPLICIT NONE
REAL(8), allocatable :: rat(:,:) ! coordinates of atoms
REAL(8), allocatable :: fat(:,:) ! force between atoms
INTEGER :: nat                   ! number of atoms
REAL(8) :: epot                  ! potential energy
REAL(8) :: ftot                  ! total force on atoms
REAL(8) :: arr(3)                ! 3x1 random numbers
CHARACTER(len=20) :: filename    ! Input data file name
CHARACTER(len=3)  :: sat         ! for reading the file of coordinates
INTEGER :: status                ! I/O status: 0 for success
INTEGER :: iat, jat              ! loop variables
!INTEGER :: filenumber

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

! Displace atoms from initial coordinates which are already in minimum
! position of potential energy
CALL randomize_r(nat,rat)

! Calling Steepest Descent algorithm
CALL SD(nat,rat,epot,fat)

! Displace atoms from minimum positions
CALL randomize_r(nat,rat)

! Calling Steepest Descent with Energy Feedback algorithm
CALL EF(nat,rat,epot,fat)

! Displace atoms from minimum positions
CALL randomize_r(nat,rat)

! Calling Steepest Descent with Gradient Feedback algorithm
CALL GR(nat,rat,epot,fat)

CONTAINS

!*****************************************************************
SUBROUTINE randomize_r(nat,rat)
!
! The purpose of this subroutine is to displace atoms from initial minimum 
! positions, interval of displacement is between [-0.05,0.05]
!
IMPLICIT NONE
REAL(8) :: rat(3,nat)
REAL(8) :: arr(3)
INTEGER :: iat, nat

DO iat =1,nat
    CALL random_number(arr)
    arr(1:3) = (arr(1:3)-0.5)*0.1
    rat(1:3,iat) = rat(1:3,iat)+ arr(1:3)
ENDDO
END SUBROUTINE randomize_r

!*****************************************************************

SUBROUTINE SD(nat,rat,epot,fat)
!
! The purpose of this subrutine is to minimize the Lennard-Jones potential energy
! of atoms by use of  Steepest Descent Method
!
IMPLICIT NONE
REAL(8) :: alph
INTEGER :: nat, it
REAL(8) :: fat(3,nat), rat(3,nat) 
REAL(8) :: epot
REAL(8) :: epot0                 ! potential energy
REAL(8) :: fmax
REAL(8) :: de    
      print *
      print *,' Steepest Descent method'
      print *

alph = 1.d-3
DO it = 0,1000
    CALL force_energy(nat,epot,rat,fat)

    if(it==0) epot0 = epot
fmax = maxval(abs(fat))
de   = epot-epot0 
write(*,'(i6,es24.15,2es24.15)') it,epot,de,fmax

IF (fmax < 5.d-3) exit
rat(1:3,1:nat) = rat(1:3,1:nat)+alph*fat (1:3,1:nat)
epot0 = epot
ENDDO

END SUBROUTINE SD 
!*****************************************************************

SUBROUTINE EF(nat,rat,epot,fat)
!
! The purpose of this subrutine is to minimize the Lennard-Jones potential energy
! of atoms by use of  Steepest Descent with energy feedback method
!
IMPLICIT NONE
REAL(8) :: alph
INTEGER :: nat, it                   ! number of atoms
REAL(8) :: fat(3,nat), rat(3,nat) 
REAL(8) :: epot
REAL(8) :: epot0           ! potential energy
REAL(8) :: fmax
REAL(8) :: de
      print *      
      print *,' Steepest Descent method with Energy Feedback'
      print *

alph = 1.d-3
DO it = 0,1000
    CALL force_energy(nat,epot,rat,fat)

    if(it==0) epot0 = epot
    fmax = maxval(abs(fat))
    de   = epot-epot0 
    write(*,'(i6,es24.15,2es24.15)') it,epot,de,fmax

    IF (fmax < 5.d-3) exit
    IF (de < 0) alph = alph * (1+5.d-2) 
    IF (de > 0) alph = alph * 0.5
    rat(1:3,1:nat) = rat(1:3,1:nat)+alph* fat (1:3,1:nat)
    epot0 = epot
ENDDO
END SUBROUTINE EF 

!*****************************************************************

SUBROUTINE GR(nat,rat,epot,fat)
!
! The purpose of this subrutine is to minimize the Lennard-Jones potential energy
! of atoms by use of  Steepest Descent with gradient feedback method
!
IMPLICIT NONE
INTEGER :: nat, it
REAL(8) :: alph
REAL(8) :: fat(3,nat), rat(3,nat), fat0(3,nat)
REAL(8) :: epot
REAL(8) :: fmax
REAL(8) :: t

      print *      
      print *,' Steepest Descent method with Gradient Feedback'
      print *
      print *,'    I      Potential Energy         COS(theta)               Maximum Force '
      print *
alph = 1.d-3
DO it = 0,1000
    CALL force_energy(nat,epot,rat,fat)

    if(it==0) fat0 = fat
    fmax = maxval(abs(fat))
    write(*,'(i6,es24.15,2es24.15)') it,epot,t,fmax
    IF (fmax < 5.d-3) exit 
    t  = cos_thet(fat,fat0,nat)
    IF (t < 0.5) alph =  alph * 0.5
    IF (t > 0.5) alph =  alph * (1.05) 
    rat(1:3,1:nat) = rat(1:3,1:nat) + alph*fat(1:3,1:nat)
    fat0 = fat
ENDDO
END SUBROUTINE GR

!*****************************************************************
Real FUNCTION cos_thet (x,y,nat)
!
! This function evaluate the cos(theta) between two vector x and y
!
IMPLICIT NONE
INTEGER :: nat
REAL(8), INTENT(IN) :: x(3,nat), y(3,nat)
REAL(8) :: denom1, denom2, thet

thet   = sum(x(1:3,1:nat)*y(1:3,1:nat))
denom1 = sum(x(1:3,1:nat)*x(1:3,1:nat))
denom2 = sum(y(1:3,1:nat)*y(1:3,1:nat))

cos_thet = thet/(sqrt(denom1)*sqrt(denom2))

END FUNCTION cos_thet

!*****************************************************************

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
END PROGRAM optim_SD
