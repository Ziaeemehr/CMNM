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


filename ="posinp_1000.xyz"

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

! Calling Conjugate Gradient Algorithm
CALL CG(nat,rat,epot,fat)

CONTAINS

!*****************************************************************
!*****************************************************************

SUBROUTINE CG(nat,rat,epot,fat)

IMPLICIT NONE
INTEGER :: nat, it
REAL(8) :: rat(3,nat), rat1(3,nat), h(3,nat)
REAL(8) :: fat(3,nat), fat1(3,nat), fat2(3,nat), df(3,nat)
REAL(8) :: epot, epot1,epot2, s, b, a
REAL(8) :: de, alphmin, alph_t, fmax, bet

! number of intervals

      print *
      print *,'  Conjugate Gradient method'
      print *
      write(*,'(3x,a4,5x,a10,12x,a10,18x,a8,17x,a10)')'iter','epot2','de','fmax','alphamin'

! Conjugate Gradient Algorithm
CALL force_energy(nat,epot,rat,fat)
h = fat
alph_t=2.d-3
! Main loop
DO it=1,1000
  rat1=rat+alph_t*fat
  CALL force_energy(nat,epot1,rat1,fat1)
  b  = sum(fat*fat)
  s  = sum(fat1*fat)
  a  = (s-b)/(2.d0*alph_t)
  alphmin = -b/(2.d0*a)
  rat=rat+alphmin*h
  CALL force_energy(nat,epot2,rat,fat2)
  df=fat2-fat
  bet=sum(df*fat2)/sum(fat*fat)
  h=fat2+bet*h
  de=epot2-epot
  fmax=maxval(abs(fat2))
  fat=fat2            !!! I forgot this
  epot=epot2
  write(*,'(i6,es24.15,2es24.15,2es24.15)') it,epot2,de,alphmin,fmax
  IF (fmax < 5.d-3) exit
ENDDO
END SUBROUTINE CG

!*****************************************************************
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

END PROGRAM LJ_CG