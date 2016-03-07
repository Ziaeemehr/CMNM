program main
! purpose:
!       This program minimize the potential energy of a system of
!       atoms which interact with each other via lennard-
!       jones potential U=a0(1/r^12-1/1^6)
!
IMPLICIT NONE
real(8), allocatable :: rat(:,:) ! coordinates of atoms
real(8), allocatable :: fat(:,:) ! force between atoms
INTEGER :: nat                   ! number of atoms
real(8) :: epot                  ! potential energy
real(8) :: ftot                  ! total force on atoms
real(8) :: arr(3)                ! 3x1 random numbers
CHARACTER(len=20) :: filename    ! Input data file name
CHARACTER(len=3)  :: sat         ! for reading the file of coordinates
INTEGER :: status                ! I/O status: 0 for success
INTEGER :: iat, jat              ! loop variables

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

CALL randr(nat,rat)
CALL cg(nat,rat,epot,fat)

CONTAINS

subroutine cg(nat,rat,epot,fat)
implicit none
real(8),dimension(3,nat) :: rat, rat0
real(8),dimension(3,nat) :: h, g, gold
real(8):: ggold
real(8), dimension(3,nat) :: fat 
integer, INTENT(IN) :: nat
integer :: it, i
real(8) :: epot, epot0
call lj(nat, rat, epot, fat)

h = fat
rat0 = rat
do it = 1, 1000
    if ( it == 1 ) epot = epot0
    
    call line( nat, rat, epot, fat, alph )
    rat =  rat0 + alph * h
    rat0 = rat
    gold = g
    ggold =  sum(g(1:3, 1,nat) * g(1:3, 1,nat))
    call lj(nat, rat, epot, fat)
    fmax = maxval(abs(fat))
    de = epot - epot0
    write(*,'(1x,I5,3es24.15)') it , epot, fmax, de
    if ( fmax < 5.d-3 ) exit
    g   =  -fat
    dg  =  g - gold
    bet = 0.d0
    do i = 1, nat
        bet = bet + dg( 1:3,i ) * g( 1:3, i )
    enddo
    bet = bet / ggold
    epot0 = epot
enddo
end subroutine

subroutine line(nat,rat,epot,fat)

IMPLICIT NONE
INTEGER :: nat, it
REAL(8) :: alphmin,alpht
REAL(8) :: fat(3,nat), rat(3,nat), fat_t(3,nat),rat_t(3,nat)
REAL(8) :: epot, epot0, a, b, de

alpht = 2.d-3
CALL lj(nat,epot0,rat,fat)
rat_t(1:3,1:nat) = rat(1:3,1:nat)+alpht*fat(1:3,1:nat)
CALL lj(nat,epot,rat_t,fat_t)
de = epot-epot0
b  = sum(fat(1:3,1:nat)*fat(1:3,1:nat))
a  = sum(fat_t(1:3,1:nat)*fat(1:3,1:nat))
a  = (a-b)/(2*alpht)
alphmin = -b/(2*a)    
end subroutine
SUBROUTINE lj(nat,epot,rat,fat)
!
! The purpose of this subrutine is to evaluate the Lennard-Jones potential energy
! and force between atoms.
!
IMPLICIT NONE
real(8), INTENT(IN) :: rat(3,nat) 
real(8), INTENT(OUT):: fat(3,nat) 
INTEGER, INTENT(IN) :: nat
CHARACTER(len=3)    :: sat
real(8),INTENT(out) :: epot
INTEGER :: iat, jat, k  
real(8) :: r2, r6, r8, r12, r14, dx, dy, dz, d
real(8) :: ftot

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
END SUBROUTINE lj

SUBROUTINE randr(nat,rat)

IMPLICIT NONE
real(8) :: rat(3,nat)
real(8) :: arr(3)
INTEGER :: iat, nat

DO iat =1,nat
    CALL random_number(arr)
    arr(1:3) = (arr(1:3)-0.5)*0.1
    rat(1:3,iat) = rat(1:3,iat)+ arr(1:3)
ENDDO
END SUBROUTINE randr

end program



