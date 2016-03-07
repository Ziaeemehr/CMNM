program main
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
real(8),dimension(3,nat) :: rat, rat1, rat0
real(8),dimension(3,nat) :: h
real(8),dimension(3,nat) :: fat, fat1, fat2, df
integer :: nat,it
real(8) :: epot, epot1, epot2, fmax, de, bet, alphmin, alph_t, s, a, b
call lj(nat, epot, rat, fat)
h = fat
alph_t = 2.d-3
do it = 1, 1000
    rat1 = rat +  alph_t * fat
    call lj(nat, epot1, rat1, fat1)
    b = sum(fat*fat)
    s =  sum(fat1*fat)
    a = ( s-b)/(2.d0*alph_t)
    alphmin = -b/(2.d0*a)
    rat =  rat + alphmin * h
    call lj(nat, epot2,rat,fat2)
    df = fat2 - fat
    bet = sum (df * fat2)/sum(fat*fat)
    h = fat2 + bet * h
    de = epot2 - epot
    fmax = maxval(abs(fat2))
    fat = fat2
    epot =  epot2
    write(*,'(I6,4es24.15)') it , epot2,de, alphmin, fmax
    if ( fmax < 5.d-5 ) exit
enddo
end subroutine

SUBROUTINE lj(nat,epot,rat,fat)
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



