program energy_ll
! cell linked list for Lennard Jones potential
USE constants
IMPLICIT NONE
REAL(8), allocatable :: rat(:,:)      ! coordinates of atoms
INTEGER :: nat,it                     ! number of atoms
REAL(8) :: epot                       ! potential energy
CHARACTER(len=20) :: filename         ! Input data file name
CHARACTER(len=3)  :: sat              ! for reading the file
INTEGER :: status                     ! I/O status: 0 for success
INTEGER :: i,j,k,iat                  ! local variables
INTEGER :: nx,ny,nz                   !number of cells per dimension
REAL(8) :: hx,hy,hz                   !size of a cell
REAL(8) :: Lx,Ly,Lz                   !length of the box
REAL(8) :: xmin,xmax,ymin,ymax,zmin,zmax
INTEGER :: cellx,celly,cellz          !id of the cell for certain atom
INTEGER,allocatable :: link(:)
INTEGER, allocatable :: cell(:,:,:)
INTEGER :: clck_counts_beg, clck_counts_end, clck_rate

filename ="posinp_1000.xyz"

! Open input data file.  
OPEN ( UNIT=21, FILE=filename, status='OLD', ACTION='READ', &
       IOSTAT=status )
READ (21,*) nat
READ (21,*)

allocate (rat(3,nat))
allocate (link(nat))

fileopen: IF ( status == 0 ) THEN    ! Open successful
    DO iat = 1, nat
        READ (21,*) sat, rat(1, iat), rat(2, iat), rat(3, iat)
    END DO
    CLOSE(21)   
END if fileopen
call system_clock ( clck_counts_beg, clck_rate )
!set up the cell list
xmin=minval(rat(1,1:nat)); ymin=minval(rat(2,1:nat));zmin=minval(rat(3,1:nat))
xmax=maxval(rat(1,1:nat)); ymax=maxval(rat(2,1:nat));zmax=maxval(rat(3,1:nat))
Lx=6+xmax-xmin
Ly=6+ymax-ymin
Lz=6+zmax-zmin
nx=int(Lx/rcut2)
ny=int(Ly/rcut2)
nz=int(Lz/rcut2)
!write(*,*) nx,ny,nz
hx=Lx/nx
hy=Ly/ny
hz=Lz/nz
do j=1,nat
     rat(1,j)=rat(1,j)+3-xmin
     rat(2,j)=rat(2,j)+3-ymin
     rat(3,j)=rat(3,j)+3-zmin
end do
allocate (cell(1:nx,1:ny,1:nz))
call cell_list (rat,nat,cell,nx,ny,nz,hx,hy,hz,link)
call force_energy (rat,nat,cell,link,nx,ny,nz,hx,hy,hz,epot)
call system_clock ( clck_counts_end, clck_rate )
write(*,*) epot
write(*,*) 'Elapsed wall clock time in seconds'
write (*,'(es12.5)' )  (clck_counts_end - clck_counts_beg) / real (clck_rate)

contains
! =====================================================================================
SUBROUTINE cell_list (rat,nat,cell,nx,ny,nz,hx,hy,hz,link)

USE constants
IMPLICIT NONE
INTEGER :: nat
REAL(8) :: rat(3,nat)
INTEGER,intent(inout) :: nx,ny,nz !number of cells per dimension
REAL(8),intent(inout) :: hx,hy,hz !size of a cell
REAL(8) :: Lx,Ly,Lz !length of the box
REAL(8) :: xmin,xmax,ymin,ymax,zmin,zmax
INTEGER :: iat !label of particles
INTEGER :: j,k !loop control
INTEGER :: cellx,celly,cellz !id of the cell for certain atom
INTEGER :: link(nat)
INTEGER :: cell(1:nx,1:ny,1:nz)

cell = 0
link = 0
DO iat = 1,nat
    cellx = int(rat(1,iat)/hx)+1
    celly = int(rat(2,iat)/hy)+1
    cellz = int(rat(3,iat)/hz)+1    
    if (cell(cellx,celly,cellz).ne.0) then
         link(iat) = cell(cellx,celly,cellz)
    endif
    cell(cellx,celly,cellz)=iat
END DO
END SUBROUTINE cell_list
! =====================================================================================
SUBROUTINE force_energy (rat,nat,cell,link,nx,ny,nz,hx,hy,hz,epot)
USE constants
IMPLICIT NONE
INTEGER :: nat
REAL(8) :: rat(3,nat)                       !position matrix
INTEGER :: link(nat)                        !link list
INTEGER :: nx,ny,nz                         !number of cells per dimension
REAL(8) :: hx,hy,hz                         !size of cells
REAL(8) :: force(3,nat)                     !force matrix
INTEGER :: iat                              !labels for atoms
INTEGER :: cellx,celly,cellz                !id of the cell for certain atom
REAL(8) :: rat1(3)                          !the positon of the other atom that interacts with the current atom
INTEGER :: i,j,k,m,id                       !loop control
REAL(8) :: distance(3)                      !distance in each dimensions
REAL(8) :: radius,radius2,radius6,radius6i  !square of the distance
REAL(8) :: epot,uij                         !potential energy for each atom
INTEGER :: cell(1:nx,1:ny,1:nz)             !cell list
force=0
epot=0
DO iat = 1,nat                                 
    cellx = int(rat(1,iat)/hx)+1            !find which cell this atom sits in
    celly = int(rat(2,iat)/hy)+1
    cellz = int(rat(3,iat)/hz)+1 
                                            !scan the nearest 27 cells
    DO i = max(cellx-1,1),min(cellx+1,nx)
        DO j = max(celly-1,1),min(celly+1,ny)
            DO k = max(cellz-1,1),min(cellz+1,nz)
                id=cell(i,j,k)
                DO WHILE (id .NE. 0)        !if not reached the tail of the chain
                    if (id==iat) then
                            id=link(id)
                            cycle
                    endif
                    rat1(1)=rat(1,id)
                    rat1(2)=rat(2,id)
                    rat1(3)=rat(3,id)
                    distance(1)=rat(1,iat)-rat1(1)
                    distance(2)=rat(2,iat)-rat1(2)
                    distance(3)=rat(3,iat)-rat1(3)
                    !if the atom is outside the cut-off radius, continue to the next atom
                    radius2 = distance(1)*distance(1)+distance(2)*distance(2)+distance(3)*distance(3)
                    radius  = sqrt(radius2)
                    IF (radius2 > rcut2*rcut2) THEN
                        id = link(id)                        
                        CYCLE
                    ELSEIF ((radius2 >= rcut1*rcut1) .and. (radius2 <= rcut2*rcut2)) THEN
                                            !Fifth order swithing polynomial in horner's form 
                        uij= a0+radius*(a1+radius*(a2+radius*(a3+radius*(a4+radius*a5))))
                        epot = epot + uij
                        id = link(id)
                    ELSE                    !(radius < rcut1)
                        radius6 = radius2**3
                        radius6i = 1.0/radius6
                                            !calculate potential energy
                        uij =4.d0*(radius6i**2-radius6i) !4.0*radius6i*(radius6i-1)
                        epot = epot + uij
                        id = link(id)       !move on to the next atom in the chain                        
                    ENDIF                    
                END DO                
            END DO
        END DO
    END DO
END DO
epot=epot/2
END SUBROUTINE force_energy
! =====================================================================================
end program energy_ll
