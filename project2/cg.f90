PROGRAM CG_poisson
!
! purpose:
! Solving 1D poisson equation from Conjugate Gradient method
!
IMPLICIT NONE
REAL(8), allocatable :: v(:),v_ex(:)          ! Potentials at grid points
REAL(8), allocatable :: g(:)          ! 
REAL(8), allocatable :: rho(:), h(:)
INTEGER :: n, ip, i, it               ! loop variables
REAL(8) :: pi , alph, x, norm_g, dt,gg, ggold, bet, denom

! number of intervals
n=100
allocate (v(n-1),g(n-1),rho(n-1),h(n-1),v_ex(n-1))

      print *
      print *,'  Conjugate Gradient method'
      print *
! grid length scale
dt    = 1.d0/(n)
pi = 4*atan(1.d0)

DO ip = 1,n-1
  x = ip*dt
  rho(ip) = 10.d0* sin(10.d0*sin(pi*x))*(cos(pi*x))**2+sin(pi*x)*cos(10.d0*sin(pi*x))
  V_ex(ip) = 1.d-1*(sin(10.d0*sin(pi*x)))/pi**2+x
ENDDO
! Conjugate Gradient Algorithm
g = -dt**2*rho 
g(n-1)=-1 -dt**2*rho(n-1)

h = -g
! Boundary Conditions
v(1:n-1)= 0.d0
! Main loop
DO it=1,1000
  gg   = DOT_PRODUCT(g,g)
  denom = h(1)*(2*h(1)-h(2))              ! for i = 1 (the first and last term are calculated manually)
  DO i = 2,n-2
  denom = denom+(h(i)*(-h(i-1)+2*h(i)-h(i+1)))
  ENDDO 
  denom=denom + h(n-1)*(2*h(n-1)-h(n-2))  !for i = n-1
  alph = gg/denom
  v = v + alph*h
  ggold = gg
  g(1)=g(1)+ alph*(2*h(1)-h(2))           !for i = 1
  DO i = 2,n-2
    g(i)= g(i) + alph*(-h(i-1)+2*h(i)-h(i+1))
  ENDDO
  g(n-1)=g(n-1)+alph*(2*h(n-1)-h(n-2))    !for i = n-1
  gg=DOT_PRODUCT(g,g)
  bet = gg/ggold
  h = -g + bet*h
  norm_g = sqrt(sum((g(1:n-1))**2))
  write(*,*) it, norm_g
  IF(norm_g < 1.d-7) exit
ENDDO
! Printing the results
print *
print *,'         It  	Norm of the gradient'
write (*,*) it, norm_g

DO ip = 1,n-1
write(12,*)v(ip),v_ex(ip),v(ip)-v_ex(ip)
enddo
END PROGRAM CG_poisson