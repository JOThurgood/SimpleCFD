program test1

  use multigrid

  implicit none

  integer, parameter :: num=selected_real_kind(p=15)
  real(num) :: pi = 4.0_num * ATAN(1.0_num)

  integer :: nx, ny, ix, iy
  integer :: nlevels

  real(num) :: dx, dy, L2, x_min, x_max, y_min, y_max

  real(num), dimension(:), allocatable :: xc, yc

  real(num), dimension(:,:), allocatable :: utrue, vtrue
  real(num), dimension(:,:), allocatable :: upol, vpol
  real(num), dimension(:,:), allocatable :: divu, phi

!  integer :: bc_xmin , bc_xmax, bc_ymin, bc_ymax
  integer, parameter :: periodic = 0
  
  ! setup a test problem 

  type(mg_input) :: input 

  nx = 16
  ny = nx
  x_min = 0.0_num
  x_max = 1.0_num
  y_min = x_min
  y_max = x_max

  print *,'Test1: Projection test, periodic, constant coeffs Poission'
  print *,'This tests the ability of MG to solve L(phi) = div(U)'
  print *,'in order to clean divergence from a polluted velocity field' 

  allocate(xc(-1:nx+2))
  allocate(yc(-1:ny+2))

  dx = (x_max - x_min) / real(nx,num)
  xc(-1) = x_min - 3.0_num * dx / 2.0_num 
  do ix = 0,nx+2
    xc(ix) = xc(ix-1) + dx
  enddo

  dy = (y_max - y_min) / real(ny,num)
  yc(-1) = y_min - 3.0_num * dy / 2.0_num 
  do iy = 0,ny+2
    yc(iy) = yc(iy-1) + dy
  enddo

  allocate(phi(-1:nx+2,-1:ny+2)) ! again, redundant ghosts but for comparison to CCAPS
  phi = 0.0_num 


  allocate(utrue(0:nx+1,0:ny+1)) !no need for ghosts in the test case
  allocate(vtrue(0:nx+1,0:ny+1))
  allocate(upol(0:nx+1,0:ny+1))
  allocate(vpol(0:nx+1,0:ny+1))
  allocate(divu(1:nx,1:ny))
 
  do iy = 0, ny+1
  do ix = 0, nx+1
    utrue(ix,iy) = - sin(1.0_num * pi * xc(ix))**2 & 
      & * sin(2.0_num * pi * yc(iy))
    vtrue(ix,iy) =   sin(1.0_num * pi * yc(iy))**2 &
      & * sin(2.0_num * pi * xc(ix))
 
    upol(ix,iy) = utrue(ix,iy) - & 
      & 0.2_num * pi * cos(2.0_num * pi * yc(iy)) &
      & * sin(2.0_num * pi * xc(ix)) 
    vpol(ix,iy) = vtrue(ix,iy) - & 
      & 0.2_num * pi * sin(2.0_num * pi * yc(iy)) & 
      & * cos(2.0_num * pi * xc(ix))
  enddo
  enddo

  do iy = 1, ny  
  do ix = 1, nx  
    divu(ix,iy) = (upol(ix+1,iy) - upol(ix-1,iy)) / 2.0_num / dx + & 
            & (vpol(ix,iy+1) - vpol(ix,iy-1)) / 2.0_num / dy
  enddo
  enddo
 
  print *, 'max divu before cleaning',maxval(divu)
  L2 = sqrt(sum(abs(utrue(1:nx,1:ny)-upol(1:nx,1:ny))**2) / real(nx*ny,num))
  print *, 'L2 (utrue vs upol) before cleaning',L2
  L2 = sqrt(sum(abs(vtrue(1:nx,1:ny)-vpol(1:nx,1:ny))**2) / real(nx*ny,num))
  print *, 'L2 (vtrue vs vpol) before cleaning',L2

  ! create the input state object

  input = mg_input(tol = 1e-12_num, nx=nx, ny = ny, dx=dx, dy=dy, f = divu, phi = phi, &
            & bc_xmin = 'periodic', bc_ymin='periodic', bc_xmax='periodic', bc_ymax = 'periodic')

  call new_mg_interface(input)
STOP

  ! solve for phi
  nlevels = -1 ! auto
  call mg_interface(f = divu, phi=phi, tol = 1e-12_num, &
                  & nx = nx, ny = ny, dx = dx, dy = dy, &
                  &  nlevels = nlevels, &
  & bc_xmin_in = periodic, bc_ymin_in = periodic, bc_xmax_in = periodic, bc_ymax_in=periodic)

  ! correct it 
   
  do iy = 1, ny
  do ix = 1, nx
    upol(ix,iy) = upol(ix,iy) - 0.5_num*(phi(ix+1,iy)-phi(ix-1,iy))/dx
    vpol(ix,iy) = vpol(ix,iy) - 0.5_num*(phi(ix,iy+1)-phi(ix,iy-1))/dy
  enddo
  enddo
   
  ! apply periodic bc to the polluted vels and then re calc divu
  upol(0,:) = upol(nx,:)
  upol(nx+1,:) = upol(1,:)
  upol(:,0) = upol(:,ny)
  upol(:,ny+1) = upol(:,1)
  vpol(0,:) = vpol(nx,:)
  vpol(nx+1,:) = vpol(1,:)
  vpol(:,0) = vpol(:,ny)
  vpol(:,ny+1) = vpol(:,1)
   
  do iy = 1, ny  
  do ix = 1, nx  
    divu(ix,iy) = (upol(ix+1,iy) - upol(ix-1,iy)) / 2.0_num / dx + & 
            & (vpol(ix,iy+1) - vpol(ix,iy-1)) / 2.0_num / dy
  enddo
  enddo
   
  print *, 'max divu after cleaning',maxval(divu)
  L2 = sqrt(sum(abs(utrue(1:nx,1:ny)-upol(1:nx,1:ny))**2) / real(nx*ny,num))
  print *, 'L2 (utrue vs upol) after cleaning',L2
  L2 = sqrt(sum(abs(vtrue(1:nx,1:ny)-vpol(1:nx,1:ny))**2) / real(nx*ny,num))
  print *, 'L2 (vtrue vs vpol) after cleaning',L2



end program test1
