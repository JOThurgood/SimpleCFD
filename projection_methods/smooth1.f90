! Simple illustration of projection / test
! Prescribe a divergence free velocity field, 
! polute it with prescribed divergence, then
! use a relaxation method to reconstruct 
! the initial divergence free field.

! See e.g. 

program smooth1

  implicit none

  ! declarations

  integer, parameter :: num = selected_real_kind(p=15)
  real(num), parameter :: pi = 4.0_num * ATAN(1.0_num)
  integer :: out_unit = 10 

  integer :: nx, ny 
  integer :: ix,iy
  integer(num) :: step, nsteps 

  real(num) :: u, v
  real(num) :: time = 0.0_num, t_end, dt
  real(num) :: x_min, x_max 
  real(num) :: y_min, y_max 
  real(num) :: dx, dy

  real(num), dimension(:), allocatable :: x, y

  real(num), dimension(:,:), allocatable :: u_orig !origional 
  real(num), dimension(:,:), allocatable :: u_pol  !polluted 
  real(num), dimension(:,:), allocatable :: u_reconstructed
  real(num), dimension(:,:), allocatable :: v_orig !origional 
  real(num), dimension(:,:), allocatable :: v_pol  !polluted 
  real(num), dimension(:,:), allocatable :: v_reconstructed

  ! simulation / domain / grid control parameters

  nsteps = 1
  t_end = 1
  nx = 64
  ny = nx ! currently hardcoded for dx = dy later in the code
  x_min = 0.0_num 
  x_max = 1.0_num
  y_min = x_min !
  y_max = x_max

  ! array allocation and initialisation

  ! -- Grid 

  allocate(x(0:nx+1))  ! cell center (cc) - counts from 1 to nx
  allocate(y(0:ny+1))  ! cell center (cc) - counts from 1 to ny
                        ! ONE GUARD CELL FOR THIS PROBLEM

  dx = (x_max - x_min) / REAL(nx-1,num)
  x(0) = x_min - dx
  do ix = 1, nx+1
    x(ix) = x(ix-1) + dx
  enddo

  dy = (y_may - y_min) / REAL(ny-1,num)
  y(0) = y_min - dy
  do iy = 1, ny+1
    y(iy) = y(iy-1) + dy
  enddo

  if (dy /= dx) then
    print *, 'Warning: dy /= dx set. Has not yet being generalised.'
    print *, 'Terminating'
    STOP
  endif 
 
  ! -- Velocity 

  allocate(u_orig(0,nx+1,0:ny+1))
  allocate(u_pol(0,nx+1,0:ny+1))
  allocate(u_reconstructed(0,nx+1,0:ny+1))
  allocate(v_orig(0,nx+1,0:ny+1))
  allocate(v_pol(0,nx+1,0:ny+1))
  allocate(v_reconstructed(0,nx+1,0:ny+1))
 

  do iy = 0, ny + 1
  do ix = 0, nx + 1
    u_orig(ix,iy) = - sin(2.0_num * pi * x(ix))**2.0_num & 
      & * sin(2.0_num * pi * y(iy))
    v_orig(ix,iy) =   sin(2.0_num * pi * y(iy))**2.0_num &
      & * sin(2.0_num * pi * x(ix))
    u_pol(ix,iy) = u_orig(ix,iy) + 
    v_pol(ix,iy) = v_orig(ix,iy) + 
    
  enddo
  enddo  

 
  ! reconstruct the divergence free field

  ! generate output 
 
  !call execute_command_line("rm -rf *.dat") 
    !dont want to stream to existing files
  
  open(out_unit, file="xc.dat", access="stream")
  write(out_unit) xc
  close(out_unit)

  open(out_unit, file="yc.dat", access="stream")
  write(out_unit) yc
  close(out_unit)
  
  open(out_unit, file="u_orig.dat", access="stream")
  write(out_unit) u_orig
  close(out_unit)
  

end
