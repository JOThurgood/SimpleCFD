! test sln of 2d diffusion equation
! under Crank-Nicholson discritisation
! and MMG 

program test5

  use multigrid

  implicit none

  integer, parameter :: num = selected_real_kind(p=15)
  integer :: out_unit = 10

  integer :: nx = 256
  integer :: ny = 256

  real(num) :: x_min = 0.0_num
  real(num) :: x_max = 1.0_num
  real(num) :: y_min = 0.0_num
  real(num) :: y_max = 1.0_num

  real(num) :: time = 0.0_num
  real(num) :: t_end = 0.02_num

  real(num) :: k = 1.0_num
  real(num) :: phi2 = 2.0_num !max val
  real(num) :: phi1 = 1.0_num !baseline

  real(num) :: CFL = 2.0_num

  real(num) :: tol = 1e-12_num

  integer :: tstep = 0
  integer :: rstep = 0
  integer :: ix, iy
  real(num) :: dx, dy
  real(num) :: r,xx,yy

  real(num) :: dt
  real(num), dimension(:), allocatable :: x, y
  real(num), dimension(:,:), allocatable :: phi, phi_an, f

  real(num) :: L_phi
  real(num) :: alpha, beta

  real(num) :: L2

  type(mg_input) :: input

  real(num) :: start
  real(num) :: finish
   
  allocate(x(0:nx+1))  ! cell center (cc) - counts from 1 to nx
  allocate(y(0:ny+1))  ! cell center (cc) - counts from 1 to ny
                        ! ONE GUARD CELL FOR THIS PROBLEM
  dx = (x_max - x_min) / REAL(nx,num)
  do ix = 1, nx+1
    x(ix) = x(ix-1) + dx
  end do
   
  dy = (y_max - y_min) / REAL(ny,num)
  y(0) = y_min - dy/2.0_num
  do iy = 1, ny+1
  y(iy) = y(iy-1) + dy
  end do
   
  if (dy /= dx) then
    print *, 'Warning: dy /= dx set. Has not yet being generalised.'
    print *, 'Terminating'
    STOP
  end if 
   

  allocate(phi(-1:nx+2,-1:ny+2))
  allocate(phi_an(-1:nx+2,-2:ny+2))

  do ix = 0, nx + 1
  do iy = 0, ny + 1
    xx = x(ix) - 0.5_num
    yy = y(iy) - 0.5_num
    r = sqrt(xx**2 + yy**2)
    phi(ix,iy) = (phi2-phi1) * exp(-0.25_num * r**2 / k / 0.001_num) + phi1
  enddo
  enddo

!  GOTO 1 ! uncomment to just dump initial condiiton and check it, bypassing the loop

  allocate(f(1:nx,1:ny))

  ! evolution loop

  time = 0.0_num
  call cpu_time(start)
  do
    tstep = tstep + 1
    if (time > t_end) exit

    dt = CFL * dx**2 / k

    alpha = 1.0_num
    beta = dt * k / 2.0_num

    do ix = 1, nx
    do iy = 1, ny

      ! function based on laplacian of phi at the previous time level
      ! - dont update during relaxation

      L_phi = (phi(ix+1,iy) - 2.0_num*phi(ix,iy) + phi(ix-1,iy)) / dx**2 + &
        & (phi(ix,iy+1) - 2.0_num*phi(ix,iy) + phi(ix,iy-1)) / dy**2 
       
      f(ix,iy) = phi(ix,iy) + 0.5_num * dt * k * L_phi

    enddo
    enddo

    input = mg_input(tol = tol, nx=nx, ny = ny, dx=dx, dy=dy, f = f, phi = phi, &
            & bc_xmin = 'zero_gradient', bc_ymin='zero_gradient', &
            & bc_xmax='zero_gradient', bc_ymax = 'zero_gradient', &
!            & deallocate_after = .true., &
            & quiet = .true., & 
            & const_helmholtz = .true., ch_alpha = alpha, ch_beta = beta)
    call mg_interface(input)
    phi = input%phi

    print *,'time', time,'tstep',tstep, 'done'

    time = time + dt
  enddo

1 CONTINUE

  call cpu_time(finish)
  print '(" test 5 total cpu_time: ",f20.3," seconds.")',finish-start


  !the analytical solution
  do ix = 1,nx
  do iy = 1,ny
    xx = x(ix) - 0.5_num
    yy = y(iy) - 0.5_num
    r = sqrt(xx**2 + yy**2)
    phi_an(ix,iy) = (phi2-phi1) * (0.001_num / (time+0.001_num)) &
      & * exp(-0.25_num * r**2 / k / (0.001_num + time)) + phi1
  enddo
  enddo

  !calculate the L2 norm against the analytical solution

  L2 = 0.0_num

  do ix = 1, nx
  do iy = 1, ny
    L2 = L2 + abs(phi_an(ix,iy)-phi(ix,iy))**2 
  enddo
  enddo
  L2 = sqrt( L2 / REAL(nx*ny,num))

  print *,'L2 error norm agains analtical solution is',L2



  call execute_command_line("rm -rf *.dat test5_*.png")
    !dont want to stream to existing files
 
  open(out_unit, file="x.dat", access="stream")
  write(out_unit) x(1:nx)
  close(out_unit)
 
  open(out_unit, file="y.dat", access="stream")
  write(out_unit) y(1:ny)
  close(out_unit)
 
  open(out_unit, file="phi.dat", access="stream")
  write(out_unit) phi(1:nx,1:ny)
  close(out_unit)

  open(out_unit, file="phi_an.dat", access="stream")
  write(out_unit) phi_an(1:nx,1:ny)
  close(out_unit)

  call execute_command_line("python test5_plot.py")
  call execute_command_line("rm -rf *.dat")




end program test5
