! Solve 
! phi_t = k phi_xx  
! using explicit method
! 2nd order centred in space, 1st order fwd in time

program explicit1d

  implicit none

  integer, parameter :: num=selected_real_kind(p=15) !p=6 or 15
  real(num) :: x_min, x_max, x_c, dx
  real(num) :: time, t_end , dt, dtk
  integer :: step, nsteps
  integer :: ix
  integer :: nx 
  real(num) :: C, k, phi1, phi2, t0
  real(num), allocatable, dimension(:) :: phi, x, phi_new, phi_analytic
  integer :: out_unit = 10 

  !user control
  t_end = 0.008_num
  nx = 16 
  x_min = 0.0_num
  x_max = 1.0_num
  x_c = (x_max - x_min) / 2.0_num
  nsteps = 100 !max number of steps
  C = 0.8_num ! CFL-like factor in time step, le 1 required

  !initial condition optin

  k = 1.0 ! diffusion coefficient
  t0 = 0.001_num !constant in initial condition, NOT the initial time
  time = 0.0_num !initial time
  phi2 = 2.0_num 
  phi1 = 1.0_num


  ! initialise the arrays and fill grid, 1 ghost cell needed

  allocate(x(0:nx+1))
  allocate(phi(0:nx+1))
  allocate(phi_new(0:nx+1))

  dx = (x_max - x_min) / REAL(nx-1,num) 
  x(0) = x_min - dx
  do ix = 1, nx+1
    x(ix) = x(ix-1) + dx
  enddo

  ! set the initial condition on phi
  phi = 0.0_num
  do ix = 1, nx 
    phi(ix) = (phi2-phi1)  * (t0/(time + t0))**0.5_num * &
      & exp(-(x(ix)-x_c)**2 / 4.0_num / k / (time+t0)) + phi1
  enddo 


  !main loop 

  step = 0
  dt = 0.5_num * C * dx**2 / k 
  dtk = dt * k 
  do
    if ((step >= nsteps .and. nsteps >= 0) .or. (time >= t_end)) exit

    step = step+1
    time = time + dt
    print *, step,time

    !set bc (simple Dirichlet)
    phi(0) = phi1
    phi(nx+1) = phi1

    !loop through domain
    do ix = 1, nx 
      phi_new(ix) = phi(ix) + & 
        & dtk * (phi(ix+1) - 2.0_num * phi(ix) + phi(ix-1))/dx**2 
    enddo !ix = 1, nx

    phi = phi_new

  enddo

  ! Output

  call execute_command_line("rm -rf x*.dat phi*.dat")
     !^ dont stream to existing *dat files

  ! output numerical solution
  open(out_unit, file="x.dat", access="stream")
  write(out_unit) x(1:nx)
  close(out_unit)
  open(out_unit, file="phi.dat", access="stream")
  write(out_unit) phi(1:nx)
  close(out_unit)

  ! calculate and output analytical sln and corresponding grid,
  ! at several times the resolution to better-sample the continuous sln 

  if (nx <= 2**14) then
    deallocate(x)
    nx = nx * 16 
    allocate(x(0:nx+1))
    dx = (x_max - x_min) / REAL(nx-1,num) 
    x(0) = x_min - dx
    do ix = 1, nx+1
      x(ix) = x(ix-1) + dx
    enddo
  endif 

  allocate(phi_analytic(0:nx+1))
  do ix = 0,nx+1
    phi_analytic(ix) = (phi2-phi1)  * (t0/(time + t0))**0.5_num * &
      & exp(-(x(ix)-x_c)**2 / 4.0_num / k / (time+t0)) + phi1
  enddo

  open(out_unit, file="x_analytic.dat", access="stream")
  write(out_unit) x(1:nx)
  close(out_unit)

  open(out_unit, file="phi_analytic.dat", access="stream")
  write(out_unit) phi_analytic(1:nx)
  close(out_unit)
 
  call execute_command_line("python plot1.py")


end program explicit1d


