! solve the 2d linear advection equation, in conservative form:
!
!   a_t + [f(a)]_x + [g(a)]_y=0 
!   f(a) = u*a
!   g(a) = v*a 
!   a=a(x,y,t) (user defined) 
!   u=constant (user defined)
!   v=constant (user defined)
!
!  using either
!   - dimensionally split 2nd order finite volume predictor-corrector 
!   - unsplit with "Corner Transport Upwind" (CTU) method (coming soon)
!  with a choice of gradient limiters
!
! https://github.com/JOThurgood/SimpleCFD/wiki/advect_2d.f90
!

module shared_data

  implicit none

  integer, parameter :: num=selected_real_kind(p=15)
  integer :: nx, ny !number of cells
  integer :: ix,iy
  integer(num) :: step, nsteps
  real(num) :: u, v
  real(num) :: time=0.0_num, t_end, dt
  real(num) :: x_min, x_max 
  real(num) :: y_min, y_max 
  real(num) :: dx, dy
  real(num), dimension(:), allocatable :: xb, xc
  real(num), dimension(:), allocatable :: yb, yc
  real(num), dimension(:,:), allocatable :: a,a1x,a1y, a1Lx, a1Rx, a1Ly,a1Ry
    !solution and states at half time step indicated by 1
  real(num) :: CFL !cfl modifier

  !keyword codes
  integer :: limiter
  integer, parameter :: lim_unlimited = 1, lim_minmod = 2
  integer :: method
  integer, parameter :: method_split = 1, method_unsplit = 2

  contains

end module shared_data

module setup

  use shared_data

  implicit none

  contains

  subroutine initial_setup
    call control
    call setup_2dfvgrid
    call initial_conditions
  end subroutine initial_setup

  subroutine control
    !user control variables
    u = 1.0_num
    v = 1.0_num
    nx = 64 
    ny = 64  
    x_min = 0.0_num
    x_max = 1.0_num
    y_min = x_min
    y_max = x_max
    t_end = 1.0_num
    CFL = 0.7_num !cfl multiplier <= 1
    nsteps = -1  !set to -1 to run till t_end
    limiter = lim_minmod
!    limiter = lim_unlimited
    method = method_split 
  end subroutine control

  subroutine initial_conditions
    !user set initial profile on a
    a = 0.0_num
!   call x_top_hat
!    call y_top_hat
    call square 
!    call gaussian
  end subroutine initial_conditions

  subroutine x_top_hat
    a = 0.0_num
    do ix = 0,nx+1
      if ((xc(ix) .ge. 1.0_num/3.0_num) .and. (xc(ix) .le. 2.0_num/3.0_num)) &
        & a(ix,:) = 1.0_num
    enddo
  end subroutine x_top_hat

  subroutine square 
    real(num) :: lo,hi
    lo = 1.0_num/3.0_num
    hi = 2.0_num*lo
    a = 0.0_num
    do iy = 0,ny+1
    do ix = 0,nx+1
      if ((xc(ix) .ge. lo) & 
        & .and. (xc(ix) .le. hi) &
        & .and. (yc(iy) .ge. lo) &
        & .and. (yc(iy) .le. hi) ) a(ix,iy) = 1.0_num
    enddo
    enddo
  end subroutine square

  subroutine y_top_hat
    a = 0.0_num
    do iy = 0,ny+1
      if ((yc(iy) .ge. 1.0_num/3.0_num) .and. (yc(iy) .le. 2.0_num/3.0_num)) &
        & a(:,iy) = 1.0_num
    enddo
  end subroutine y_top_hat

  subroutine gaussian
    print *, 'gaussian currently a stub'
  end subroutine gaussian

  subroutine setup_2dfvgrid !setup the grid and allocate other arrays 

    allocate(xc(-1:nx+2))  !cell center (cc) - counts from 1 to nx
    allocate(xb(-2:nx+2))  !cell boundary (cb) - counts from 0 to nx
    allocate(yc(-1:ny+2))  !cell center (cc) - counts from 1 to ny
    allocate(yb(-2:ny+2))  !cell boundary (cb) - counts from 0 to ny
                           !all above have two guards  either side

    allocate(a(-1:nx+2,-1:ny+2))   !cc + 2 guards
    allocate(a1x(0:nx,1:ny))  !cb,cc - no guards needed 
    allocate(a1Lx(0:nx,1:ny)) !cb,cc 
    allocate(a1Rx(0:nx,1:ny)) !cb,cc
    allocate(a1y(1:nx,0:ny))  !cc,cb - no guards needed 
    allocate(a1Ly(1:nx,0:ny)) !cc,cb 
    allocate(a1Ry(1:nx,0:ny)) !cc,cb

    dx = (x_max - x_min) / REAL(nx,num)
    do ix = -2, nx+2
      xb(ix) = REAL(ix,num) * dx 
      if (ix /= -2) xc(ix) = xb(ix) - dx/2.0_num
    enddo

    dy = (y_max - y_min) / REAL(ny,num)
    do iy = -2, ny+2
      yb(iy) = REAL(iy,num) * dy 
      if (iy /= -2) yc(iy) = yb(iy) - dy/2.0_num
    enddo

  end subroutine setup_2dfvgrid  

end module setup

module split_solver 

  use shared_data

  implicit none

  private

  public :: solve_split

  contains 

  subroutine solve_split

    call set_dt 

    if (MOD(step,2) == 1) then !strang splitting 
      call x_bcs
      call x_interface_state
      call x_reimann
      call x_update
      call y_bcs
      call y_interface_state
      call y_reimann
      call y_update
    else
      call y_bcs
      call y_interface_state
      call y_reimann
      call y_update
      call x_bcs
      call x_interface_state
      call x_reimann
      call x_update
    endif

    time = time + dt

    !if plotting ghost cells may want to call bcs before plot
    !to update (order here will have them one time before?)
    !no effect on sln / correct for what you intend
  end subroutine solve_split

  subroutine set_dt
    real(num) ::dtx,dty
    dtx = CFL * dx / abs(u)
    dty = CFL * dy / abs(v)
    dt = MIN(dtx,dty)
    if ((step==1) .and. (abs(u-v) .gt. 1e-2_num)) then
      print *,'u not ~= v detected ; even if CFL=1 perfect translation not poss'  
      print *,'effective CFL in x drn',dt/dtx * CFL
      print *,'effective CFL in y drn',dt/dty * CFL
    endif
!    if (step==1) print *,'effective CFL in x drn',dt/dtx
!    if (step==1) print *,'effective CFL in y drn',dt/dty
  end subroutine set_dt

  subroutine x_bcs !periodic
    a(-1,:) = a(nx-1,:)
    a(0,:) = a(nx,:)
    a(nx+1,:) = a(1,:)
    a(nx+2,:) = a(2,:)
  end subroutine x_bcs

  subroutine y_bcs !periodic
    a(:,-1) = a(:,ny-1)
    a(:,0) = a(:,ny)
    a(:,ny+1) = a(:,1)
    a(:,ny+2) = a(:,2)
  end subroutine y_bcs

  subroutine x_interface_state
    !be careful with cc vs cb here
    real(num) dadx
    do iy = 1, ny !inefficient, refactor to vector operation
    do ix = 0, nx 
      !left state
      if (limiter == lim_minmod) then
        dadx = minmod( ( a(ix,iy)-a(ix-1,iy) ) / dx, ( a(ix+1,iy)-a(ix,iy) ) / dx )
      else 
        dadx = ( a(ix+1,iy) - a(ix-1,iy) ) / 2.0_num / dx
      endif
      a1Lx(ix,iy) = a(ix,iy) + 0.5_num * dx * (1.0_num - dt * u / dx) *  dadx
      !right state
      if (limiter == lim_minmod) then
        dadx = minmod( ( a(ix+1,iy)-a(ix,iy) )/dx , ( a(ix+2,iy)-a(ix+1,iy) ) / dx )
      else
        dadx = ( a(ix+2,iy) - a(ix,iy) ) / 2.0_num / dx
      endif

      a1Rx(ix,iy) = a(ix+1,iy) - 0.5_num * dx * (1.0_num + dt * u / dx) *  dadx
    enddo
    enddo
  end subroutine x_interface_state

  subroutine x_reimann
    !worlds simplest reimann problem
    if (u < 0.0_num) then
      a1x = a1Rx
    else if (u > 0.0_num) then
      a1x = a1Lx
    endif
  end subroutine x_reimann


  subroutine x_update
    do ix = 1,nx
      a(ix,1:ny)  = a(ix,1:ny)  - dt / dx * (u*a1x(ix,1:ny)-u*a1x(ix-1,1:ny))
    enddo
  endsubroutine x_update

  subroutine y_interface_state
    !consider a1L as "bottom" and a1R as "top" conceptually
    ! (or most negative -> most positive y)
    !be careful with cc vs cb here
    real(num) dady
    do ix = 1, nx !refactor to vector operation
    do iy = 0, ny 
      !left state
      if (limiter == lim_minmod) then
        dady = minmod( ( a(ix,iy)-a(ix,iy-1) ) / dy, ( a(ix,iy+1)-a(ix,iy) ) / dy )
      else 
        dady = ( a(ix,iy+1) - a(ix,iy-1) ) / 2.0_num / dy
      endif
      a1Ly(ix,iy) = a(ix,iy) + 0.5_num * dy * (1.0_num - dt * v / dy) *  dady
      !right state
      if (limiter == lim_minmod) then
        dady = minmod( ( a(ix,iy+1)-a(ix,iy) )/dy , ( a(ix,iy+2)-a(ix,iy+1) ) / dy )
      else
        dady = ( a(ix,iy+2) - a(ix,iy) ) / 2.0_num / dy
      endif
      a1Ry(ix,iy) = a(ix,iy+1) - 0.5_num * dy * (1.0_num + dt * v / dy) *  dady
    enddo
    enddo
  end subroutine y_interface_state

  subroutine y_reimann
    if (v < 0.0_num) then
      a1y = a1Ry
    else if (v > 0.0_num) then
      a1y = a1Ly
    endif
  end subroutine y_reimann

  subroutine y_update
    do iy = 1,ny
      a(1:nx,iy)  = a(1:nx,iy)  - dt / dy * (v*a1y(1:nx,iy)-v*a1y(1:nx,iy-1))
    enddo
  endsubroutine y_update

  real(num) function minmod(a,b)  ! refactor - "a" is an unfortunate choice
    real(num), intent(in) :: a, b 
    if ( (abs(a) < abs(b)) .and. (a*b > 0) ) then
      minmod = a
    else if ( (abs(a) > abs(b)) .and. (a*b > 0) ) then
      minmod = b
    else 
      minmod = 0.0_num
    endif 
  end function minmod

end module split_solver

module diagnostics

  use shared_data

  implicit none 

  contains

  subroutine do_io

    integer :: out_unit =10

    call execute_command_line("rm -rf xc.dat yc.dat a.dat") 
      !dont want to stream to existing files

    open(out_unit, file="xc.dat", access="stream")
    write(out_unit) xc
    close(out_unit)

    open(out_unit, file="yc.dat", access="stream")
    write(out_unit) yc
    close(out_unit)

    open(out_unit, file="a.dat", access="stream")
    write(out_unit) a
    close(out_unit)
 
    call execute_command_line("python plot_advect_2d.py")


  end subroutine do_io

end module diagnostics

program advection2d !main driver

  use shared_data
  use setup
  use split_solver
  use diagnostics

  implicit none

  call initial_setup

  do 
   if ((step >= nsteps .and. nsteps >= 0) .or. (time >= t_end)) exit
    step = step + 1
    if (method == method_split) call solve_split
  enddo

  call do_io

  print *, 'Done in',step,'steps', 'with CFL',CFL

end program advection2d
