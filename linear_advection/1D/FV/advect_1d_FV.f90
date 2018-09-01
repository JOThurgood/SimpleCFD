! solve the linear advection equation, in conservative form:
!
!   a_t + [f(a)]_x=0 
!   f(a) = u*a 
!   a=a(x,t) (user defined) 
!   u=constant (user defined)
!
! using a second order finite volume predictor-corrector scheme and optional gradient limiters
! Boundaries are periodic.
!
! user controled parameters in subroutines "control" and "initial_conditions"
!
! When CFL<1, if limiters off shows spurious oscillations in addition to numerical diffusion
! (as opposed to just diffusion in first order method w/ appropriate upwinding?)
! This highlights issues 
! relating to Godunov's theorem and that in order to have a second-order 
! accurate solution that is monotonic (doesn't introduce new minima and maxima)
! that the method itself must be nonlinear. Easiest to see for tophat initial condition.
! 
! Outputs "output.png" of the solution at t_end. Best to set u and x_min and x_max
! (e.g. u=1, x_min=0,x_max=1)
! so that t = 1 corresponds to one crossing time (recycled through periodic boundary)
!
! The code is broken into more subroutines and modules than neccesary 
! for such a simple problem. The logic of the subroutines mimics a 
! more sophistic hydrodynamics code whihc would warrent such an approach 
! e.g., set_dt is trivial herem but would be more involved a variable speed problem
!
! The plot uses pyplot to simplify IO significantly.
! Requires pyplot_module.f90 to be linked. 
! Easiest to place both in "src", then automagicaly build with FoBis.py  
!
! Jonathan Thurgood , 2019-06-28 
! 

module shared_data
  implicit none
  integer, parameter :: num=selected_real_kind(p=15)
  integer :: nx !number of cells
  integer :: ix
  integer(num) :: step, nsteps
  real(num) :: u
  real(num) :: time=0.0_num, t_end, dt
  real(num) :: x_min, x_max 
  real(num) :: dxb, dxc !redundant, dx=dxb=dxc=constant here
  real(num), dimension(:), allocatable :: xb, xc
  real(num), dimension(:), allocatable :: a,a1, a1L, a1R
    !solution and states at half time step indicated by 1
  real(num) :: CFL !cfl modifier
  !keyword codes
  integer :: limiter
  integer, parameter :: lim_unlimited = 1, lim_minmod = 2
  contains
end module shared_data

module setup

  use shared_data

  implicit none

  contains

  subroutine initial_setup
    call control
    call setup_1dfvgrid
    call initial_conditions
  end subroutine initial_setup

  subroutine control !user control variables
    u = 1.0_num
    nx = 64 
    x_min = 0.0_num
    x_max = 1.0_num
    t_end = 1.0_num
    CFL = 0.7_num
    nsteps = -1   !set to -1 to run till t_end
    limiter = lim_unlimited !lim_unlimited, lim_minmod
  end subroutine control

  subroutine initial_conditions !choose by commenting or add your own
    call top_hat
    !call gaussian
  end subroutine initial_conditions

  subroutine top_hat
    a = 0.0_num
    do ix = 0,nx+1 
      if ((xc(ix) .ge. 1.0_num/3.0_num) .and. (xc(ix) .le. 2.0_num/3.0_num)) a(ix) = 1.0_num
    enddo
  end subroutine top_hat

  subroutine gaussian
    real(num) :: fwhm, sigma,x0,xx
    a = 0.0_num
    fwhm = 0.2_num
    x0 = 0.5_num 
    sigma = fwhm / 2.35482_num 
    do ix = 0,nx+1 
      xx = xc(ix) - x0
      a(ix) = 1.0_num * exp(-xx**2 / 2.0_num / sigma**2 ) 
    enddo
  end subroutine gaussian

  subroutine setup_1dfvgrid !setup grid  + allocate other vars
    allocate(xc(-1:nx+2))  !cell center (cc) - counts from 1 to nx, with 2 guard
    allocate(xb(-2:nx+2))  !cell boundary (cb) - counts from 0 to nx, with 2 gds
    allocate(a(-1:nx+2))   !cc + 2 guards
    allocate(a1(0:nx))  !cb - then used to calc a cc var in corrector step 
    allocate(a1L(0:nx)) !cb   (no guards needed)
    allocate(a1R(0:nx)) !cb
    dxb = (x_max - x_min) / REAL(nx,num)
    dxc = dxb !redundant for uniform grids dxb=dxc=dx=constant
    do ix = -2, nx+2
      xb(ix) = REAL(ix,num) * dxb
      if (ix /= -2) xc(ix) = xb(ix) - dxc/2.0_num
    enddo
  end subroutine setup_1dfvgrid  

end module setup

module solver 

  use shared_data

  implicit none

  contains

  subroutine update
    call boundary
    call set_dt 
    call predictor
    call corrector
    time = time + dt
  end subroutine update

  subroutine boundary !periodic
    a(-1) = a(nx-1)
    a(0) = a(nx)
    a(nx+1) = a(1)
    a(nx+2) = a(2) 
  end subroutine boundary

  subroutine set_dt 
    dt = CFL * dxc / abs(u) 
  end subroutine set_dt

  subroutine predictor
    call interface_states 
    call reimann
  end subroutine predictor

  subroutine corrector
    do ix = 1,nx 
      a(ix)  = a(ix)  - dt / dxb * (u*a1(ix)-u*a1(ix-1))
    enddo
  endsubroutine corrector

  subroutine interface_states
    !be careful with cc vs cb here
    real(num) dadx
    do ix = 0, nx ! states on boundary defined at xb 

      !left state
      if (limiter == lim_minmod) then
        dadx = minmod( ( a(ix)-a(ix-1) )/dxb, ( a(ix+1)-a(ix) )/dxb )
      else 
        dadx = ( a(ix+1) - a(ix-1) ) / 2.0_num / dxb
      endif
      a1L(ix) = a(ix) + 0.5_num * dxb * (1.0_num - dt * u / dxb) *  dadx

      !right state
      if (limiter == lim_minmod) then
        dadx = minmod( ( a(ix+1)-a(ix) )/dxb, ( a(ix+2)-a(ix+1) )/dxb )
      else 
        dadx = ( a(ix+2) - a(ix) ) / 2.0_num / dxb
      endif

      a1R(ix) = a(ix+1) - 0.5_num * dxb * (1.0_num + dt * u / dxb) *  dadx
    enddo
  end subroutine interface_states

  subroutine reimann 
    !worlds simplest Reimann problem
    if (u < 0.0_num) then 
      a1 = a1R
    else if (u > 0.0_num) then
      a1 = a1L
    endif 
  end subroutine reimann 

  real(num) function minmod(a,b)  ! "a" is an unforunate choice of var name!
    real(num), intent(in) :: a, b !left and right difference so change to l+r?
    if ( (abs(a) < abs(b)) .and. (a*b > 0) )then
      minmod = a
    else if ( (abs(a) > abs(b)) .and. (a*b > 0) ) then
      minmod = b
    else 
      minmod = 0.0_num
    endif 
  end function minmod

end module solver

module diagnostics

  use shared_data

  implicit none 

  contains

  subroutine do_io
    integer :: out_unit =10
    open(out_unit, file="xc.dat", access="stream")
    write(out_unit) xc
    close(out_unit)
    open(out_unit, file="a.dat", access="stream")
    write(out_unit) a
    close(out_unit)
    call execute_command_line("python plot_advect_1d_FV.py")
  end subroutine do_io

end module diagnostics

program advection2 !main driver

  use shared_data
  use setup
  use solver
  use diagnostics

  implicit none

  call initial_setup

  do 
   if ((step >= nsteps .and. nsteps >= 0) .or. (time >= t_end)) exit
    step = step + 1
    call update    
  enddo

  call do_io

  print *, 'Done in',step,'steps', 'with CFL',CFL

end program advection2

