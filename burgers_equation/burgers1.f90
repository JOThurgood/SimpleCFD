! solve the 1d inviscid burgers equation, in conservative form:
!
!   u_t + [0.5 * u**2 ]_x=0 
!   u=u(x,t) (user defined at t=0)
!
! using a second order finite volume predictor-corrector scheme and optional
! gradient limiters. Boundaries are (minimally-tested) zero gradient.
!
! user controled parameters in subroutines "control" and "initial_conditions"
!
! https://github.com/JOThurgood/SimpleCFD/wiki/burgers1.f90

module shared_data
  implicit none
  integer, parameter :: num=selected_real_kind(p=15)
  integer :: nx !number of cells
  integer :: ix
  integer(num) :: step, nsteps
  real(num) :: time=0.0_num, t_end, dt
  real(num) :: x_min, x_max 
  real(num) :: dx
  real(num), dimension(:), allocatable :: xb, xc
  real(num), dimension(:), allocatable :: u,u1, u1L, u1R
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

  private 

  public :: initial_setup

  contains

  subroutine initial_setup
    call control
    call setup_1dfvgrid
    call initial_conditions
    call warnings
  end subroutine initial_setup

  subroutine control
    !user control variables
    nx = 256
    x_min = 0.0_num
    x_max = 1.0_num
    t_end = 0.10_num
    CFL = 0.8_num
    nsteps = -1    !set to -1 to run till t_end
!    limiter = lim_unlimited
    limiter = lim_minmod
  end subroutine control

  subroutine initial_conditions
    !user set initial profile on a
!    call sinusoid1
    call jump
!    call gaussian
  end subroutine initial_conditions

  subroutine jump 
    real(num) :: left, right, x_diaphragm 
    x_diaphragm = 0.5 
    left = 1.0_num
    right = 2.0_num

    u = left
    do ix = -1,nx+2
      if (xc(ix) .GE. x_diaphragm) u(ix) = right 
    enddo 
  end subroutine jump 

  subroutine sinusoid1
    real(num) :: lo, hi, pi
    pi = 4.0_num * ATAN(1.0_num)
    lo = 1.0_num / 3.0_num
    hi = 2.0_num / 3.0_num
    u = 1.0_num
    do ix = -1,nx+2
      if ( (xc(ix) .ge. lo) .and. (xc(ix) .le. hi) ) then
        u(ix) = u(ix) + 0.5_num * &
          & SIN( (2.0_num * pi * (xc(ix) - 1.0_num/3.0_num)) / &
          & (1.0_num/3.0_num))
      endif
    enddo 
  end subroutine sinusoid1

  subroutine setup_1dfvgrid !set up grid and other alloctes

    !with 2 guards below
    allocate(xc(-1:nx+2))  !cell center (cc) - coutns from 1 to nx
    allocate(xb(-2:nx+2))  !cell boundary (cb) - counts from 0 to nx
    allocate(u(-1:nx+2))   !cc + 2 guards

    !no guards needed below
    allocate(u1(0:nx))  !cb - then used to calc a cc var in corrector step  
    allocate(u1L(0:nx)) !cb 
    allocate(u1R(0:nx)) !cb
  
    dx = (x_max - x_min) / REAL(nx,num)
    do ix = -2, nx+2
    xb(ix) = x_min + REAL(ix,num) * dx
    if (ix /= -2) xc(ix) = xb(ix) - dx/2.0_num
    enddo

  end subroutine setup_1dfvgrid  

  subroutine warnings
    !stub
  end subroutine warnings

  end module setup

module solver 

  use shared_data

  implicit none

  private 

  public :: update

  contains 

  subroutine update
    call boundary
    call set_dt 
    call interface_states !i_s and riem  are "predictor" step
    call riemann 
    call corrector
    time = time + dt
  end subroutine update

  subroutine boundary
  !zero-gradient (untested)
    u(0) = u(1) 
    u(-1) = u(2) 
    u(nx+1) = u(nx)
    u(nx+2) = u(nx-1)
  end subroutine boundary

  subroutine set_dt 
    dt = 1e15_num
    do ix = 0,nx
      dt = MIN(CFL * dx / abs(u(ix)), dt)
    enddo
  !    print *,'dt debug: step/selected dt', dt  
  end subroutine set_dt


  subroutine corrector
  !sign right in eq 6.17 ? maybe just due to rearangement of left and right flux
    real(num) :: fl, fr ! flux
    do ix = 1,nx !cc
      fr = 0.5_num * u1(ix) ** 2
      fl = 0.5_num * u1(ix-1) ** 2
      u(ix)  = u(ix)  - dt / dx * (fr-fl)
    enddo
  endsubroutine corrector

  subroutine interface_states
    real(num) :: slope
    do ix = 0, nx ! states on cb !careful with cc vs cb
      !left state
      if (limiter == lim_minmod) then
        slope = minmod( ( u(ix)-u(ix-1) )/dx, ( u(ix+1)-u(ix) )/dx )
      else 
        slope = ( u(ix+1) - u(ix-1) ) / 2.0_num / dx
      endif
      u1L(ix) = u(ix) + 0.5_num * dx * (1.0_num - dt * u(ix) / dx) *  slope
      !right state
      if (limiter == lim_minmod) then
        slope = minmod( ( u(ix+1)-u(ix) )/dx, ( u(ix+2)-u(ix+1) )/dx )
      else 
        slope = ( u(ix+2) - u(ix) ) / 2.0_num / dx
      endif
      u1r(ix) = u(ix+1) - 0.5_num * dx * (1.0_num + dt * u(ix) / dx) *  slope
    enddo
  end subroutine interface_states

  subroutine riemann 
    real(num) :: S, ul, ur,us
    !keep until sure about their notation
    do ix = 0,nx !a1 + interface states on cb
      ul = u1l(ix) 
      ur = u1r(ix) 
      S = 0.5_num * (ul + ur) 
      if ( ul > ur ) then ! compression
        if ( S > 0.0_num) then
          us = ul
        else if (S < 0.0_num) then
          us = ur
        else
          us = 0.0_num
        endif
      else !uniform flow or expansion 
        !(this is normal upwinding,  no assumption on sign of u) 
        if (ul > 0.0_num) then
          us = ul
        else if (ur < 0.0_num) then
          us = ur
        else 
          us = 0.0_num
        endif
      endif   

      u1(ix) = us 
   enddo

  end subroutine riemann 

  real(num) function minmod(a,b)  ! "a" is an unforunate choice of var name 
    real(num), intent(in) :: a, b !change to l and r?
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
    call execute_command_line("rm -rf xc.dat u.dat") !dont stream to existing
    open(out_unit, file="xc.dat", access="stream")
    write(out_unit) xc(1:nx)
    close(out_unit)
    open(out_unit, file="u.dat", access="stream")
    write(out_unit) u(1:nx)
    close(out_unit)
    call execute_command_line("python plot_burgers1.py")
  end subroutine do_io

end module diagnostics

program burgers1 !main driver

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

end program burgers1

