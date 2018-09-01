! solve the 1d inviscid burgers equation, in conservative form:
!
!   u_t + [0.5 * u**2 ]_x=0 
!   u=u(x,t) (user defined at t=0)
!
! using a second order finite volume predictor-corrector scheme and optional gradient limiters
! Boundaries are (minimally-tested) zero gradient.
!
! user controled parameters in subroutines "control" and "initial_conditions"
!
! Visually seems to capture correct shock formation (e.g. use IC sinusoid1) and rarefaction
! (e.g. IC jump with uright>uleft). Not done a more quantatative comparison to analytics etc 
! 
! I was a bit surprised in jump (shock) problems with no gradient limiter that any 
! overshoot / gibbs phenomena downstream is pretty minimal - inherent advantage of
! involving a R-solve in the solution scheme?
!
! outputs a .png of u(x) excluding ghosts at t_end
!
! The code is broken into more subroutines and modules than neccesary 
! for such a simple problem. The logic of the subroutines mimics a 
! more sophistic hydrodynamics code whihc would warrent such an approach 
!
! The plot uses pyplot to simplify IO significantly.
! Requires pyplot_module.f90 to be linked. 
! Easiest to place both in "src", then automagicaly build with FoBiS.py  
!
! Jonathan Thurgood , 2019-07-04 
! 

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
    t_end = 1.00_num
    CFL = 0.7_num
    nsteps = -1    !set to -1 to run till t_end
    limiter = lim_unlimited
!    limiter = lim_minmod
  end subroutine control

  subroutine initial_conditions
    !user set initial profile on a
!    call sinusoid1
    call jump
!    call gaussian
  end subroutine initial_conditions

  subroutine jump 
    real(num) :: left, right, x_diaphragm 
    x_diaphragm = 0.1 
    left = 1.0_num
    right = 0.1_num

    u = left
    do ix = -1,nx+2
      if (xc(ix) .GE. 0.1_num) u(ix) = right 
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
        u(ix) = u(ix) + 0.5_num * SIN( (2.0_num * pi * (xc(ix) - 1.0_num/3.0_num)) / &
          & (1.0_num/3.0_num))
      endif
    enddo 
  end subroutine sinusoid1

  subroutine setup_1dfvgrid
  !setup the grid and also allocate other quantities 
    allocate(xc(-1:nx+2))  !cell center (cc) - counts from 1 to nx, with 2 guards
    allocate(xb(-2:nx+2))  !cell boundary (cb) - counts from 0 to nx, with 2 guards
    allocate(u(-1:nx+2))   !cc + 2 guards
    allocate(u1(0:nx))  !cb - then used to calc a cc var in corrector step -  no guards needed 
    allocate(u1L(0:nx)) !cb 
    allocate(u1R(0:nx)) !cb
  
    dx = (x_max - x_min) / REAL(nx,num)
    do ix = -2, nx+2
    xb(ix) = REAL(ix,num) * dx
    if (ix /= -2) xc(ix) = xb(ix) - dx/2.0_num
    enddo

  end subroutine setup_1dfvgrid  

  subroutine warnings
! might want to stop on some condition e.g.
!    if (condition) then
!      print *,'ic have set <bad condition> - incompatible with assumptions - STOP'
!      STOP
!    endif
  end subroutine warnings

  end module setup

  module solver 

  use shared_data

  implicit none

  private 

  public :: update

  contains !number of subroutines obviously a bit of over-kill for this problem

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
    !be careful with cc vs cb here
    do ix = 0, nx ! states on cb
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
    real(num), intent(in) :: a, b !left and right difference so maybe change to l and r
!    real(num), intent(out) :: minmod
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

  use pyplot_module, only : pyplot
  use shared_data

  implicit none 

  contains

  subroutine plot1d 
    use, intrinsic :: iso_fortran_env, only : wp => real64
    integer :: istat
    type(pyplot) :: plt

    call plt%initialize(grid=.true.,xlabel='x',figsize=[20,10],&
                        title='plot test',legend=.true.,axis_equal=.false.,&
                        tight_layout=.true.)
    call plt%add_plot(xc(1:nx),u(1:nx),label='u(x)',linestyle='b-o',markersize=5,linewidth=2,istat=istat)
!    call plt%add_plot(xc,u,label='u(x)',linestyle='b-o',markersize=5,linewidth=2,istat=istat)
    call plt%savefig('output.png', pyfile='plottest.py',istat=istat)!
  endsubroutine plot1d

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

  call plot1d

  print *, 'Done in',step,'steps', 'with CFL',CFL

end program burgers1

