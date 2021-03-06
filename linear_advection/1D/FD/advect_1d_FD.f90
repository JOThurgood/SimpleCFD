!  advect_1d_FD.f90
!
!  Solve a_t + u a_x = 0 with finite-difference
!  with either FTCS or upwinding. 
!
!  https://github.com/JOThurgood/SimpleCFD/wiki

module shared_data

  implicit none

  integer, parameter  :: num=selected_real_kind(p=15)
  integer (num)  :: nx, ix
  integer (num) :: step = 0,nsteps
  real(num) :: u
  real (num) :: x_min, x_max, dxb
  real (num), dimension(:), allocatable :: a
  real(num), dimension(:), allocatable :: xb
  real (num) :: t_end, time = 0.0_num, dt, cfl 
  logical :: ftcs,upwind, verbose
  !anticipating staggered grid in later codes
  !so 'c' for cell center, 'b' for boundary 

  contains

end module shared_data

module setup

  use shared_data

  implicit none

  contains

  subroutine user_control
    nx = 64 
    x_min = 0.0_num 
    x_max = 1.0_num
    t_end = 1.0_num
    u = 1.0_num !char speed of linear advection equation
    cfl = 1.0_num !cfl number
    nsteps = -1 !<0 to run to t_end 
    ftcs = .false.
    upwind = .true.
    verbose = .true.
  end subroutine user_control 

  subroutine initial_conditions
    ! step function
    a = 0.0_num
    do ix = -1,nx+1 
      if ((xb(ix) >= 1.0_num/3.0_num) .and. (xb(ix) <= 2.0_num / 3.0_num))  then
        a(ix) = 1.0_num
      endif
    end do 
  endsubroutine initial_conditions

  subroutine setup_1d_fd_grid
    ! setup 1d finite difference grid for var a of size nx 
    ! plus one ghost either side
    ! here 0 is leftmost, nx is right, -1 and nx+1 are ghosts
    ! and allocate "a" etc

    allocate(xb (-1:nx+1) )
    allocate(a (-1:nx+1) )

    dxb = (x_max - x_min) / real(nx,num) 
    do ix = -1,nx+1 
      xb(ix) = x_min + real(ix,num) * dxb
    end do 
  endsubroutine setup_1d_fd_grid 

  subroutine init
    call user_control
    call setup_1d_fd_grid
    call initial_conditions
  end subroutine init

end module setup

module solver

  use shared_data

  implicit none

  contains

  subroutine update
    dt = cfl * dxb / u
    call periodic_bc

    if (ftcs) then
!    do ix=0,nx 
!      a(ix) =  a(ix) - cfl * (a(ix+1) - a(ix-1)) / 2.0_num
!    enddo
! NB: above will not work - need previous state so any direction will break it
      a = a - cfl * ( cshift(a,1) - cshift(a,-1)) / 2.0_num
    endif 

    if (upwind) then
    ! upwind - must either store old solution separately,
    ! or fill from right to left
      do ix=nx,0,-1
        a(ix) =  a(ix) - cfl * (a(ix) - a(ix-1)) 
      enddo
    endif
    time = time + dt    
  end subroutine update

  subroutine periodic_bc
    a(-1) = a(nx) 
    a(nx+1) = a(0)
  endsubroutine periodic_bc

end module solver

module diagnostics

  use shared_data

  implicit none

  contains

  subroutine do_io
    integer :: out_unit =10
    open(out_unit, file="xb.dat", access="stream")
    write(out_unit) xb
    close(out_unit)
    open(out_unit, file="a.dat", access="stream")
    write(out_unit) a
    close(out_unit)
    call execute_command_line("python plot_advect_1d_FD.py")
  end subroutine do_io

end module diagnostics

program advect_1d_FD

  use shared_data
  use setup
  use solver
  use diagnostics

  implicit none

  call init

  do 
    if ((step >= nsteps .and. nsteps >= 0) .or. (time >= t_end)) exit
    step = step + 1
    call update
    if (verbose) print *,'step',step, 'time',time,'dt',dt
  end do 

  call do_io

  print *, 'done in',step,'steps', 'with cfl',cfl

end program advect_1d_FD
