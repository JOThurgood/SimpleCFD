!  finite-difference implementation of upwind / FTCS for linear advection.
! 
!  We are solving a_t + u a_x = 0
! 
!  The FTCS discretization is: anew = aold + (C/2) (aold_{i+1} - aold_{i-1})
! 
!  where C is the CFL number
! 
!  J. Thurgood 2018-06-27 
!
!  Alogorithm follows  M. Zingale's lecturenotes/textbook
!  on astro CFD. Independently written in Fortran. 
!
!  needs the pyplot stuff removed and own independent python plotting script 
!  added and called from within code, to make as portable / dependence free as
!  possible

MODULE shared_data

  implicit none
  integer, parameter  :: num=selected_real_kind(p=15)
  integer (num)  :: nx, ix
  integer (num) :: step = 0,nsteps
  real(num) :: u
  real (num) :: x_min, x_max, dxb
  real (num), dimension(:), allocatable :: a
  real(num), dimension(:), allocatable :: xb
  real (num) :: t_end, time = 0.0_num, dt, CFL 
  LOGICAL :: FTCS,upwind
  !anticipating staggered grid in late exercises
  !so 'c' for cell centred points, 'b' for boundary 

  CONTAINS

END MODULE shared_data

MODULE setup

  USE shared_data

  implicit none

  CONTAINS

  SUBROUTINE user_control
    nx = 64   
    x_min = 0.0_num 
    x_max = 1.0_num
    t_end = 1.00_num
    u = 1.0_num !char speed of linear advection equation
    CFL = 0.1_num !CFL number
    nsteps = -1
    FTCS = .FALSE.
    upwind = .TRUE.
  END SUBROUTINE user_control 

  SUBROUTINE initial_conditions
    !step function
    a = 0.0_num
    do ix = -1,nx+1 
      if ((xb(ix) >= 1.0_num/3.0_num) .and. (xb(ix) <= 2.0_num / 3.0_num))  then
        a(ix) = 1.0_num
      endif
    end do 
  ENDSUBROUTINE initial_conditions

  SUBROUTINE setup_1d_fd_grid
    !setup 1d finite difference grid for var a of size nx and has one ghost either side
    !Here 0 is leftmost, nx is right, -1 and nx+1 are ghosts
    !And allocate "a" etc

    ALLOCATE(xb (-1:nx+1) )
    ALLOCATE(a (-1:nx+1) )

    dxb = (x_max - x_min) / REAL(nx,num) 
    DO ix = -1,nx+1 
      xb(ix) = REAL(ix,num) * dxb
    END DO 
  ENDSUBROUTINE setup_1d_fd_grid 

  SUBROUTINE init
    call user_control
    call setup_1d_fd_grid
    call initial_conditions
  END SUBROUTINE init

END MODULE setup

MODULE solver
  USE shared_data
  implicit none

  CONTAINS

  SUBROUTINE update
    dt = CFL * dxb / u
    call periodic_bc

    IF (FTCS) THEN
!    do ix=0,nx  WILL NOT WORK  - need previous state so any direction will break it
!      a(ix) =  a(ix) - CFL * (a(ix+1) - a(ix-1)) / 2.0_num
!    enddo
      a = a - CFL * ( CSHIFT(a,1) - CSHIFT(a,-1)) / 2.0_num
    ENDIF 

    IF (upwind) THEN
    !upwind - must either store old solution separately or fill from right to left
      do ix=nx,0,-1
        a(ix) =  a(ix) - CFL * (a(ix) - a(ix-1)) 
      enddo
!      a = a - CFL * (a - CSHIFT(a,-1)) <- havent tested yet
    ENDIF
    time = time + dt    
  END SUBROUTINE update

  SUBROUTINE periodic_bc
    a(-1) = a(nx) 
    a(nx+1) = a(0)
  ENDSUBROUTINE periodic_bc

END MODULE solver

MODULE diagnostics
  USE pyplot_module, only : pyplot
  USE shared_data
  implicit none

  CONTAINS

  SUBROUTINE do_io
    integer :: out_unit =10
    open(out_unit, file="xb.dat", access="stream")
    write(out_unit) xb
    close(out_unit)
    open(out_unit, file="a.dat", access="stream")
    write(out_unit) a
    close(out_unit)
  END SUBROUTINE do_io

  SUBROUTINE do_a_plot
    use, intrinsic :: iso_fortran_env, only : wp => real64
    integer :: istat
    type(pyplot) :: plt
    call plt%initialize(grid=.true.,xlabel='x',figsize=[20,10],&
                        title='plot test',legend=.true.,axis_equal=.false.,&
                        tight_layout=.true.)
    call plt%add_plot(xb,a,label='a',linestyle='b-o',markersize=5,linewidth=2,istat=istat)
    call plt%savefig('output.png', pyfile='plottest.py',istat=istat)!
  END SUBROUTINE do_a_plot

END MODULE diagnostics

PROGRAM fdadvect

  USE shared_data
  USE setup
  USE solver
  USE diagnostics
  implicit none

  call init

  DO 
    IF ((step >= nsteps .AND. nsteps >= 0) .OR. (time >= t_end)) EXIT
    step = step + 1
    call update
!    print *,'step',step, 'time',time,'dt',dt
  END DO 

!  call do_io
  call do_a_plot

  print *, 'Done in',step,'steps', 'with CFL',CFL
END PROGRAM fdadvect
