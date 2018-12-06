module shared_data

  implicit none

  integer :: out_unit = 10
  integer, parameter :: num = selected_real_kind(p=15) 
  real(num), parameter :: pi2 = 8.0_num * atan(1.0_num)

  integer :: nx, ix
  real(num) :: x_min,x_max, dx 

  real(num) :: A, k

  real(num), dimension(:), allocatable :: xc ! x at cell centers
  real(num), dimension(:), allocatable :: phi ! phi at cell centers
  real(num), dimension(:), allocatable :: f ! f at cell centers

  contains 

  subroutine initial_setup ! allocate and initialise

    allocate(xc(0:nx+1)) ! 1 ghost cell either side needed
    allocate(f(1:nx))
    allocate(phi(0:nx+1))

    dx = (x_max - x_min) / real(nx,num)

    xc(0) = x_min - dx/2. 
    do ix = 1,nx+1
      xc(ix) = xc(ix-1) + dx
    enddo

    do ix = 1, nx ! choose f from RHS of  d^2/dx^2(phi) = f
      f(ix) = A * sin(k * pi2 * xc(ix) / (x_max-x_min))
    enddo

  end subroutine initial_setup

end module shared_data

module diagnostics

  use shared_data

  implicit none

  public :: check_vs_analytic

  contains

  subroutine check_vs_analytic

    real(num) :: analytic, diff, C1, C2, d
    real(num) :: L2 

    d = x_max-x_min

    C1 = A * d / k / pi2
    C2 = - C1 * d


    L2 = 0.0_num
    do ix = 1, nx
      analytic = - A * sin(k * pi2 * xc(ix) / d) * (d / k / pi2)**2 &
                & + C1 * xc(ix) + C2
      diff = phi(ix) - analytic
!print *, phi(ix),analytic, abs(diff)
      L2 = L2 + abs(diff)**2
    enddo

    L2 = sqrt(L2 / real(nx,num))

    print *, '****** L2 on computed vs analytic',L2

  end subroutine check_vs_analytic

end module diagnostics

module gauss_seidel

  use shared_data

  implicit none
  private
  public :: solve_gs

  contains

  subroutine solve_gs(phigs,rhs,tol) 

    real(num), dimension(:), allocatable, intent(inout) :: phigs
    real(num), dimension(:), allocatable, intent(in) :: rhs
    real(num), intent(in) :: tol
    real(num) :: L_phigs ! discritised laplacian of phigs
    real(num) :: L2_old, L2 ! norms 
    integer :: ir

    real(num) :: start, finish

    call cpu_time(start)

    L2_old = 1e6_num 

    ir = 0

    call bcs(phigs)

    do
      ir = ir + 1

      do ix = 1, nx
        phigs(ix) = 0.5_num * (phigs(ix-1) + phigs(ix+1) - dx**2 * rhs(ix)) 
      enddo

      call bcs(phigs)
      
      L2 = 0.0_num
      do ix = 1, nx
        L_phigs = ( phigs(ix-1) - 2.0_num * phigs(ix) + phigs(ix+1) ) / dx**2
        L2 = L2 + abs(rhs(ix) - L_phigs)**2
      enddo
      L2 = sqrt(L2 / real(nx,num))

      ! verbose
!      print *,'step',ir,'L2',L2

      ! exit criterion
      if (abs(L2-L2_old) <= tol) exit
!      if (ir >= 50) exit
      L2_old = L2
    enddo

    call cpu_time(finish)

    print *, '*** solve_gs completed:'
    print *, '****** steps:',ir
    print '(" ****** cpu_time: ",f6.3," seconds.")',finish-start
    print *, '****** tolerance (L2-L2_old <=tol for exit)',tol 
    print *, '****** L2 norm on residual', L2

  end subroutine solve_gs

  subroutine bcs(phigs)

    real(num), dimension(:), allocatable, intent(inout) :: phigs

    ! Neumann (dphi/dx =0 at x=0) BC
    ! i.e. zero gradient
    phigs(0) = phigs(1)
    ! Dirichlet (phi = 0 at x = x_max) 
    phigs(nx+1) = -phigs(nx)

  end subroutine bcs

end module gauss_seidel

module mg_2level

  use shared_data

  implicit none
  private
  public :: solve_mg

  contains

  subroutine solve_mg(myphi,rhs,tol) 

    real(num), dimension(:), allocatable, intent(inout) :: myphi
    real(num), dimension(:), allocatable, intent(in) :: rhs
    real(num), intent(in) :: tol

  end subroutine solve_mg

  subroutine bcs
  end subroutine bcs

end module mg_2level

program poission_2dgrid_test

  use shared_data
  use diagnostics
  use gauss_seidel
  use mg_2level

  implicit none

  ! user parameters here

  x_min = 0.0_num
  x_max = 1.0_num 
  nx = 32

  A = 10.0_num
  k = 4.0_num

  call initial_setup

  ! GS Relaxation (for comparisons sake)
  print *, 'Begining GS Relaxation'
  phi = 0.0_num ! initial guess, arbitary
  call solve_gs(phigs = phi, rhs = f, tol = 1e-16_num) 
  call check_vs_analytic
  print *, '*** GS Relaxation complete'

  ! Two Level Multigrid
  print *, 'Begining two-level Multi Grid Relaxation'
  phi = 0.0_num ! reset to same initial guess


end program poission_2dgrid_test
