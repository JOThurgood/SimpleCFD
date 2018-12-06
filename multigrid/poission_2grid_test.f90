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

    ! assume nx is divisible by 2
    if (mod(nx,2) /= 0) then
      print *,'idiot'
      STOP
    endif

    allocate(xc(0:nx+1)) ! 1 ghost cell either side needed
    allocate(f(1:nx))
    allocate(phi(0:nx+1))

    dx = (x_max - x_min) / real(nx,num)

    xc(0) = x_min - dx/2.0_num 
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

  subroutine solve_gs(myphi,rhs,tol) 

    real(num), dimension(:), allocatable, intent(inout) :: myphi
    real(num), dimension(:), allocatable, intent(in) :: rhs
    real(num), intent(in) :: tol
    real(num) :: L_myphi ! discritised laplacian of myphi
    real(num) :: L2_old, L2 ! norms 
    integer :: ir

    real(num) :: start, finish

    call cpu_time(start)

    L2_old = 1e6_num 

    ir = 0

    call bcs(myphi)

    do
      ir = ir + 1

      do ix = 1, nx
        myphi(ix) = 0.5_num * (myphi(ix-1) + myphi(ix+1) - dx**2 * rhs(ix)) 
      enddo

      call bcs(myphi)
      
      L2 = 0.0_num
      do ix = 1, nx
        L_myphi = ( myphi(ix-1) - 2.0_num * myphi(ix) + myphi(ix+1) ) / dx**2
        L2 = L2 + abs(rhs(ix) - L_myphi)**2
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

  subroutine bcs(myphi)

    real(num), dimension(:), allocatable, intent(inout) :: myphi

    ! Neumann (dphi/dx =0 at x=0) BC
    ! i.e. zero gradient
    myphi(0) = myphi(1)
    ! Dirichlet (phi = 0 at x = x_max) 
    myphi(nx+1) = -myphi(nx)

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

    real(num), dimension(1:nx) :: residual

    real(num) :: Lap ! Laplacian
    real(num) :: L2, L2_old ! norms

    integer :: nsteps ! number of steps through the algoritm (main do loop)
    integer :: nop ! number of operations through a 1D 
    integer :: c ! counter for sub iterations


    ! initialise a coarse grid and data

    ! the main cycle

    L2_old = 1e6_num ! init as huge
    nsteps = 0
    nop = 0 
    do ! main do
      nsteps = nsteps + 1

      ! 3 iterations of GS relaxation on initial (finest) grid
      do c = 1, 3
        call bcs(myphi) 
        do ix = 1, nx
          myphi(ix) = 0.5_num * (myphi(ix-1) + myphi(ix+1) - dx**2 * rhs(ix)) 
        enddo
      enddo 

      ! see if converged to tolerance on this finest grid and exit 
      ! if satisfied
  
      call bcs(myphi)

      L2 = 0.0_num
      do ix = 1, nx
        Lap = ( myphi(ix-1) - 2.0_num * myphi(ix) + myphi(ix+1) ) / dx**2
        residual(ix) = Lap - rhs(ix)
      enddo
      L2 = sqrt(sum(abs(residual)**2)/real(nx,num))

      if (abs(L2-L2_old) <= tol) exit
      L2_old = L2

      ! restrict to a coarse grid
      
    enddo !main do

    ! output

    print *, '*** solve_mg completed:'
    print *, '****** ncycles:',nsteps
!    print '(" ****** cpu_time: ",f6.3," seconds.")',finish-start
    print *, '****** tolerance (L2-L2_old <=tol for exit)',tol 
    print *, '****** L2 norm on residual', L2

  end subroutine solve_mg

  subroutine bcs(myphi)

    real(num), dimension(:), allocatable, intent(inout) :: myphi

    ! Neumann (dphi/dx =0 at x=0) BC
    ! i.e. zero gradient
    myphi(0) = myphi(1)
    ! Dirichlet (phi = 0 at x = x_max) 
    myphi(nx+1) = -myphi(nx)

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
  print *,'Gauss-Seidel Iteration'
  print *,'*******************'
  print *, '***Start GS relaxation'
  phi = 0.0_num ! initial guess, arbitary
  call solve_gs(myphi = phi, rhs = f, tol = 1e-16_num) 
  print *, '*** Comparing to analytical sln'
  call check_vs_analytic

  ! Two Level Multigrid
  print *,''
  print *, 'Two-level Multi Grid Relaxation'
  print *,'*******************'
  print *, '***Start MG2 relaxation'
  phi = 0.0_num ! reset to same initial guess
  call solve_mg(myphi = phi, rhs = f, tol = 1e-16_num)
  print *, '*** Comparing to analytical sln'
  call check_vs_analytic

end program poission_2dgrid_test
