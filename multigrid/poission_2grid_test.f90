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
    print '(" ****** cpu_time: ",f20.3," seconds.")',finish-start
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

  subroutine solve_mg(myphi,rhs,tol,ncoarse,nfine) 

    real(num), dimension(:), allocatable, intent(inout) :: myphi
    real(num), dimension(:), allocatable, intent(in) :: rhs
    real(num), intent(in) :: tol
    integer, intent(in) :: ncoarse, nfine

    real(num), dimension(1:nx) :: residual
    real(num), dimension(1:nx/2) :: residual_c!oarse
    real(num), dimension(:), allocatable :: epsi

    real(num) :: Lap ! Laplacian
    real(num) :: L2, L2_old ! norms

    integer :: nsteps ! number of steps through the algoritm (main do loop)
    integer :: nop ! number of iterative operations through a 1D array
    integer :: c ! counter for sub iterations

    integer :: ixc !coarse x index

    real(num) :: start, finish

    call cpu_time(start)

    ! initialise a coarse grid and data
    if (.not.(allocated(epsi))) allocate(epsi(0:nx/2+1))
    epsi = 0.0_num
    residual = 0.0_num

    ! the main cycle

    L2_old = 1e6_num ! init as huge
    nsteps = 0
    nop = 0 
    do ! main do
      nsteps = nsteps + 1

      ! 3 iterations of GS relaxation on initial (finest) grid
      do c = 1, ncoarse
        nop = nop + 1
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

      ! restrict residue to a coarse grid
      ix = 1
      do ixc = 1, nx/2
        residual_c(ixc) = 0.5_num * (residual(ix) + residual(ix+1))
        ix = ix + 2
      enddo

      ! GS relaxation on the coarse grid to estimate a solution to
      ! Lap * epsi = residual_c 

      do c = 1, nfine
        if (mod(c,2)==0) nop = nop + 1 !only do every 2 steps since coarse array is half the size of fine
        call coarse_bcs(epsi)
        do ixc = 1, nx /2
          epsi(ixc) = 0.5_num * (epsi(ixc-1) + epsi(ixc+1) - (dx*2.0_num)**2 * residual_c(ixc))
        enddo
      enddo

      ! prolongation of epsi to update the estimate of phi on the fine grid
      ix = 1
      do ixc = 1, nx/2
        ! direct injection - replace with 2nd order method down line?
        ! phi => phi - epsi
        myphi(ix) = myphi(ix) - epsi(ixc)
        myphi(ix+1) = myphi(ix+1) - epsi(ixc)
        ix = ix + 2
      enddo
     
    enddo !main do

    call cpu_time(finish)

    ! output

    print *, '*** solve_mg completed:'
    print *, '****** ncycles:',nsteps
    print *, '****** relative number of iterative steps through array(1:nx)',nop
    print '(" ****** cpu_time: ",f20.3," seconds.")',finish-start
    print *, '****** tolerance (L2-L2_old <=tol for exit)',tol 
    print *, '****** ncoarse iterations per cycle',ncoarse
    print *, '****** nfine iterations per cycle',nfine
    print *, '****** L2 norm on residual', L2

    ! deallocate the coarse grid

  end subroutine solve_mg

  subroutine bcs(myphi)

    real(num), dimension(:), allocatable, intent(inout) :: myphi

    ! Neumann (dphi/dx =0 at x=0) BC
    ! i.e. zero gradient
    myphi(0) = myphi(1)
    ! Dirichlet (phi = 0 at x = x_max) 
    myphi(nx+1) = -myphi(nx)

  end subroutine bcs

  subroutine coarse_bcs(epsi)

    real(num), dimension(:), allocatable, intent(inout) :: epsi

    ! since solving a residual eqn, always homogeneous
    ! for this test case, the fine BCs (i.e. bcs on phi) are also homogeneous
    ! so these are the same

    ! exept upper limit nx => nx/2 . In other words this subroutine
    ! exists purely for this reason -> should do a more compact solution

    ! Neumann (dphi/dx =0 at x=0) BC
    ! i.e. zero gradient
    epsi(0) = epsi(1)
    ! Dirichlet (phi = 0 at x = x_max) 
    epsi(nx/2+1) = -epsi(nx/2)

  end subroutine coarse_bcs

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
  nx = 128

  A = 10.0_num
  k = 4.0_num

  call initial_setup

  ! Some output to remind you of your parameters
  print *,''
  print *,'Multilevel relaxation tests on Poission eqn'
  print *,'*******************'
  print *,'nx =',nx
  

  ! Two Level Multigrid
  print *,''
  print *, 'Two-level Multi Grid Relaxation'
  print *,'*******************'
  print *, '***Start MG2 relaxation'
  phi = 0.0_num
  call solve_mg(myphi = phi, rhs = f, tol = 1e-16_num, ncoarse = 3, nfine =50)
  print *, '*** Comparing to analytical sln'
  call check_vs_analytic

  ! GS Relaxation (for comparisons sake)
  print *,''
  print *,'Gauss-Seidel Iteration'
  print *,'*******************'
  print *, '***Start GS relaxation'
  phi = 0.0_num 
  call solve_gs(myphi = phi, rhs = f, tol = 1e-16_num) 
  print *, '*** Comparing to analytical sln'
  call check_vs_analytic

end program poission_2dgrid_test
