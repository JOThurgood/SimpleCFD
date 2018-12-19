program test6

  use multigrid

  implicit none

  integer, parameter :: num=selected_real_kind(p=15)
  real(num) :: pi = 4.0_num * ATAN(1.0_num)
  real(num) :: pi2 = 8.0_num * ATAN(1.0_num)

  integer :: nx, ny, ix, iy

  real(num) :: dx, dy, L2, x_min, x_max, y_min, y_max

  real(num), dimension(:), allocatable :: xc, yc
  real(num), dimension(:,:), allocatable :: f,analytic, phi
  real(num), dimension(:), allocatable :: L2_arr, n_arr

  real(num) :: x, y, A, B, k, C1, C2 

  integer :: power, power_min, power_max

  type(mg_input) :: input

  ! setup a test problem 



  A = 1.0_num
  B = 1.0_num
  k = 4.0_num


  
  print *,'Test6: Lphi= A * sin(k * 2pi *x)'
  print *,'x_min = zerograd'
  print *,'x_max = fixed = B = ',B
  print *,'2D => 1D via ybcs = periodic'
  print *,'also can test variable nx for convergence properties'

  power_min = 5
  power_max = 9 

  allocate(L2_arr(1:1+power_max-power_min))
  allocate(n_arr(1:1+power_max-power_min))
  L2_arr = 1e6_num
  n_arr = 0

  different_resolutions: do power= power_min, power_max

    nx = 2**power 
    ny = nx
    x_min = 0.0_num
    x_max = 1.0_num ! dont change these, I've hardcoded len = 1 
    y_min = x_min
    y_max = x_max
  
    allocate(xc(-1:nx+2))
    allocate(yc(-1:ny+2))
  
    dx = (x_max - x_min) / real(nx,num)
    xc(-1) = x_min - 3.0_num * dx / 2.0_num 
    do ix = 0,nx+2
      xc(ix) = xc(ix-1) + dx
    enddo
  
    dy = (y_max - y_min) / real(ny,num)
    yc(-1) = y_min - 3.0_num * dy / 2.0_num 
    do iy = 0,ny+2
      yc(iy) = yc(iy-1) + dy
    enddo
  
    allocate(phi(-1:nx+2,-1:ny+2)) ! again, redundant ghosts but for comparison to CCAPS
    phi = 0.0_num 
  
    allocate(analytic(1:nx,1:ny)) ! again, redundant ghosts but for comparison to CCAPS
    analytic = 0.0_num 
  
  
    allocate(f(1:nx,1:ny))
   
    do iy = 1, ny
    do ix = 1, nx
      x = xc(ix)
      y = yc(ix) 
      f(ix,iy) = A * sin(k * pi2 * x) 
    enddo
    enddo
  
  
 ! Solve for phi

  input = mg_input(tol = 1e-12_num, nx=nx, ny = ny, dx=dx, dy=dy, f = f, phi = phi, &
            & bc_xmin = 'zero_gradient', bc_xmax = 'fixed', phi_bc_xmax = 1.0_num, &
            & bc_ymin = 'periodic', bc_ymax = 'periodic', &
            & deallocate_after = .true.)

  call mg_interface(input)

  phi = input%phi
  
    ! test against analytical solution
     
    do iy = 1, ny
    do ix = 1, nx
      x = xc(ix)
      y = yc(ix)
      C1 = A / k / pi2
      C2 = B - C1 
      analytic(ix,iy) = - A * sin(k * pi2 * x) *(1.0_num/k/pi2)**2 + C1*x + C2
    enddo
    enddo
  
    L2 = sqrt(sum(abs(analytic(1:nx,1:ny)-phi(1:nx,1:ny))**2) / real(nx*ny,num))
    
    print *,'L2',L2
 
    L2_arr(power-power_min+1) = L2
    n_arr(power-power_min+1) = real(nx,num)

 
    ! deallocate all so can do again
  
    deallocate(xc)
    deallocate(yc)
    deallocate(f)
    deallocate(phi)
    deallocate(analytic)

  enddo different_resolutions

  print *, 'L2 _arr',L2_arr
  print *, 'n _arr',n_arr

  call execute_command_line("rm -rf test6_l2.dat")
  call execute_command_line("rm -rf test6_nx.dat")

  open(10, file="test6_l2.dat", access="stream")
  write(10) L2_arr
  close(10)

  open(10, file="test6_nx.dat", access="stream")
  write(10) n_arr
  close(10)

  call execute_command_line("python test6_plots.py")
  call execute_command_line("rm test6_l2.dat")
  call execute_command_line("rm test6_nx.dat")


end program test6
