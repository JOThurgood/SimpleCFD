program test4

  use multigrid

  implicit none

  integer, parameter :: num=selected_real_kind(p=15)
  real(num) :: pi = 4.0_num * ATAN(1.0_num)
  real(num) :: pi2 = 8.0_num * ATAN(1.0_num)

  integer :: nx, ny, ix, iy

  real(num) :: x,y,  dx, dy, L2, x_min, x_max, y_min, y_max

  real(num), dimension(:), allocatable :: xc, yc
  real(num), dimension(:,:), allocatable :: f,analytic, phi, eta
  real(num), dimension(:), allocatable :: L2_arr, n_arr

  integer :: power, power_min, power_max

  type(mg_input) :: input

  ! setup a test problem 

  print *,'Test4: This does something'

  power_min = 3
  power_max = 10

  allocate(L2_arr(1:1+power_max-power_min))
  allocate(n_arr(1:1+power_max-power_min))
  L2_arr = 1e6_num
  n_arr = 0

  different_resolutions: do power= power_min, power_max

    nx = 2**power 
    ny = nx
    x_min = 0.0_num
    x_max = 1.0_num
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
      y = yc(iy)
      f(ix,iy) = -16.0_num * pi**2 * (cos(pi2*x)*cos(pi2*y) +1.0_num) * sin(pi2*x) * sin(pi2*y)
!      f(ix,iy) = -4.0_num * pi**2 * ( 4.0_num * sin(pi2 * x) * sin(pi2 * y) + sin(2.0_num * pi2 *x) * sin(2.0_num * pi2 *y))
    enddo
    enddo
  
    allocate(eta(1:nx,1:ny)) ! again, redundant ghosts but for comparison to CCAPS
    do iy = 1, ny
    do ix = 1, nx
     eta(ix,iy) = 2.0_num + cos(2.0_num * pi * xc(ix)) * cos(2.0_num * pi * yc(iy)) 
    enddo
    enddo 
    ! solve for phi
  


  input = mg_input(tol = 1e-12_num, nx=nx, ny = ny, dx=dx, dy=dy, f = f, phi = phi, &
            & bc_xmin = 'periodic', bc_ymin='periodic', bc_xmax='periodic', bc_ymax = 'periodic', &
            eta = eta, eta_present = .true., & 
            & eta_bc_xmin = 'periodic', eta_bc_ymin='periodic', eta_bc_xmax='periodic', eta_bc_ymax = 'periodic')

  call mg_interface(input)

  phi = input%phi

    ! test against analytical solution
     
    do iy = 1, ny
    do ix = 1, nx
      x = xc(ix)
      y = yc(iy)
      analytic(ix,iy) = sin(pi2*x)*sin(pi2*y)
    enddo
    enddo

    ! subtract the average of phi since there is nothing to anchor the solution
    phi = phi- sum(phi(1:nx,1:ny)) / real(nx*ny,num)  
    
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
    deallocate(eta)

  enddo different_resolutions

  print *, 'L2 _arr',L2_arr
  print *, 'n _arr',n_arr

  call execute_command_line("rm -rf test4_l2.dat")
  call execute_command_line("rm -rf test4_nx.dat")

  open(10, file="test4_l2.dat", access="stream")
  write(10) L2_arr
  close(10)

  open(10, file="test4_nx.dat", access="stream")
  write(10) n_arr
  close(10)

  call execute_command_line("python test4_plots.py")
  call execute_command_line("rm test4_l2.dat")
  call execute_command_line("rm test4_nx.dat")


end program test4
