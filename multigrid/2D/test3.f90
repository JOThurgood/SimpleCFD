program test3

  use multigrid

  implicit none

  integer, parameter :: num=selected_real_kind(p=15)
!  real(num) :: pi = 4.0_num * ATAN(1.0_num)

  integer :: nx, ny, ix, iy
  integer :: nlevels

  real(num) :: dx, dy, L2, x_min, x_max, y_min, y_max

  real(num), dimension(:), allocatable :: xc, yc
  real(num), dimension(:,:), allocatable :: f,analytic, phi, eta
  real(num), dimension(:), allocatable :: L2_arr, n_arr

!  integer :: bc_xmin , bc_xmax, bc_ymin, bc_ymax
  integer, parameter :: fixed = 2
  
  integer :: power, power_min, power_max

  ! setup a test problem 

  
  print *,'Test3: This does something'

  power_min = 5
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
      f(ix,iy) = -2.0_num * ((1.0_num-6.0_num*xc(ix)**2)*(yc(iy)**2)*(1.0_num-yc(iy)**2) &
        & + (1.0_num-6.0_num*yc(iy)**2)*(xc(ix)**2)*(1.0_num-xc(ix)**2))
    enddo
    enddo
  
    allocate(eta(1:nx,1:ny)) ! again, redundant ghosts but for comparison to CCAPS
    eta = 1.0_num 
  
    ! solve for phi
    nlevels = -1!-1 ! auto
    call mg_interface(f = f, phi=phi, eta=eta, tol = 1e-12_num, &
                    & nx = nx, ny = ny, dx = dx, dy = dy, &
                    &  nlevels = nlevels, &
    & bc_xmin_in = fixed, bc_ymin_in = fixed, bc_xmax_in = fixed, bc_ymax_in = fixed)
  
    ! test against analytical solution
     
    do iy = 1, ny
    do ix = 1, nx
      analytic(ix,iy) = (xc(ix)**2-xc(ix)**4) * (yc(iy)**4-yc(iy)**2)
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
    deallocate(eta)

  enddo different_resolutions

  print *, 'L2 _arr',L2_arr
  print *, 'n _arr',n_arr

  call execute_command_line("rm -rf test3_l2.dat")
  call execute_command_line("rm -rf test3_nx.dat")

  open(10, file="test3_l2.dat", access="stream")
  write(10) L2_arr
  close(10)

  open(10, file="test3_nx.dat", access="stream")
  write(10) n_arr
  close(10)

  call execute_command_line("python test3_plots.py")
  call execute_command_line("rm test3_l2.dat")
  call execute_command_line("rm test3_nx.dat")


end program test3
