program test2_constHelmholtz

  use multigrid

  implicit none

  integer, parameter :: num=selected_real_kind(p=15)
!  real(num) :: pi = 4.0_num * ATAN(1.0_num)

  integer :: nx, ny, ix, iy

  real(num) :: dx, dy, L2, x_min, x_max, y_min, y_max

  real(num), dimension(:), allocatable :: xc, yc
  real(num), dimension(:,:), allocatable :: f,analytic, phi
  real(num), dimension(:), allocatable :: L2_arr, n_arr

  integer :: power, power_min, power_max

  type(mg_input) :: input

  ! setup a test problem 

  
  print *,' As test 2 but passed to the const helmholtz mode with alpha=0 beta = -1 '
  print *,' (i.e., test the relaxation for CH mode reduces to poissions equation)'


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
      f(ix,iy) = -2.0_num * ((1.0_num-6.0_num*xc(ix)**2)*(yc(iy)**2)*(1.0_num-yc(iy)**2) &
        & + (1.0_num-6.0_num*yc(iy)**2)*(xc(ix)**2)*(1.0_num-xc(ix)**2))
    enddo
    enddo
  
  
!    ! Solve for phi
!    nlevels = -1 ! auto
!    call mg_interface(f = f, phi=phi, tol = 1e-12_num, &
!                    & nx = nx, ny = ny, dx = dx, dy = dy, &
!                    &  nlevels = nlevels, &
!    & bc_xmin_in = fixed, bc_ymin_in = fixed, bc_xmax_in = fixed, bc_ymax_in = fixed)


  input = mg_input(tol = 1e-12_num, nx=nx, ny = ny, dx=dx, dy=dy, f = f, phi = phi, &
            & bc_xmin = 'fixed', bc_ymin='fixed', bc_xmax='fixed', bc_ymax = 'fixed', &
            & const_helmholtz = .true., ch_alpha = 0.0_num, ch_beta = -1.0_num, &
            & deallocate_after = .true.)

  call mg_interface(input)

  phi = input%phi
  
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

  enddo different_resolutions

  print *, 'L2 _arr',L2_arr
  print *, 'n _arr',n_arr

  call execute_command_line("rm -rf test2_l2.dat")
  call execute_command_line("rm -rf test2_nx.dat")

  open(10, file="test2_l2.dat", access="stream")
  write(10) L2_arr
  close(10)

  open(10, file="test2_nx.dat", access="stream")
  write(10) n_arr
  close(10)

  call execute_command_line("python test2_plots.py")
  call execute_command_line("rm test2_l2.dat")
  call execute_command_line("rm test2_nx.dat")


end program test2_constHelmholtz
