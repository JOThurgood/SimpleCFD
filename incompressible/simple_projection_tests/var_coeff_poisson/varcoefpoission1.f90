program varcoeffpoisson1

  implicit none

  ! declarations

  integer, parameter :: num = selected_real_kind(p=15)
  real(num), parameter :: pi = 4.0_num * ATAN(1.0_num)
  real(num), parameter :: pi2 = 8.0_num * ATAN(1.0_num)
  integer :: out_unit = 10 

  logical :: verbose=.false.

  logical :: gsrb=.true.

  integer :: nx, ny 
  integer :: ix,iy
  integer(num) :: step = 0
  integer(num) :: nsteps 

  real(num) :: tol
  real(num) :: x_min, x_max 
  real(num) :: y_min, y_max 
  real(num) :: dx, dy
  real(num) :: L_phi ! discritised laplacian 

  real(num), dimension(:), allocatable :: x, y

  ! all parameterss at cc

  real(num), dimension(:,:), allocatable :: phi, phi_true
  real(num) :: phi_avg
  real(num), dimension(:,:), allocatable :: f !test function on RHS
  real(num), dimension(:,:), allocatable :: eta !variable coefficient
  real(num) :: eta_ip, eta_im, eta_jp, eta_jm ! coefficient eval at i+1/2 (ip), j-1/2 (jm) etc

  real(num) :: L2, old_L2

  ! simulation / domain / grid control parameters

  nsteps = int(1e6) ! max no of steps before crashing out
  nx = 128
  ny = nx ! currently hardcoded for dx = dy later in the code
  x_min = 0.0_num 
  x_max = 1.0_num
  y_min = x_min !
  y_max = x_max
  tol = 1e-16_num !1e-16 appropriate if using L2-L2_old type condition


  ! init out
  print *,'nx=',nx,'ny=',ny
  print *,'exit tolerance on L2(residual) - L2(residual)_old =',tol

  ! array allocation and initialisation

  ! -- Grid 

  allocate(x(0:nx+1))  ! cell center (cc) - counts from 1 to nx
  allocate(y(0:ny+1))  ! cell center (cc) - counts from 1 to ny
                        ! ONE GUARD CELL FOR THIS PROBLEM

  dx = (x_max - x_min) / REAL(nx,num)
  x(0) = x_min - dx/2.0_num
  do ix = 1, nx+1
    x(ix) = x(ix-1) + dx
  end do

  dy = (y_max - y_min) / REAL(ny,num)
  y(0) = y_min - dy/2.0_num
  do iy = 1, ny+1
    y(iy) = y(iy-1) + dy
  end do

  if (dy /= dx) then
    print *, 'Warning: dy /= dx set. Has not yet being generalised.'
    print *, 'Terminating'
    STOP
  end if 
 
  ! -- test f (RHS) and eta (variable cooefficient), and true solution of phi

  allocate(f(0:nx+1,0:ny+1))
  allocate(eta(0:nx+1,0:ny+1))
  allocate(phi_true(0:nx+1,0:ny+1))
 

  do iy = 0, ny + 1
  do ix = 0, nx + 1
    f(ix,iy) = -16.0_num * pi**2 * (cos(pi2*x(ix))*cos(pi2*y(iy))+1.0_num) * &
              & sin(pi2*x(ix))*sin(pi2*y(iy))
    eta(ix,iy) = 2.0_num + cos(pi2*x(ix))*cos(pi2*y(iy))
    phi_true(ix,iy) = sin(pi2*x(ix))*sin(pi2*y(iy))
  end do
  end do 

  ! -- initial guess on phi 

  allocate(phi(0:nx+1,0:ny+1))
  phi = 0.0_num 

  ! Solve variable coefficient Poisson equation via relaxation for phi

  print *, 'begining relaxation to solve for phi.'
  print *, 'this can take a while, use VERBOSE if you want to monitor stepping'

  do
    step = step + 1


    ! if not using redblack order
    if (.not. gsrb) then 
      do iy = 1, ny 
      do ix = 1, nx
        !eta_ip actually also has the deltas built in 
        eta_ip = 0.5_num * (eta(ix,iy)+eta(ix+1,iy)) / dx**2
        eta_im = 0.5_num * (eta(ix,iy)+eta(ix-1,iy)) / dx**2 
        eta_jp = 0.5_num * (eta(ix,iy)+eta(ix,iy+1)) / dy**2
        eta_jm = 0.5_num * (eta(ix,iy)+eta(ix,iy-1)) / dy**2

        phi(ix,iy) = eta_ip * phi(ix+1,iy) + eta_im*phi(ix-1,iy) &
          & + eta_jp * phi(ix,iy+1) + eta_jm * phi(ix,iy-1) - f(ix,iy)
        phi(ix,iy) = phi(ix,iy) / ( eta_ip + eta_im + eta_jp + eta_jm )
      end do
      end do 
   
    else !use red black

      ! odd iteration
      do iy = 1, ny 
      do ix = 1, nx 
        if (modulo(ix+iy,2) == 1) then
          eta_ip = 0.5_num * (eta(ix,iy)+eta(ix+1,iy)) / dx**2
          eta_im = 0.5_num * (eta(ix,iy)+eta(ix-1,iy)) / dx**2 
          eta_jp = 0.5_num * (eta(ix,iy)+eta(ix,iy+1)) / dy**2
          eta_jm = 0.5_num * (eta(ix,iy)+eta(ix,iy-1)) / dy**2

          phi(ix,iy) = eta_ip * phi(ix+1,iy) + eta_im*phi(ix-1,iy) &
            & + eta_jp * phi(ix,iy+1) + eta_jm * phi(ix,iy-1) - f(ix,iy)
          phi(ix,iy) = phi(ix,iy) / ( eta_ip + eta_im + eta_jp + eta_jm )
        endif
      end do
      end do 

      ! even iteration

      do iy = 1, ny 
      do ix = 1, nx 
        if (modulo(ix+iy,2) == 0) then
          eta_ip = 0.5_num * (eta(ix,iy)+eta(ix+1,iy)) / dx**2
          eta_im = 0.5_num * (eta(ix,iy)+eta(ix-1,iy)) / dx**2 
          eta_jp = 0.5_num * (eta(ix,iy)+eta(ix,iy+1)) / dy**2
          eta_jm = 0.5_num * (eta(ix,iy)+eta(ix,iy-1)) / dy**2

          phi(ix,iy) = eta_ip * phi(ix+1,iy) + eta_im*phi(ix-1,iy) &
            & + eta_jp * phi(ix,iy+1) + eta_jm * phi(ix,iy-1) - f(ix,iy)
          phi(ix,iy) = phi(ix,iy) / ( eta_ip + eta_im + eta_jp + eta_jm )
        endif
      end do
      end do 

    endif




    ! Apply periodic boundary conditions on phi's ghost cells

    phi(0,:) = phi(nx,:)
    phi(nx+1,:) = phi(1,:)
    phi(:,0) = phi(:,ny)
    phi(:,ny+1) = phi(:,1)

    !residual = -1e6_num ! reset the maximum value of the residual to <0 
    L2 = 0.0_num
    do iy = 1, ny 
    do ix = 1, nx 
      eta_ip = 0.5_num * (eta(ix,iy)+eta(ix+1,iy))
      eta_im = 0.5_num * (eta(ix,iy)+eta(ix-1,iy))
      eta_jp = 0.5_num * (eta(ix,iy)+eta(ix,iy+1))
      eta_jm = 0.5_num * (eta(ix,iy)+eta(ix,iy-1))

      L_phi = ( eta_ip * (phi(ix+1,iy)-phi(ix,iy)) &
            &    -  eta_im * (phi(ix,iy)-phi(ix-1,iy)) ) / dx**2 + &
            & ( eta_jp * (phi(ix,iy+1)-phi(ix,iy)) &
            &    -  eta_jm * (phi(ix,iy)-phi(ix,iy-1)) ) / dy**2  

      !residual = max(residual, abs(D_pol(ix,iy) - L_phi))
      L2 = L2 + abs(f(ix,iy)-L_phi)**2
    end do
    end do 
      L2 = sqrt( L2 / REAL(nx*ny,num))


    if (verbose) &
      & print *, 'Step',step,'complete. L2 is',L2,'|L2-old_L2| is',abs(L2-old_L2),'tol is',tol

    if ((step >= nsteps .and. nsteps >= 0) .or. (abs(L2-old_L2) <= tol)) exit
    old_L2 = L2

  end do

  print *, 'Relaxation completed in',step,'steps'
  print *, 'residual (numerical) L2 norm on f-Lphi',L2



  ! calculate the l2 norm of difference between phi and analytical
  ! field

  phi_avg = sum(phi(1:nx,1:ny)) / REAL(nx*ny,num)
  L2 = 0.0_num
  
  do ix = 1, nx
  do iy = 1, ny
    L2=  L2 + abs(phi(ix,iy)-phi_avg-phi_true(ix,iy) )**2
  enddo
  enddo
  L2 = sqrt(L2 / REAL(2*nx*ny,num)) 

  print *,'L2 norm of analytical differences',L2

  ! calculate L_phi for output /

  call execute_command_line("rm -rf *.dat *.png") 
    !dont want to stream to existing files
  
  open(out_unit, file="x.dat", access="stream")
  write(out_unit) x(1:nx)
  close(out_unit)

  open(out_unit, file="y.dat", access="stream")
  write(out_unit) y(1:ny)
  close(out_unit)
  
  open(out_unit, file="f.dat", access="stream")
  write(out_unit) f(1:nx,1:ny)
  close(out_unit)

  open(out_unit, file="eta.dat", access="stream")
  write(out_unit) eta(1:nx,1:ny)
  close(out_unit)

  open(out_unit, file="phi.dat", access="stream")
  write(out_unit) phi(1:nx,1:ny)-phi_avg
  close(out_unit)

  call execute_command_line("python plot.py")
  call execute_command_line("rm *.dat")
end program
