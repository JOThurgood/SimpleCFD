! test sln of 2d diffusion equation
! under Crank-Nicholson discritisation
! and Gauss Seidel relaxation 

program implicit1

  implicit none

  integer, parameter :: num = selected_real_kind(p=15)
  integer :: out_unit = 10

  integer :: nx = 128
  integer :: ny = 128

  real(num) :: x_min = 0.0_num
  real(num) :: x_max = 1.0_num
  real(num) :: y_min = 0.0_num
  real(num) :: y_max = 1.0_num

  real(num) :: time = 0.0_num
  real(num) :: t_end = 0.02_num

  real(num) :: k = 1.0_num
  real(num) :: phi2 = 2.0_num !max val
  real(num) :: phi1 = 1.0_num !baseline

  real(num) :: CFL = 1.0_num

  real(num) :: tol = 1e-12_num

  integer :: tstep = 0
  integer :: rstep = 0
  integer :: ix, iy
  real(num) :: dx, dy
  real(num) :: r,xx,yy

  real(num) :: dt
  real(num), dimension(:), allocatable :: x, y
  real(num), dimension(:,:), allocatable :: phi, phi_an, f

  real(num) :: L_phi
  real(num) :: alpha, beta

  real(num) :: L2
   
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
   

  allocate(phi(0:nx+1,0:ny+1))
  allocate(phi_an(0:nx+1,0:ny+1))

  do ix = 0, nx + 1
  do iy = 0, ny + 1
    xx = x(ix) - 0.5_num
    yy = y(iy) - 0.5_num
    r = sqrt(xx**2 + yy**2)
    phi(ix,iy) = (phi2-phi1) * exp(-0.25_num * r**2 / k / 0.001_num) + phi1
! 1.0_num + amp * exp(-r**2/ (2.0_num *(fwtm/4.29193_num)**2))
  enddo
  enddo

!  GOTO 1 ! uncomment to just dump initial condiiton and check it, bypassing the loop

  allocate(f(1:nx,1:ny))

  ! evolution loop

  time = 0.0_num

  do
    if (time > t_end) exit

    dt = CFL * dx**2 / k

    alpha = 1.0_num
    beta = dt * k / 2.0_num

    do ix = 1, nx
    do iy = 1, ny

      ! function based on lagrangian of phi at the previous time level
      ! - dont update during relaxation

      L_phi = (phi(ix+1,iy) - 2.0_num*phi(ix,iy) + phi(ix-1,iy)) / dx**2 + &
        & (phi(ix,iy+1) - 2.0_num*phi(ix,iy) + phi(ix,iy-1)) / dy**2 
       
      f(ix,iy) = phi(ix,iy) + 0.5_num * dt * k * L_phi

    enddo
    enddo

    !forget initial phi ?
    !phi = 0.0_num

    !GS relaxation loop

    do
      tstep = tstep + 1
      rstep = rstep + 1      
       
      ! odd iteration
      do iy = 1, ny 
      do ix = 1, nx 
        if (modulo(ix+iy,2) == 1) then

          ! set to the adjacent phi values first, then 
          phi(ix,iy) = phi(ix+1,iy) + phi(ix-1,iy) + phi(ix,iy+1) + phi(ix,iy-1)
          phi(ix,iy) = (f(ix,iy) + beta * phi(ix,iy) / dx**2) / &
            & (alpha + 4.0_num * beta / dx**2)

        endif
      end do
      end do 

      ! even iteration

      do iy = 1, ny 
      do ix = 1, nx 
        if (modulo(ix+iy,2) == 0) then
          phi(ix,iy) = phi(ix+1,iy) + phi(ix-1,iy) + phi(ix,iy+1) + phi(ix,iy-1)
          phi(ix,iy) = (f(ix,iy) + beta * phi(ix,iy) / dx**2) / &
            & (alpha + 4.0_num * beta / dx**2)
        endif
      end do
      end do 



      ! Apply boundary conditions on phi's ghost cells here
      ! homogenous Neumann BC

      phi(0,:) = phi(1,:)
      phi(nx+1,:) = phi(nx,:)
      phi(:,0) = phi(:,1)
      phi(:,ny+1) = phi(:,ny)


      ! Apply periodic boundary conditions on phi's ghost cells

      !phi(0,:) = phi(nx,:)
      !phi(nx+1,:) = phi(1,:)
      !phi(:,0) = phi(:,ny)
      !phi(:,ny+1) = phi(:,1)
 

 
      L2 = 0.0_num
      do iy = 1, ny 
      do ix = 1, nx 

        !Lagrangian of phi at new time level / this relaxaed state

        L_phi = (phi(ix+1,iy) - 2.0_num*phi(ix,iy) + phi(ix-1,iy)) / dx**2 + &
          & (phi(ix,iy+1) - 2.0_num*phi(ix,iy) + phi(ix,iy-1)) / dy**2 

     
        !residual
        L2 = L2 + abs(f(ix,iy)-(alpha*phi(ix,iy)-beta*L_phi))**2
      end do
      end do 
      L2 = sqrt( L2 / REAL(nx*ny,num))
 
!      print *,'time',time,'rstep',rstep,'L2',L2
   
      if (L2 <= tol) exit
    enddo

    print *,'time', time,'tstep',tstep
    print *,'took', rstep,'iterations to converge to (residual) L2=',L2


    time = time + dt
  enddo



  !the analytical solution
  do ix = 1,nx
  do iy = 1,ny
    xx = x(ix) - 0.5_num
    yy = y(iy) - 0.5_num
    r = sqrt(xx**2 + yy**2)
    phi_an(ix,iy) = (phi2-phi1) * (0.001_num / (time+0.001_num)) &
      & * exp(-0.25_num * r**2 / k / (0.001_num + time)) + phi1
  enddo
  enddo

  !calculate the L2 norm against the analytical solution

  L2 = 0.0_num

  do ix = 1, nx
  do iy = 1, ny
    L2 = L2 + abs(phi_an(ix,iy)-phi(ix,iy))**2 
  enddo
  enddo
  L2 = sqrt( L2 / REAL(nx*ny,num))

  print *,'L2 error norm agains analtical solution is',L2

1 CONTINUE

  call execute_command_line("rm -rf *.dat *.png")
    !dont want to stream to existing files
 
  open(out_unit, file="x.dat", access="stream")
  write(out_unit) x(1:nx)
  close(out_unit)
 
  open(out_unit, file="y.dat", access="stream")
  write(out_unit) y(1:ny)
  close(out_unit)
 
  open(out_unit, file="phi.dat", access="stream")
  write(out_unit) phi(1:nx,1:ny)
  close(out_unit)

  open(out_unit, file="phi_an.dat", access="stream")
  write(out_unit) phi_an(1:nx,1:ny)
  close(out_unit)

  !calculate the difference of phi from the analytical solution
  ! (just call the difference phi to avoid allocating a new array)
!  do ix = 1, nx
!  do iy = 1, ny
!    phi(ix,iy) = phi(ix,iy) - 
!  enddo
!  enddo

  call execute_command_line("python plot.py")
  call execute_command_line("rm -rf *.dat")




end program implicit1


