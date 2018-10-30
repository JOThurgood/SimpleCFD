! Simple illustration of projection / test
! Prescribe a divergence free velocity field, 
! polute it with prescribed divergence, then
! use a relaxation method to reconstruct 
! the initial divergence free field.

! See e.g. 

program smooth1

  implicit none

  ! declarations

  integer, parameter :: num = selected_real_kind(p=15)
  real(num), parameter :: pi = 4.0_num * ATAN(1.0_num)
  integer :: out_unit = 10 

  integer :: nx, ny 
  integer :: ix,iy
  integer(num) :: step = 0
  integer(num) :: nsteps 

  real(num) :: tol
  real(num) :: x_min, x_max 
  real(num) :: y_min, y_max 
  real(num) :: dx, dy
  real(num) :: L_phi ! discritised laplacian 
  real(num) :: residual

  real(num), dimension(:), allocatable :: x, y

  real(num), dimension(:,:), allocatable :: u_orig !origional 
  real(num), dimension(:,:), allocatable :: u_pol  !polluted 
  real(num), dimension(:,:), allocatable :: u_reconstructed
  real(num), dimension(:,:), allocatable :: v_orig !origional 
  real(num), dimension(:,:), allocatable :: v_pol  !polluted 
  real(num), dimension(:,:), allocatable :: v_reconstructed
  real(num), dimension(:,:), allocatable :: phi
  real(num), dimension(:,:), allocatable:: D_pol ! numerical div of polluted velocity field
  real(num) :: L2, old_L2

  ! simulation / domain / grid control parameters

  nsteps = 10000 ! max no of steps before crashing out
  nx = 64 
  ny = nx ! currently hardcoded for dx = dy later in the code
  x_min = 0.0_num 
  x_max = 1.0_num
  y_min = x_min !
  y_max = x_max
  tol = 1e-14_num !less than 1e-2 takes too long -> clear limitation of the method and slow convergence?

  ! array allocation and initialisation

  ! -- Grid 

  allocate(x(0:nx+1))  ! cell center (cc) - counts from 1 to nx
  allocate(y(0:ny+1))  ! cell center (cc) - counts from 1 to ny
                        ! ONE GUARD CELL FOR THIS PROBLEM

  dx = (x_max - x_min) / REAL(nx-1,num)
  x(0) = x_min - dx
  do ix = 1, nx+1
    x(ix) = x(ix-1) + dx
  end do

  dy = (y_max - y_min) / REAL(ny-1,num)
  y(0) = y_min - dy
  do iy = 1, ny+1
    y(iy) = y(iy-1) + dy
  end do

  if (dy /= dx) then
    print *, 'Warning: dy /= dx set. Has not yet being generalised.'
    print *, 'Terminating'
    STOP
  end if 
 
  ! -- Velocity 

  allocate(u_orig(0:nx+1,0:ny+1))
  allocate(u_pol(0:nx+1,0:ny+1))
  allocate(u_reconstructed(0:nx+1,0:ny+1))
  allocate(v_orig(0:nx+1,0:ny+1))
  allocate(v_pol(0:nx+1,0:ny+1))
  allocate(v_reconstructed(0:nx+1,0:ny+1))
 

  do iy = 0, ny + 1
  do ix = 0, nx + 1

    u_orig(ix,iy) = - sin(1.0_num * pi * x(ix))**2 & 
      & * sin(2.0_num * pi * y(iy))
    v_orig(ix,iy) =   sin(1.0_num * pi * y(iy))**2 &
      & * sin(2.0_num * pi * x(ix))

   u_pol(ix,iy) = u_orig(ix,iy) - &
     & 0.2_num * pi * cos(2.0_num * pi * y(iy)) &
     & * sin(2.0_num * pi * x(ix)) 
   v_pol(ix,iy) = v_orig(ix,iy) - &
     & 0.2_num * pi * sin(2.0_num * pi * y(iy)) & 
     & * cos(2.0_num * pi * x(ix))
    
  end do
  end do 

  ! -- initial guess on phi 

  allocate(phi(0:nx+1,0:ny+1))


!  do iy = 0, ny + 1
!  do ix = 0, nx + 1
!    phi(ix,iy) = (u_pol(ix+1,iy) - u_pol(ix-1,iy)) / 2.0_num / dx + &
!      & (v_pol(ix,iy+1) - v_pol(ix,iy-1)) / 2.0_num / dy
!  end do
!  enddo 

  phi = 0.0_num  

  ! Solve Poisson equation via relaxation for phi

  print *, 'begining relaxation'


  allocate(D_pol(0:nx,0:ny))

  do iy = 1, ny 
  do ix = 1, nx 
    D_pol(ix,iy) = (u_pol(ix+1,iy) - u_pol(ix-1,iy)) / 2.0_num / dx + &
            & (v_pol(ix,iy+1) - v_pol(ix,iy-1)) / 2.0_num / dy
  end do
  end do


  do
    step = step + 1

    ! Periodic boundary conditions on phi

    phi(0,:) = phi(nx,:)
    phi(nx+1,:) = phi(1,:)
    phi(:,0) = phi(:,ny)
    phi(:,ny+1) = phi(:,1)

    do iy = 1, ny 
    do ix = 1, nx 

      phi(ix,iy) = 0.25_num * ( &
        & phi(ix+1,iy) + phi(ix-1,iy) + phi(ix,iy+1) + phi(ix,iy-1) &
        - dx**2 * D_pol(ix,iy) ) 

    end do
    end do 

    ! re-apply bc ? perversely seems to have a -ve effect on error

    !phi(0,:) = phi(nx,:)
    !phi(nx+1,:) = phi(1,:)
    !phi(:,0) = phi(:,ny)
    !phi(:,ny+1) = phi(:,1)

    ! reset the maxmimum val of residual
    residual = -1e6_num 
    L2 = 0.0_num
    do iy = 1, ny 
    do ix = 1, nx 
      L_phi = (phi(ix+1,iy) - 2.0_num*phi(ix,iy) + phi(ix-1,iy)) &
        & / dx**2 + &
        & (phi(ix,iy+1) - 2.0_num*phi(ix,iy) + phi(ix,iy-1)) / dy**2 

      residual = max(residual, abs(D_pol(ix,iy) - L_phi))
      L2 = L2 + abs(D_pol(ix,iy)-L_phi)**2
    end do
    end do 
      L2 = L2 / REAL(nx*ny,num)
    print *, 'Step',step,'completed, residual error is',residual
    print *, 'Step',step,'completed, L2 is',residual
 
   if ((step >= nsteps .and. nsteps >= 0) .or. (abs(L2-old_L2) <= tol)) exit
   old_L2 = L2
  end do

  print *, 'Relaxation completed in',step,'steps'


  ! reconstruct velocity field 

  print *, 'Reconstructing the velocity field' 

  phi(0,:) = phi(nx,:) !need to update BC here
  phi(nx+1,:) = phi(1,:)
  phi(:,0) = phi(:,ny)
  phi(:,ny+1) = phi(:,1)

  do iy = 1, ny
  do ix = 1, nx

    u_reconstructed(ix,iy) = u_pol(ix,iy) - &
      & (phi(ix+1,iy) - phi(ix-1,iy))/2.0_num/dx 
 
    v_reconstructed(ix,iy) = v_pol(ix,iy) - &
      & (phi(ix,iy+1) - phi(ix,iy-1))/2.0_num/dy



  enddo
  enddo

  u_reconstructed(0,:) = u_reconstructed(nx,:)
  u_reconstructed(nx+1,:) = u_reconstructed(1,:)
  v_reconstructed(:,0) = v_reconstructed(:,ny)
  v_reconstructed(:,ny+1) = v_reconstructed(:,1)


  print *, 'Reconstruction complete'
  print *, 'Maximum absolute difference between origional and', & 
    & ' reconstructed field is', &
    & max(maxval(abs(u_orig-u_reconstructed)), &
        & maxval(abs(v_orig-v_reconstructed)))

  ! calculate the l2 norm ? 

  ! generate output 
 
  call execute_command_line("rm -rf *.dat *.png") 
    !dont want to stream to existing files
  
  open(out_unit, file="x.dat", access="stream")
  write(out_unit) x(1:nx)
  close(out_unit)

  open(out_unit, file="y.dat", access="stream")
  write(out_unit) y(1:ny)
  close(out_unit)
  
  open(out_unit, file="u_orig.dat", access="stream")
  write(out_unit) u_orig(1:nx,1:ny)
  close(out_unit)
  
  open(out_unit, file="v_orig.dat", access="stream")
  write(out_unit) v_orig(1:nx,1:ny)
  close(out_unit)

  open(out_unit, file="u_pol.dat", access="stream")
  write(out_unit) u_pol(1:nx,1:ny)
  close(out_unit)
  
  open(out_unit, file="v_pol.dat", access="stream")
  write(out_unit) v_pol(1:nx,1:ny)
  close(out_unit)

  open(out_unit, file="u_rec.dat", access="stream")
  write(out_unit) u_reconstructed(1:nx,1:ny)
  close(out_unit)

  open(out_unit, file="v_rec.dat", access="stream")
  write(out_unit) v_reconstructed(1:nx,1:ny)
  close(out_unit)

  call execute_command_line("python plot_smooth1.py")

end program
