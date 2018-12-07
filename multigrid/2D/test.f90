program test


  use shared_data
  use multigrid_2d

  implicit none

  real(num), dimension(:,:), allocatable :: utrue, vtrue
  real(num), dimension(:,:), allocatable :: upol, vpol
  real(num), dimension(:,:), allocatable :: divu

  real(num) :: L2

  ! 0. User options

  nx = 64
  ny = nx
  x_min = 0.0_num
  x_max = 1.0_num
  y_min = x_min
  y_max = x_max

  call setup !setup the shared data

  ! setup a test problem 

  allocate(utrue(0:nx+1,0:ny+1)) !no need for ghosts in the test case
  allocate(vtrue(0:nx+1,0:ny+1))
  allocate(upol(0:nx+1,0:ny+1))
  allocate(vpol(0:nx+1,0:ny+1))
  allocate(divu(1:nx,1:ny))

  do iy = 0, ny+1
  do ix = 0, nx+1
    utrue(ix,iy) = - sin(1.0_num * pi * xc(ix))**2 & 
      & * sin(2.0_num * pi * yc(iy))
    vtrue(ix,iy) =   sin(1.0_num * pi * yc(iy))**2 &
      & * sin(2.0_num * pi * xc(ix))

    upol(ix,iy) = utrue(ix,iy) - & 
      & 0.2_num * pi * cos(2.0_num * pi * yc(iy)) &
      & * sin(2.0_num * pi * xc(ix)) 
    vpol(ix,iy) = vtrue(ix,iy) - & 
      & 0.2_num * pi * sin(2.0_num * pi * yc(iy)) & 
      & * cos(2.0_num * pi * xc(ix))
  enddo
  enddo


  do iy = 1, ny 
  do ix = 1, nx 
    divu(ix,iy) = (upol(ix+1,iy) - upol(ix-1,iy)) / 2.0_num / dx + & 
            & (vpol(ix,iy+1) - vpol(ix,iy-1)) / 2.0_num / dy
  enddo
  enddo

  print *, 'max divu before cleaning',maxval(divu)
  L2 = sqrt(sum(abs(utrue(1:nx,1:ny)-upol(1:nx,1:ny))**2) / real(nx*ny,num))
  print *, 'L2 (utrue vs upol) before cleaning',L2
  L2 = sqrt(sum(abs(vtrue(1:nx,1:ny)-vpol(1:nx,1:ny))**2) / real(nx*ny,num))
  print *, 'L2 (vtrue vs vpol) before cleaning',L2

  ! solve it

  call solve_mg_2d(mgphi = phi, f =  divu, tol = 1e-12_num)
  
  ! correct it 

  do iy = 1, ny
  do ix = 1, nx
    upol(ix,iy) = upol(ix,iy) - 0.5_num*(phi(ix+1,iy)-phi(ix-1,iy))/dx
    vpol(ix,iy) = vpol(ix,iy) - 0.5_num*(phi(ix,iy+1)-phi(ix,iy-1))/dy
  enddo
  enddo

  ! apply periodic bc to the polluted vels and then re calc divu
  upol(0,:) = upol(nx,:)
  upol(nx+1,:) = upol(1,:)
  upol(:,0) = upol(:,ny)
  upol(:,ny+1) = upol(:,1)
  vpol(0,:) = vpol(nx,:)
  vpol(nx+1,:) = vpol(1,:)
  vpol(:,0) = vpol(:,ny)
  vpol(:,ny+1) = vpol(:,1)

  do iy = 1, ny 
  do ix = 1, nx 
    divu(ix,iy) = (upol(ix+1,iy) - upol(ix-1,iy)) / 2.0_num / dx + & 
            & (vpol(ix,iy+1) - vpol(ix,iy-1)) / 2.0_num / dy
  enddo
  enddo

  print *, 'max divu after cleaning',maxval(divu)
  L2 = sqrt(sum(abs(utrue(1:nx,1:ny)-upol(1:nx,1:ny))**2) / real(nx*ny,num))
  print *, 'L2 (utrue vs upol) after cleaning',L2
  L2 = sqrt(sum(abs(vtrue(1:nx,1:ny)-vpol(1:nx,1:ny))**2) / real(nx*ny,num))
  print *, 'L2 (vtrue vs vpol) after cleaning',L2

end program test
