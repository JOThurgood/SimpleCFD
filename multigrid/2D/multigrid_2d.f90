module multigrid_2d

  use shared_data

  implicit none

  private

  public :: solve_mg_2d

  contains

  subroutine solve_mg_2d(mgphi, f, tol) ! outside facing interface

    real(num), dimension(:,:), allocatable, intent(inout) :: mgphi
    real(num), dimension(1:nx,1:ny), intent(in) :: f
    real(num), intent(in) :: tol
    ! add an argument to select type of MG solver when more than one exists

    call mg_2level(mgphi, f, tol)

    print *, 'called solve_mg_2d'
  end subroutine solve_mg_2d

  ! different solvers

  subroutine mg_2level(mgphi, f, tol)

    real(num), dimension(:,:), allocatable, intent(inout) :: mgphi
    real(num), dimension(1:nx,1:ny), intent(in) :: f
    real(num), intent(in) :: tol

    integer :: nsteps
    integer :: ncoarse = 3
    integer :: nfine = 50
    integer :: c

    real(num), dimension(1:nx,1:ny) :: residual
    real(num) :: L2, L2_old
    real(num) :: Lap
    real(num) :: start
    real(num) :: finish

    ! coarse stuff
    real(num), dimension(1:nx/2,1:ny/2) :: residual_c
    real(num), dimension(:,:), allocatable :: epsi
    integer :: ixc, iyc
    if (.not.allocated(epsi)) allocate(epsi(0:nx/2+1,0:ny/2+1))

    print *, 'called mg_2level'
    call cpu_time (start)

    ! main cycle 
    nsteps = 0
    L2_old = 1e6_num
    do 
      nsteps = nsteps + 1
      ! ncoarse iterations of GS relaxation on finest grid
      do c = 1, ncoarse

        call bcs_fine(mgphi)

        ! redblack, odd
        do iy = 1, ny  
        do ix = 1, nx  
          if (modulo(ix+iy,2) == 1) then
            mgphi(ix,iy) = 0.25_num * ( & 
              & mgphi(ix+1,iy) + mgphi(ix-1,iy) + mgphi(ix,iy+1) + mgphi(ix,iy-1) &
              - dx**2 * f(ix,iy) ) 
          endif
        end do
        end do 
        ! even iteration
        do iy = 1, ny  
        do ix = 1, nx  
          if (modulo(ix+iy,2) == 0) then
            mgphi(ix,iy) = 0.25_num * ( & 
              & mgphi(ix+1,iy) + mgphi(ix-1,iy) + mgphi(ix,iy+1) + mgphi(ix,iy-1) &
              - dx**2 * f(ix,iy) ) 
          endif
        end do
        end do 
      enddo
      
      ! see if converged to tolerance and exit if satisfied

      call bcs_fine(mgphi)

      L2 = 0.0_num
      do iy = 1, ny
      do ix = 1, nx
        Lap = (mgphi(ix+1,iy) - 2.0_num*mgphi(ix,iy) + mgphi(ix-1,iy)) / dx**2 + & 
            & (mgphi(ix,iy+1) - 2.0_num*mgphi(ix,iy) + mgphi(ix,iy-1)) / dy**2 
        residual(ix,iy) = Lap - f(ix,iy)
      enddo
      enddo
      L2 = sqrt(sum(abs(residual)**2)/real(nx*ny,num))

      if (abs(L2-L2_old) <= tol) exit
      L2_old = L2

! skip for rest purposes
!GOTO 9999

      ! restrict the residue down to a grid with half the resolution in each
      ! dimension  
      iy=  1
      do iyc = 1, ny/2
      ix = 1
      do ixc = 1, nx/2
        residual_c(ixc,iyc) = 0.25_num * (residual(ix,iy) + residual(ix+1,iy) &
          & + residual(ix,iy+1) + residual(ix+1,iy+1) ) 
        ix = ix + 2
      enddo
      iy = iy + 2
      enddo

      ! perform smoothing on the coarser grid for Lap * epsi = residual_c

      do c = 1, nfine

        call bcs_coarse(epsi)

        ! redblack, odd
        do iy = 1, ny/2 
        do ix = 1, nx/2
          if (modulo(ix+iy,2) == 1) then
            epsi(ix,iy) = 0.25_num * ( & 
              & epsi(ix+1,iy) + epsi(ix-1,iy) + epsi(ix,iy+1) + epsi(ix,iy-1) &
              - (2.0_num*dx)**2 * residual_c(ix,iy) ) 
          endif
        end do
        end do 
        ! even iteration
        do iy = 1, ny/2  
        do ix = 1, nx/2  
          if (modulo(ix+iy,2) == 0) then
            epsi(ix,iy) = 0.25_num * ( & 
              & epsi(ix+1,iy) + epsi(ix-1,iy) + epsi(ix,iy+1) + epsi(ix,iy-1) &
              - (2.0_num*dx)**2 * residual_c(ix,iy) ) 
          endif
        end do
        end do 

      enddo

      ! prolongate via direct injection back into fine residual

      iy=  1
      do iyc = 1, ny/2
      ix = 1
      do ixc = 1, nx/2
        mgphi(ix,iy) = mgphi(ix,iy) - epsi(ixc,iyc)
        mgphi(ix+1,iy) = mgphi(ix+1,iy) - epsi(ixc,iyc)
        mgphi(ix,iy+1) = mgphi(ix,iy+1) - epsi(ixc,iyc)
        mgphi(ix+1,iy+1) = mgphi(ix+1,iy+1) - epsi(ixc,iyc)
        ix = ix + 2
      enddo
      iy = iy + 2
      enddo
9999 CONTINUE

    enddo ! exit main clycle

    call cpu_time (finish)

    print *, 'nsteps', nsteps
    print *,' L2', L2
    print '(" ****** cpu_time: ",f20.3," seconds.")',finish-start


  end subroutine mg_2level

  subroutine bcs_fine(mgphi)

    real(num), dimension(:,:), allocatable, intent(inout) :: mgphi

    mgphi(0,:) = mgphi(nx,:)
    mgphi(nx+1,:) = mgphi(1,:)
    mgphi(:,0) = mgphi(:,ny)
    mgphi(:,ny+1) = mgphi(:,1)

  end subroutine bcs_fine

  subroutine bcs_coarse(epsi)

    real(num), dimension(:,:), allocatable, intent(inout) :: epsi

    epsi(0,:) = epsi(nx/2,:)
    epsi(nx/2+1,:) = epsi(1,:)
    epsi(:,0) = epsi(:,ny/2)
    epsi(:,ny/2+1) = epsi(:,1)

  end subroutine bcs_coarse

end module multigrid_2d
