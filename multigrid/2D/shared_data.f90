module shared_data

  implicit none

  ! intrinsics

  integer :: out_unit = 10
  integer, parameter :: num=selected_real_kind(p=15) 
  real(num) :: pi = 4.0_num * ATAN(1.0_num)
 ! core solver
 
  integer :: nx, ny
  integer :: ix, iy

  real(num) :: x_min, x_max
  real(num) :: y_min, y_max
  real(num) :: dx, dy
 
  real(num), dimension(:), allocatable :: xc, yc !xb, yc, yb

  real(num), dimension(:,:), allocatable :: phi

  contains

  subroutine setup

    ! two ghosts are allocated because thats what we have in CCAPS

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
  end subroutine setup


end module shared_data
