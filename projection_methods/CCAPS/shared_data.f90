module shared_data

  implicit none

  integer, parameter :: num=selected_real_kind(p=15) 
  integer :: nx, ny
  integer :: ix, iy
  integer(num) :: step, nsteps

  real(num) :: time = 0.0_num, t_end, dt 
  real(num) :: CFL
  real(num) :: x_min, x_max
  real(num) :: y_min, y_max
  real(num) :: dx, dy

  real(num), dimension(:), allocatable :: xc, xb, yc, yb

  real(num), dimension(:,:), allocatable :: u, v, p
  real(num), dimension(:,:), allocatable :: uhxl, uhxr, vhxl, vhxr
  real(num), dimension(:,:), allocatable :: uhyl, uhyr, vhyl, vhyr

  real(num), dimension(:,:), allocatable :: uha, vha

  real(num), dimension(:,:), allocatable :: ul, ur, vl, vr 

end module shared_data
