module shared_data

  implicit none

  integer :: out_unit = 10

  integer, parameter :: num=selected_real_kind(p=15) 
  integer :: nx, ny
  integer :: ix, iy
  integer :: step, nsteps

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
  real(num), dimension(:,:), allocatable :: uhx, vhx, uhy, vhy


  real(num), dimension(:,:), allocatable :: uxl, uxr, vxl, vxr
  real(num), dimension(:,:), allocatable :: uyl, uyr, vyl, vyr

  real(num), dimension(:,:), allocatable :: ua, va
  real(num), dimension(:,:), allocatable :: macu, macv

  real(num), dimension(:,:), allocatable :: ux, vx, uy, vy !u on x face, v on x face..
  real(num), dimension(:,:), allocatable :: divu ! divergence of MAC, cc

  real(num), dimension(:,:), allocatable :: phi !cc 

  real(num), parameter :: pi = 4.0_num * ATAN(1.0_num)

  real(num), dimension(:,:), allocatable :: ustar, vstar


  real(num), dimension(:,:), allocatable :: gradp_x, gradp_y ! pressure gradient, cc

end module shared_data
