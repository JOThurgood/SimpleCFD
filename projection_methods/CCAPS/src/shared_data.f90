module shared_data

  implicit none

  integer, parameter :: num=selected_real_kind(p=15) 
  integer :: nx, ny
  integer(num) :: step, nsteps

  real(num) :: time = 0.0_num, t_end, dt 
  real(num) :: x_min
  real(num) :: x_max
  real(num) :: y_min
  real(num) :: y_max

  real(num), dimension(:), allocatable :: xc, xb, yc, yb

end module shared_data
