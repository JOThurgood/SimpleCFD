module shared_data

  implicit none

  integer, parameter :: num=selected_real_kind(p=15) 
  integer :: nx, ny
  integer(num) :: step, nsteps

  real(num) :: time = 0.0_num, t_end, dt 

end module shared_data
