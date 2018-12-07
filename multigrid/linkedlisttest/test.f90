program test

  use multigrid

  implicit none
  integer, parameter :: num=selected_real_kind(p=15)

  call mg_interface(nx = 16, ny = 16, dx = 1.0_num, dy =1.0_num ,  nlevels = 5)

end program test
