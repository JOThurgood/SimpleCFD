program test

  use multigrid

  implicit none

  call mg_interface(nx = 16, ny = 16, nlevels = 5)

end program test
