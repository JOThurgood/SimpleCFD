! CCAPS - Cell-centered approximate projection solver

program ccaps

  use shared_data
  use setup
  use advance_module

  implicit none 

  call initial_setup

  ! do some sort of bootstrap here? 
  step = step + 1
  call advance_dt

  print *, 'Code Terminated Normally'
end program ccaps
