! CCAPS - Cell-centered approximate projection solver

program ccaps

  use shared_data
  use setup
  use advance_module
  use diagnostics

  implicit none 

  call initial_setup
  call test_analytic_sln



  ! do some sort of bootstrap here? 
  step = step + 1
  call advance_dt
  call test_analytic_sln
  print *, 'Code Terminated Normally'
end program ccaps
