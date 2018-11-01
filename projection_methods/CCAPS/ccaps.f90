! CCAPS - Cell-centered approximate projection solver

program ccaps

  use shared_data
  use setup
  use advance_module
  use diagnostics
  use welcome

  implicit none 

  call welcome_msg

  call initial_setup
  call test_analytic_sln

  print *,nsteps 

  ! step zero / bootstrapping

  call advance_dt ! to initialise gradp_x, gradp_y 

  call bootstrap
!  call execute_command_line("clear")

  do
    step = step + 1
    if ((step > nsteps .and. nsteps >= 0) .or. (time >= t_end)) exit

    print *,'******************************************************************'
    print *,'Cycle: ', step, 'Time: ',time
    print *,'******************************************************************'

    call advance_dt
    call test_analytic_sln
  
    if (modulo(step,10) ==0) call sln_plots

  enddo 

  call sln_plots


  print *, 'Code Terminated Normally'
end program ccaps
