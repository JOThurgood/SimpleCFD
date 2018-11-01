! CCAPS - Cell-centered approximate projection solver

program ccaps

  use shared_data
  use setup
  use advance_module
  use diagnostics

  implicit none 

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

  enddo 

  call analytic_sln_plots

  print *, 'Code Terminated Normally'
end program ccaps
