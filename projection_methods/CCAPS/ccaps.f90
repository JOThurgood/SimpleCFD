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

  do
    step = step + 1
    if ((step > nsteps .and. nsteps >= 0) .or. (time >= t_end)) exit

    print *,'******************************************************************'
    print *,'Cycle: ', step, 'Time: ',time
    print *,'******************************************************************'

    call advance_dt
    if ( step ==0 ) call bootstrap 
!    call test_minion
  
    !if (modulo(step,10) ==0) call sln_plots
!    if (modulo(step,10) ==0) call minion_plots

  enddo 

!  call minion_plots
 call sln_plots

!    print *,'warning bootstrap turned off'
  print *, 'CCAPS Terminated Normally'
end program ccaps
