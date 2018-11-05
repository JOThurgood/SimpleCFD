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

    !special diagnostics
    if (minion_test) call test_minion

  
  enddo 

  if (minion_test) then 
    call minion_plots
  else 
    call sln_plots
  endif
!    print *,'warning bootstrap turned off'
  print *, 'CCAPS Terminated Normally'
end program ccaps
