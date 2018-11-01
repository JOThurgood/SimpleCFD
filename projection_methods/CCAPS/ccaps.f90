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

  print *,nsteps 

  do
    step = step + 1
    if ((step > nsteps .and. nsteps >= 0) .or. (time >= t_end)) exit

    print *,'******************************************************************'
    print *,'Cycle: ', step, 'Time: ',time
    print *,'******************************************************************'

    call advance_dt
    if ( step ==0 ) call bootstrap 

!    call test_analytic_sln !compare against Morrisons analytic sln (only for the relevant IC!)
  
!    if (modulo(step,10) ==0) call sln_plots

  enddo 

  call sln_plots


  print *, 'Code Terminated Normally'
end program ccaps
