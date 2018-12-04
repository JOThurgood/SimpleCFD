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

  next_dump = 0.0_num !dt_snapshot

  do
    step = step + 1
    if ((step > nsteps .and. nsteps >= 0) .or. (time >= t_end)) exit

    print *,'******************************************************************'
    print *,'Cycle: ', step, 'Time: ',time
    print *,'******************************************************************'
 
    print *, 'rho on grid',sum(rho(1:nx,1:ny)*dx*dy)

    call advance_dt
    if ( step ==0 ) call bootstrap 

    ! periodic dumps 
    if ( (modulo(step,dumpfreq) == 0) .and. (dumpfreq > 0) ) call sln_plots
    if ( time >= next_dump) then
      call sln_plots  
      next_dump = next_dump + dt_snapshot
    endif 

    ! special diagnostics
    if (minion_test) call test_minion

  enddo 

  ! final dump
  call sln_plots
  if (minion_test) call minion_plots

  print *, 'CCAPS Terminated Normally'
end program ccaps
