module control 

  use shared_data

  implicit none

  contains

  subroutine user_control

    x_min = 0.0_num
    x_max = 1.0_num
    y_min = x_min
    y_max = x_max

    nx = 32 
    ny = nx
  
    CFL = 0.8_num
    t_end = 10.0_num

    dt_min = 0.1_num ! if calc dt > this, uses this. Set to huge if you always want dynamically calced

    nsteps = 10 !set to <0 to run till t_end 
    dumpfreq = 1
    use_minmod = .false.    

    bc_xmin = periodic
    bc_xmax = periodic
    bc_ymin = no_slip
    bc_ymax = no_slip

    grav_y  = -1.0_num

    gamma = 5.0_num/3.0_num

    ! uncomment one of these to overwrite user_control
    ! and users initial condition with one of the test cases

!    minion_test = .true. 
!    shear_test = .true.
    rti1_test = .true. 
!    blob1_test = .true. 

! following not yet updated for bkgd
!    vardens_adv_test = .true.

    ! DONT set more than one of the above true

  end subroutine user_control

  subroutine shear_test_control
    x_min = 0.0_num
    x_max = 1.0_num
    y_min = x_min
    y_max = x_max
    nx = 32
    ny = nx !but do force nx=ny
    CFL = 0.8_num
    t_end = 1.0_num 
    dt_snapshot = 0.1_num
    nsteps = -1 
    use_minmod = .false.    
    bc_xmin = periodic
    bc_xmax = periodic
    bc_ymin = periodic
    bc_ymax = periodic
    dumpfreq = -1
    grav_y = 0.0_num
    dt_min = 1e6_num ! i.e. allow it to always be based on u and v and g
  end subroutine shear_test_control

  subroutine minion_test_control
    ! periodic every unit length in each direction - can
    ! use to test Lx /= Ly scenarios also (dx=dy must be maintained though)
    x_min = 0.0_num
    x_max = 1.0_num
    y_min = x_min
    y_max = x_max
    nx = 32
    ny = nx !but do force nx=ny
    CFL = 0.5_num
    t_end = 0.5_num 
    nsteps = -1 
    use_minmod = .false.    
    bc_xmin = periodic
    bc_xmax = periodic
    bc_ymin = periodic
    bc_ymax = periodic
    dumpfreq = 10000
    grav_y = 0.0_num
    dt_min = 1e6_num
  end subroutine minion_test_control

  subroutine vardens_adv_test_control
    x_min = 0.0_num
    x_max = 1.0_num
    y_min = x_min
    y_max = x_max
    !nx = 32   allow it to use whatever was set for easy overwriting in
    ! user control. Dont change it here in practice.
    ny = nx !but do force nx=ny
    CFL = 1.0_num
    t_end = 1000.0_num
    dt_snapshot = 10.0_num 
    nsteps = -1 
    use_minmod = .false.    
    bc_xmin = periodic
    bc_xmax = periodic
    bc_ymin = periodic
    bc_ymax = periodic
    dumpfreq = -1
  end subroutine vardens_adv_test_control

  subroutine rti1_control
    x_min = -0.5_num
    x_max = 0.5_num
    y_min = -2.0_num
    y_max = 2.0_num
    nx = 16
    ny = 4*nx !but do force nx=ny
    CFL = 0.5_num
    t_end = 10.0_num 
    nsteps = -1
    use_minmod = .false.    
    bc_xmin = periodic
    bc_xmax = periodic
    bc_ymin = slip_hse
    bc_ymax = outflow_hse
    dumpfreq = -1 
    dt_snapshot = 1.0_num
    grav_y = -1.0_num
  end subroutine rti1_control

  subroutine blob1_control
    x_min = 0.0_num
    x_max = 1.0_num
    y_min = x_min
    y_max = x_max
    nx = 512
    ny = nx !but do force nx=ny
    CFL = 0.8_num
    t_end = 3.0_num
    nsteps = -1
    use_minmod = .false.    
    bc_xmin = periodic
    bc_xmax = periodic
    bc_ymin = no_slip
    bc_ymax = outflow !no_slip
    dumpfreq = -1
    dt_snapshot = 0.1_num
  end subroutine blob1_control

end module control 

