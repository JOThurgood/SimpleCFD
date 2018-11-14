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
    t_end = 1.0_num

    nsteps = -1 !set to <0 to run till t_end 

    use_minmod = .false.    

    use_viscosity = .false.

    visc = 0.0_num

    bc_xmin = periodic
    bc_xmax = periodic
    bc_ymin = periodic
    bc_ymax = periodic

    ! set one of these to true to overwrite everything except nx 
    ! in control with correct setup for the test problems.
    !
    ! This will also cause any nitial_condition set to be ignored in
    ! favour of the specific problems

!    shear_test = .true.
    minion_test = .true. 
!    vortex1_test = .true.
!    drivenlid_test = .true.
    ! DONT set more than one of the above true

  end subroutine user_control

  subroutine shear_test_control
    x_min = 0.0_num
    x_max = 1.0_num
    y_min = x_min
    y_max = x_max
    !nx = 32   allow it to use whatever was set for easy overwriting in
    ! user control. Dont change it here in practice.
    ny = nx !but do force nx=ny
    CFL = 0.8_num
    t_end = 1.0_num 
    nsteps = -1 
    use_minmod = .false.    
    use_viscosity = .false.
    visc = 0.0_num
    bc_xmin = periodic
    bc_xmax = periodic
    bc_ymin = periodic
    bc_ymax = periodic
    dumpfreq = 10
  end subroutine shear_test_control

  subroutine minion_test_control
    x_min = 0.0_num
    x_max = 1.0_num
    y_min = x_min
    y_max = x_max
    !nx = 32   allow it to use whatever was set for easy overwriting in
    ! user control. Dont change it here in practice.
    ny = nx !but do force nx=ny
    CFL = 0.5_num
    t_end = 0.5_num 
    nsteps = -1 
    use_minmod = .false.    
    use_viscosity = .false.
    visc = 0.0_num
    bc_xmin = periodic
    bc_xmax = periodic
    bc_ymin = periodic
    bc_ymax = periodic
    dumpfreq = 10000
  end subroutine minion_test_control

  subroutine vortex1_test_control
    x_min = 0.0_num
    x_max = 1.0_num
    y_min = x_min
    y_max = x_max
    !nx = 32   allow it to use whatever was set for easy overwriting in
    ! user control. Dont change it here in practice.
    ny = nx !but do force nx=ny
    CFL = 0.8_num
    t_end = 40.0_num 
    nsteps = -1
    use_minmod = .false.    
    use_viscosity = .true.
    visc = 1e-3_num
    bc_xmin = no_slip
    bc_xmax = no_slip
    bc_ymin = no_slip
    bc_ymax = no_slip
  end subroutine vortex1_test_control

  subroutine driven_lid_control
    x_min = 0.0_num
    x_max = 1.0_num
    y_min = x_min
    y_max = x_max
    nx = 32 
    ny = nx
    CFL = 1.0_num
    t_end = 100.0_num  
    nsteps = -1
    use_minmod = .false.    
    use_viscosity = .true.
    visc = 1e-3_num
    bc_xmin = no_slip
    bc_xmax = no_slip
    bc_ymin = no_slip
    bc_ymax = driven
    dumpfreq = 100
  end subroutine driven_lid_control


end module control 

