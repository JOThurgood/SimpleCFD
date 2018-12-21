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

    nsteps = 10 !set to <0 to run till t_end 
    dumpfreq = 1
    use_minmod = .false.    

    use_viscosity = .false.

    visc = 0.0_num

    bc_xmin = periodic
    bc_xmax = periodic
    bc_ymin = no_slip
    bc_ymax = no_slip

    grav_y  = -1.0_num

    ! set one of these to true to overwrite everything except nx 
    ! in control with correct setup for the test problems.
    !
    ! This will also cause any nitial_condition set to be ignored in
    ! favour of the specific problems

!    shear_test = .true.
    minion_test = .true. 
!    vortex1_test = .true.
!    drivenlid_test = .true.
!!    vardens_adv_test = .true.
!    rti1_test = .true. 
!    blob1_test = .true. 

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
    dumpfreq = -1
    grav_x = 0.0_num
    grav_y = 0.0_num
! use_vardens = .true.
  end subroutine shear_test_control

  subroutine minion_test_control
    ! periodic every unit length in each direction - can
    ! use to test Lx /= Ly scenarios also (dx=dy must be maintained though)
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
    grav_x = 0.0_num
    grav_y = 0.0_num
!    use_vardens = .true.
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
    grav_x = 0.0_num
    grav_y = 0.0_num
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
    bc_ymax = dirichlet ! dirichlet is hardcoded in boundary.f90 to set ux = drive_vel
    drive_vel = 1.0_num
    dumpfreq = -1
    dt_snapshot = 5.0_num
    grav_x = 0.0_num
    grav_y = 0.0_num
    use_vardens = .false.
  end subroutine driven_lid_control

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
    nsteps = -1 
    use_minmod = .false.    
    use_viscosity = .false.
    visc = 0.0_num
    bc_xmin = periodic
    bc_xmax = periodic
    bc_ymin = periodic
    bc_ymax = periodic
    dumpfreq = -1
    use_vardens = .true.
  end subroutine vardens_adv_test_control

  subroutine rti1_control
    x_min = -0.5_num
    x_max = 0.5_num
!    y_min = -0.5_num
!    y_max = 0.5_num
    y_min = -2.0_num
    y_max = 2.0_num
    nx = 16
!    ny = nx
    ny = 4*nx !but do force nx=ny
    CFL = 0.5_num
    t_end = 10.0_num 
    nsteps = -1
    use_minmod = .false.    
    use_viscosity = .false.
    visc = 0.0_num
    bc_xmin = periodic
    bc_xmax = periodic
    bc_ymin = no_slip
    bc_ymax = no_slip
    dumpfreq = -1 
    dt_snapshot = 0.2_num
    use_vardens = .true.
    grav_x = 0.0_num
    grav_y = -1.0_num
  end subroutine rti1_control

  subroutine blob1_control
    x_min = 0.0_num
    x_max = 1.0_num
    y_min = x_min
    y_max = x_max
    nx = 32
    ny = nx !but do force nx=ny
    CFL = 1.0_num
    t_end = 1000.0_num 
    nsteps = 10 
    use_minmod = .false.    
    use_viscosity = .false.
    visc = 0.0_num
    bc_xmin = periodic
    bc_xmax = periodic
    bc_ymin = no_slip
    bc_ymax = no_slip
    dumpfreq = 1
    use_vardens = .true.
  end subroutine blob1_control

end module control 

