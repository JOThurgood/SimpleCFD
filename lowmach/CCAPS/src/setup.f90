module setup

  use control
  use shared_data
  use initial_conditions

  implicit none

  private

  public :: initial_setup, bootstrap

  contains

  subroutine initial_setup

    step = -1 !so that the bootstrap can be in the main loop (+1)

    call user_control ! assign user control parameters

    if (shear_test) call shear_test_control
    if (minion_test) call minion_test_control 
    if (vortex1_test) call vortex1_test_control
    if (drivenlid_test) call driven_lid_control
    if (vardens_adv_test) call vardens_adv_test_control
    if (rti1_test) call rti1_control
    if (blob1_test) call blob1_control

    call allocate_global_arrays

    call set_ic 
    if (shear_test) call shear_test_ic
    if (minion_test) call minion_test_ic
    if (vortex1_test) call vortex1_test_ic
    if (drivenlid_test) call driven_lid_ic
    if (vardens_adv_test) call vardens_adv_test_ic
    if (rti1_test) call rti1_ic
    if (blob1_test) call blob1_ic

    call setup_mg

    time = 0.0_num

    call setup_report

  end subroutine initial_setup

  subroutine setup_mg

    if (bc_xmin == periodic) mg_bc_xmin = 'periodic'
    if (bc_xmax == periodic) mg_bc_xmax = 'periodic'
    if (bc_ymin == periodic) mg_bc_ymin = 'periodic'
    if (bc_ymax == periodic) mg_bc_ymax = 'periodic'

!!!!    if (bc_xmin == zero_gradient) mg_bc_xmin = 'zero_gradient'
!!!!    if (bc_xmax == zero_gradient) mg_bc_xmax = 'zero_gradient'
!!!!    if (bc_ymin == zero_gradient) mg_bc_ymin = 'zero_gradient'
!!!!    if (bc_ymax == zero_gradient) mg_bc_ymax = 'zero_gradient'

    if (bc_xmin == driven) mg_bc_xmin = 'zero_gradient'
    if (bc_xmax == driven) mg_bc_xmax = 'zero_gradient'
    if (bc_ymin == driven) mg_bc_ymin = 'zero_gradient'
    if (bc_ymax == driven) mg_bc_ymax = 'zero_gradient'

    if (bc_xmin == no_slip) mg_bc_xmin = 'zero_gradient'
    if (bc_xmax == no_slip) mg_bc_xmax = 'zero_gradient'
    if (bc_ymin == no_slip) mg_bc_ymin = 'zero_gradient'
    if (bc_ymax == no_slip) mg_bc_ymax = 'zero_gradient'

    if (bc_xmin == slip) mg_bc_xmin = 'zero_gradient'
    if (bc_xmax == slip) mg_bc_xmax = 'zero_gradient'
    if (bc_ymin == slip) mg_bc_ymin = 'zero_gradient'
    if (bc_ymax == slip) mg_bc_ymax = 'zero_gradient'
 
    if (bc_xmin == slip_hse) mg_bc_xmin = 'zero_gradient'
    if (bc_xmax == slip_hse) mg_bc_xmax = 'zero_gradient'
    if (bc_ymin == slip_hse) mg_bc_ymin = 'zero_gradient'
    if (bc_ymax == slip_hse) mg_bc_ymax = 'zero_gradient'

    if (bc_xmin == driven) mg_bc_xmin = 'fixed'
    if (bc_xmax == driven) mg_bc_xmax = 'fixed'
    if (bc_ymin == driven) mg_bc_ymin = 'fixed'
    if (bc_ymax == driven) mg_bc_ymax = 'fixed' ! * questionable? see Drikkas and Rider

    if (bc_xmin == outflow_hse) mg_bc_xmin = 'fixed'
    if (bc_xmax == outflow_hse) mg_bc_xmax = 'fixed'
    if (bc_ymin == outflow_hse) mg_bc_ymin = 'fixed'
    if (bc_ymax == outflow_hse) mg_bc_ymax = 'fixed'

    if (use_vardens) then

      if (bc_xmin == periodic) mg_etabc_xmin = 'periodic'
      if (bc_xmax == periodic) mg_etabc_xmax = 'periodic'
      if (bc_ymin == periodic) mg_etabc_ymin = 'periodic'
      if (bc_ymax == periodic) mg_etabc_ymax = 'periodic'

      if (bc_xmin == no_slip) mg_etabc_xmin = 'zero_gradient'
      if (bc_xmax == no_slip) mg_etabc_xmax = 'zero_gradient'
      if (bc_ymin == no_slip) mg_etabc_ymin = 'zero_gradient'
      if (bc_ymax == no_slip) mg_etabc_ymax = 'zero_gradient'

      if (bc_xmin == slip) mg_etabc_xmin = 'zero_gradient'
      if (bc_xmax == slip) mg_etabc_xmax = 'zero_gradient'
      if (bc_ymin == slip) mg_etabc_ymin = 'zero_gradient'
      if (bc_ymax == slip) mg_etabc_ymax = 'zero_gradient'

      if (bc_xmin == slip_hse) mg_etabc_xmin = 'zero_gradient'
      if (bc_xmax == slip_hse) mg_etabc_xmax = 'zero_gradient'
      if (bc_ymin == slip_hse) mg_etabc_ymin = 'zero_gradient'
      if (bc_ymax == slip_hse) mg_etabc_ymax = 'zero_gradient'

      if (bc_xmin == outflow_hse) mg_etabc_xmin = 'zero_gradient'
      if (bc_xmax == outflow_hse) mg_etabc_xmax = 'zero_gradient'
      if (bc_ymin == outflow_hse) mg_etabc_ymin = 'zero_gradient'
      if (bc_ymax == outflow_hse) mg_etabc_ymax = 'zero_gradient'

    endif ! use_vardens

!!!!!    if (bc_ymax == outflow) then 
!!!!!      mg_bc_ymax = 'fixed'
!!!!!      mg_etabc_ymax = 'zero_gradient'
!!!!!    endif

  end subroutine setup_mg


  subroutine bootstrap

    call set_ic 
    if (shear_test) call shear_test_ic
    if (minion_test) call minion_test_ic
    if (vortex1_test) call vortex1_test_ic
    if (drivenlid_test) call driven_lid_ic
    if (vardens_adv_test) call vardens_adv_test_ic
    if (rti1_test) call rti1_ic
    if (blob1_test) call blob1_ic

    time = 0.0_num !reset to zero important

  endsubroutine bootstrap

  subroutine allocate_global_arrays 

    ! setup grid 

    ! This method requires 2 ghost cells

    allocate(xc(-1:nx+2))  !cell center (cc) - counts from 1 to nx
    allocate(xb(-2:nx+2))  !cell boundary (cb) - counts from 0 to nx
    allocate(yc(-1:ny+2))  !cell center (cc) - counts from 1 to ny
    allocate(yb(-2:ny+2))  !cell boundary (cb) - counts from 0 to ny
                           !all above have two guards  either side
                                                    
    dx = (x_max - x_min) / REAL(nx,num)             
    do ix = -2, nx+2                                
      xb(ix) = REAL(ix,num) * dx                    
      if (ix /= -2) xc(ix) = xb(ix) - dx/2.0_num 
    enddo                                           

    xb = xb + x_min
    xc = xc + x_min !translate 

!    print *,dx
!    print *,xb
!    print *,xc
!    STOP
                                                    
    dy = (y_max - y_min) / REAL(ny,num)             
    do iy = -2, ny+2                                
      yb(iy) = REAL(iy,num) * dy                    
      if (iy /= -2) yc(iy) = yb(iy) - dy/2.0_num 
    enddo                                           

    yb = yb + y_min
    yc = yc + y_min

    ! setup the cell-centered arrays of primative variables

    allocate(u(-1:nx+2,-1:ny+2))
    allocate(v(-1:nx+2,-1:ny+2))
!    allocate(p(-1:nx+2,-1:ny+2))
 
    ! time centered left and right interface states, defined on
    ! the MAC grid



    allocate(uhxl(0:nx,1:ny)) ! uhxl - u hat xface left state
    allocate(uhxr(0:nx,1:ny)) 
    allocate(vhxl(0:nx,1:ny))
    allocate(vhxr(0:nx,1:ny)) 
  
    allocate(uhyl(1:nx,0:ny)) ! uhxl - u hat yface "left" state 
    allocate(uhyr(1:nx,0:ny)) ! (aka "bottom")
    allocate(vhyl(1:nx,0:ny)) 
    allocate(vhyr(1:nx,0:ny)) 

    !uha,vha, uhx, vhx, uhy,vhy all need ghost filling
    allocate(uha(-2:nx+2,-1:ny+2)) ! uhx - u hat advective (x face) 
    allocate(vha(-1:nx+2,-2:ny+2)) ! v hat advective (y face)

    allocate(uhx(-2:nx+2,-1:ny+2)) ! u hat x face at end of 1C 
    allocate(vhx(-2:nx+2,-1:ny+2)) 
    allocate(uhy(-1:nx+2,-2:ny+2)) ! y gace at end of 1C 
    allocate(vhy(-1:nx+2,-2:ny+2)) 


    allocate(uxl(0:nx,1:ny)) ! full prediction for normal vel left state
    allocate(uxr(0:nx,1:ny)) ! in 1D and 3E (some aren't used in
    allocate(vxl(0:nx,1:ny)) !
    allocate(vxr(0:nx,1:ny)) !

    allocate(uyl(1:nx,0:ny))
    allocate(uyr(1:nx,0:ny))
    allocate(vyl(1:nx,0:ny)) 
    allocate(vyr(1:nx,0:ny))

    allocate(ua(0:nx,1:ny)) ! ua - u advective (x face)
    allocate(va(1:nx,0:ny)) ! v  advective (y face)

    allocate(macu(-2:nx+2,-1:ny+2)) ! mac velocities on faces
    allocate(macv(-1:nx+2,-2:ny+2)) 

    allocate(ux(0:nx,1:ny)) ! u on x face
    allocate(uy(1:nx,0:ny)) !u on y face
    allocate(vx(0:nx,1:ny))
    allocate(vy(1:nx,0:ny))

    allocate(divu(1:nx,1:ny)) ! divergence of mac velocities, cell centers, no ghost

    allocate(phi(-1:nx+2,-1:ny+2)) !phi on cc
    phi = 0.0_num !initial guess in step 0

    allocate(ustar(-1:nx+2,-1:ny+2))
    allocate(vstar(-1:nx+2,-1:ny+2))

    allocate(gradp_x(-1:nx+2,-1:ny+2)) !gradp_x
    allocate(gradp_y(-1:nx+2,-1:ny+2)) !gradp_y
    gradp_x = 0.0_num
    gradp_y = 0.0_num

    ! variable density

    if (use_vardens) then

      allocate(rho(-1:nx+2,-1:ny+2))
      allocate(rhohxl(0:nx,1:ny))
      allocate(rhohxr(0:nx,1:ny))
      allocate(rhohyl(1:nx,0:ny))
      allocate(rhohyr(1:nx,0:ny))
      allocate(rhohx(-2:nx+2,-1:ny+2)) !these need ghosts
      allocate(rhohy(-1:nx+2,-2:ny+2))

      allocate(rhoxl(0:nx,1:ny))
      allocate(rhoxr(0:nx,1:ny))
      allocate(rhoyl(1:nx,0:ny))
      allocate(rhoyr(1:nx,0:ny))
      allocate(rhox(-2:nx+2,-1:ny+2)) !these need ghosts
      allocate(rhoy(-1:nx+2,-2:ny+2))

    endif

  end subroutine allocate_global_arrays
  
  subroutine setup_report

    print *,'******************************************************************'
    print *,' Code initialised with the following control variables:'
    print *,''
    print *,'x_min =',x_min 
    print *,'x_max =',x_max
    print *,'y_min =',y_min
    print *,'y_max =',y_max
    print *,'nx =', nx
    print *,'ny =', ny
    print *,'CFL factor CFL = ',CFL
    print *,'Max no of cycles nsteps = ',nsteps
    print *,'Termination time t_end =', t_end
    print *,''
    print *,' This gives the derived grid variables:'
    print *,''
    print *,'dx =', dx
    print *,'dy =', dy

    print *,''
    print *,'With boundary conditions:'
    print *,''
    if (bc_xmin == periodic) then
      print *,'bc_xmin = periodic , (code =',bc_xmin,')'
    else if (bc_xmin == no_slip) then
      print *,'bc_xmin = no_slip, (code =',bc_xmin,')'
    else 
      print *,'bc_xmin = ',bc_xmin,'(parameter number - name unanticipated by setup report'
    endif 
    if (bc_xmax == periodic) then
      print *,'bc_xmax = periodic , (code =',bc_xmax,')'
    else if (bc_xmax == no_slip) then
      print *,'bc_xmax = no_slip, (code =',bc_xmax,')'
    else 
      print *,'bc_xmax = ',bc_xmax,'(parameter number - name unanticipated by setup report'
    endif 
    if (bc_ymin == periodic) then
      print *,'bc_ymin = periodic , (code =',bc_ymin,')'
    else if (bc_ymin == no_slip) then
      print *,'bc_ymin = no_slip, (code =',bc_ymin,')'
    else 
      print *,'bc_ymin = ',bc_ymin,'(parameter number - name unanticipated by setup report'
    endif 
    if (bc_ymax == periodic) then
      print *,'bc_ymax = periodic , (code =',bc_ymax,')'
    else if (bc_ymax == no_slip) then
      print *,'bc_ymax = no_slip, (code =',bc_ymax,')'
    else 
      print *,'bc_ymax = ',bc_ymax,'(parameter number - name unanticipated by setup report'
    endif 

    print *,''
    print *,'With gradient limiting options'
    print *,''
    print *,'use_minmod',use_minmod

    print *,''
    print *,'with extra physics options'
    print *,''
    print *,'use_viscosity',use_viscosity
    print *,'use_vardens',use_vardens
    
    print *,''
    print *,'with gravity options'
    print *,''
    print *,'grav_x', grav_x
    print *,'grav_y', grav_y


    print *,''
    print *,'with multigrid BC options'
    print *,''    
    print *,'mg_bc_xmin', mg_bc_xmin
    print *,'mg_bc_xmax', mg_bc_xmax
    print *,'mg_bc_ymin', mg_bc_ymin
    print *,'mg_bc_ymax', mg_bc_ymax
    print *,''    
    print *,'mg_etabc_xmin', mg_etabc_xmin
    print *,'mg_etabc_xmax', mg_etabc_xmax
    print *,'mg_etabc_ymin', mg_etabc_ymin
    print *,'mg_etabc_ymax', mg_etabc_ymax
    print *,''    
    print *,'mg_etaval_bx_xmin',  mg_etaval_bc_xmin 
    print *,'mg_etaval_bx_xmax',  mg_etaval_bc_xmin 
    print *,'mg_etaval_bx_ymin',  mg_etaval_bc_ymin 
    print *,'mg_etaval_bx_ymax',  mg_etaval_bc_ymax 

  end subroutine setup_report

end module setup

