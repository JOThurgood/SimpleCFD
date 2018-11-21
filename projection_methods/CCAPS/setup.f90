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

    time = 0.0_num

    call setup_report

  end subroutine initial_setup

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
    print *,'With gradient limiting options'
    print *,''
    print *,'use_minmod',use_minmod

    print *,''
    print *,'with extra physics options'
    print *,''
    print *,'use_viscosity',use_viscosity
    print *,'use_vardens',use_vardens
    
  end subroutine setup_report

end module setup

