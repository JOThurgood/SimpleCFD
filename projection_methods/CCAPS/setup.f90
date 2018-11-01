module setup

  use control
  use shared_data
  use initial_conditions

  implicit none

  private

  public :: initial_setup, bootstrap

  contains

  subroutine initial_setup

    call user_control ! assign user control parameters
    call allocate_global_arrays
    call set_ic 

    step = 0
    time = 0.0_num

  end subroutine initial_setup

  subroutine bootstrap

    call set_ic 
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

!    print *,dx
!    print *,xb
!    print *,xc
!    STOP
                                                    
    dy = (y_max - y_min) / REAL(ny,num)             
    do iy = -2, ny+2                                
      yb(iy) = REAL(iy,num) * dy                    
      if (iy /= -2) yc(iy) = yb(iy) - dy/2.0_num 
    enddo                                           

    ! setup the cell-centered arrays of primative variables

    allocate(u(-1:nx+2,-1:ny+2))
    allocate(v(-1:nx+2,-1:ny+2))
!    allocate(p(-1:nx+2,-1:ny+2))
 
    ! time centered left and right interface states, defined on
    ! the MAC grid
  
    allocate(uhxl(-2:nx+2,-1:ny+2)) ! uhxl - u hat xface left state
    allocate(uhxr(-2:nx+2,-1:ny+2)) 
    allocate(vhxl(-2:nx+2,-1:ny+2))
    allocate(vhxr(-2:nx+2,-1:ny+2)) 
  
    allocate(uhyl(-1:nx+2,-2:ny+2)) ! uhxl - u hat yface "left" state 
    allocate(uhyr(-1:nx+2,-2:ny+2)) ! (aka "bottom")
    allocate(vhyl(-1:nx+2,-2:ny+2)) 
    allocate(vhyr(-1:nx+2,-2:ny+2)) 

    allocate(uha(-2:nx+2,-1:ny+2)) ! uhx - u hat advective (x face)
    allocate(vha(-1:nx+2,-2:ny+2)) ! v hat advective (y face)

    allocate(uhx(-2:nx+2,-1:ny+2)) ! u hat x face at end of 1C 
    allocate(vhx(-2:nx+2,-1:ny+2)) 
    allocate(uhy(-1:nx+2,-2:ny+2)) ! y gace at end of 1C 
    allocate(vhy(-1:nx+2,-2:ny+2)) 


  

    allocate(uxl(-2:nx+2,-1:ny+2)) ! full prediction for normal vel left state
    allocate(uxr(-2:nx+2,-1:ny+2)) ! in 1D and 3E (some aren't used in
    allocate(vxl(-1:nx+2,-2:ny+2)) !
    allocate(vxr(-1:nx+2,-2:ny+2)) !

    allocate(uyl(-2:nx+2,-1:ny+2)) 
    allocate(uyr(-2:nx+2,-1:ny+2))
    allocate(vyl(-1:nx+2,-2:ny+2)) 
    allocate(vyr(-1:nx+2,-2:ny+2))

    allocate(ua(-2:nx+2,-1:ny+2)) ! ua - u advective (x face)
    allocate(va(-1:nx+2,-2:ny+2)) ! v  advective (y face)

    allocate(macu(-2:nx+2,-1:ny+2)) ! mac velocities on faces
    allocate(macv(-1:nx+2,-2:ny+2)) 

    allocate(ux(-2:nx+2,-1:ny+2)) ! u on x face
    allocate(uy(-1:nx+2,-2:ny+2)) !u on y face
    allocate(vx(-2:nx+2,-1:ny+2))
    allocate(vy(-1:nx+2,-2:ny+2))

    allocate(divu(-1:nx+1,-1:ny+2)) ! divergence of mac velocities, cell centers

    allocate(phi(-1:nx+1,-1:ny+2)) !phi on cc
    phi = 0.0_num !initial guess in step 0

    allocate(ustar(-1:nx+2,-1:ny+2))
    allocate(vstar(-1:nx+2,-1:ny+2))

    allocate(gradp_x(-1:nx+1,-1:ny+2)) !gradp_x
    allocate(gradp_y(-1:nx+1,-1:ny+2)) !gradp_y
    gradp_x = 0.0_num
    gradp_y = 0.0_num


  end subroutine allocate_global_arrays

end module setup

