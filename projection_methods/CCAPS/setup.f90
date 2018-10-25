module setup

  use control
  use shared_data
  use initial_conditions

  implicit none

  private

  public :: initial_setup

  contains

  subroutine initial_setup

    call user_control ! assign user control parameters
    call allocate_global_arrays
    call set_ic 

  end subroutine initial_setup

  subroutine allocate_global_arrays 

    ! setup grid 

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
                                                    
    dy = (y_max - y_min) / REAL(ny,num)             
    do iy = -2, ny+2                                
      yb(iy) = REAL(iy,num) * dy                    
      if (iy /= -2) yc(iy) = yb(iy) - dy/2.0_num 
    enddo                                           

    ! setup the cell-centered arrays of primative variables

    allocate(u(-1:nx+2,-1:ny+2))
    allocate(v(-1:nx+2,-1:ny+2))
    allocate(p(-1:nx+2,-1:ny+2))
 

  end subroutine allocate_global_arrays

end module setup

