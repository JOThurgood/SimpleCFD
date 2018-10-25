module advance_module

  use shared_data

  implicit none

  private 

  public :: advance

  contains

  subroutine advance

    call set_dt 
    print *,'Step', step,'dt', dt

    ! 1 - Calculate the advective velocities

    ! 1A -  calculate  time-centered interface states for the normal 
    ! velocities
    ! (list var names here once you know exactly what is needed)

    do iy = -1, ny+2 ! limits for values on x boundaries, y center
    do ix = -2, nx+2 ! i.e. x faces
      du = 
        ! ^ slope of u at relative cc, (should be limited really)
      uhxl(ix,iy) = u(ix-1,iy) + & 
        & 0.5_num * (1.0_num - dt * u(ix-1,iy) / dx ) * du
      du = 
      uhxr(ix,iy) = u(ix,iy) - & 
        & 0.5_num * (1.0_num + dt * u(ix,iy) / dx) * du 
      dv = 
      vhxl(ix,iy) = v(ix-1,iy) + & 
        & 0.5_num * (1.0_num - dt * u(ix-1,iy) / dx ) * dv
      dv = 
      vhxr(ix,iy) = v(ix,iy) - & 
        & 0.5_num * (1.0_num + dt * u(ix,iy) / dx) * dv 

    end do
    end do

  end subroutine advance

  subroutine set_dt

    real(num) :: dtx, dty 

    dtx = CFL * dx / maxval(abs(u))
    dty = CFL * dy / maxval(abs(v))
    dt = MIN(dtx,dty)

  end subroutine set_dt


end module advance_module

