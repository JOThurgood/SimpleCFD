module advance_module

  use shared_data
  use boundary_conditions

  implicit none

  private 

  public :: advance

  contains

  subroutine advance

    real(num) :: du, dv

    call set_dt 
    print *,'Step', step,'dt', dt

    call velocity_bcs
    ! 1 - Calculate the advective velocities

    ! 1A -  calculate  time-centered interface states for the normal 
    ! velocities
    ! (list var names here once you know exactly what is needed)

    ! be careful with indicies (conversion between i+1/2 type notation 
    ! and separate indicies for xc and xb values)

!    do iy =  1, ny ! limits for values on x boundaries, y center
!    do ix =  0, nx ! i.e. x faces
!
!      du = minmod((u(ix-1,iy)-u(ix-2,iy))/dx, (u(ix,iy)-u(ix-1,iy))/dx)
!
!      uhxl(ix,iy) = u(ix-1,iy) + & 
!        & 0.5_num * (1.0_num - dt * u(ix-1,iy) / dx ) * du
!      du = 0.0_num
!      uhxr(ix,iy) = u(ix,iy) - & 
!        & 0.5_num * (1.0_num + dt * u(ix,iy) / dx) * du 
!      dv = 0.0_num
!      vhxl(ix,iy) = v(ix-1,iy) + & 
!        & 0.5_num * (1.0_num - dt * u(ix-1,iy) / dx ) * dv
!      dv = 0.0_num
!      vhxr(ix,iy) = v(ix,iy) - & 
!        & 0.5_num * (1.0_num + dt * u(ix,iy) / dx) * dv 
!
!    end do
!    end do

    ! maybe easier to iterate through the cc vars and assign to offset 

    ! x faced data 

    do iy = 1, ny 
    do ix = 0, nx  !xb counts from 0 to nx, <0 and >nx are ghosts 
  
      du = minmod((u(ix,iy)-u(ix-1,iy))/dx,(u(ix+1,iy)-u(ix,iy))/dx)
      uhxl(ix,iy) = u(ix,iy) + &
        0.5_num * (1.0_num - dt * u(ix,iy) / dx ) * du

      du = minmod((u(ix+1,iy)-u(ix,iy))/dx,(u(ix+2,iy)-u(ix+1,iy))/dx)
      uhxr(ix,iy) = u(ix+1,iy) - &
        0.5_num * (1.0_num + dt * u(ix+1,iy) / dx ) * du

      dv = minmod((v(ix,iy)-v(ix-1,iy))/dx, (v(ix+1,iy)-v(ix,iy))/dx) 
      vhxl(ix,iy) = v(ix,iy) + & 
        & 0.5_num * (1.0_num - dt * u(ix,iy) / dx ) * dv

      dv = minmod((v(ix+1,iy)-v(ix,iy))/dx, (v(ix+2,iy)-v(ix+1,iy))/dx) 
      vhxr(ix,iy) = v(ix+1,iy) - & 
        & 0.5_num * (1.0_num + dt * u(ix+1,iy) / dx) * dv 
      
    enddo 
    enddo

    ! y faced data

    do iy = 0, ny 
    do ix = 1, nx

      du = minmod( (u(ix,iy)-u(ix,iy-1))/dy, (u(ix,iy+1)-u(ix,iy))/dy)
      uhyl(ix,iy) = u(ix,iy) + &
        & 0.5_num * (1.0_num - dt * u(ix,iy) / dy) * du

      du = minmod( (u(ix,iy+1)-u(ix,iy))/dy, (u(ix,iy+2)-u(ix,iy+1))/dy)
      uhyr(ix,iy) = u(ix,iy+1) - &
        & 0.5_num * (1.0_num + dt * u(ix,iy+1) / dy) * du

      dv = minmod( (v(ix,iy)-v(ix,iy-1))/dy, (v(ix,iy+1)-v(ix,iy))/dy)
      vhyl(ix,iy) = v(ix,iy) + &
        & 0.5_num * (1.0_num - dt * v(ix,iy) / dy ) * dv

      dv = minmod( (v(ix,iy+1)-v(ix,iy))/dy, (v(ix,iy+2)-v(ix,iy+1))/dy)
      vhyr(ix,iy) = v(ix,iy+1) - &
        & 0.5_num * (1.0_num + dt * v(ix,iy+1) / dy) * dv

    enddo
    enddo



    ! 1B - Use the normal velocities to calculate the advective vel
    ! by solving a Riemann problem


    do iy = 0, ny
    do ix = 0, ny
      if (iy /= 0) then !can do the xface stuff
        uha(ix,iy) = riemann(uhxl(ix,iy),uhxr(ix,iy))
      endif
      if (ix /= 0) then !can do the yface stuff
        vha(ix,iy) = riemann(vhyl(ix,iy), vhyr(ix,iy)) 
      endif
    enddo
    enddo

    ! 1C - Upwind the hat states using the advective vels 

    do iy = 0, ny
    do ix = 0, ny
      if (iy /= 0) then !can do the xface stuff
        uhx(ix,iy) = upwind(uha(ix,iy),uhxl(ix,iy),uhxr(ix,iy))
      endif
      if (ix /= 0) then !can do the yface stuff
      endif
    enddo
    enddo

  
  end subroutine advance

  subroutine set_dt

    real(num) :: dtx, dty 

    dtx = CFL * dx / maxval(abs(u))
    dty = CFL * dy / maxval(abs(v))
    dt = MIN(dtx,dty)

  end subroutine set_dt

  real(num) function minmod(a,b)  
    real(num), intent(in) :: a, b 
    if ( (abs(a) < abs(b)) .and. (a*b > 0.0_num) ) then
      minmod = a        
    else if ( (abs(a) > abs(b)) .and. (a*b > 0.0_num) ) then
      minmod = b        
    else                
      minmod = 0.0_num
    endif               
  end function minmod

  real(num) function riemann(a,b)  
    real(num), intent(in) :: a, b 
    if ( (a > 0.0_num) .and. (a+b>0.0_num) ) then
      riemann = a
    else if ( (a <= 0.0_num) .and. (b >= 0.0_num) ) then
      riemann = 0.0_num
    else 
      riemann = b
    endif 
  end function riemann

  real(num) function upwind(sadv, a, b)
    real(num), intent(in) :: sadv,a,b
    if (sadv > 0.0_num) then
      upwind = a
    else if (sadv < 0.0_num) then
      upwind = b
    else
      upwind = 0.0_num
    endif
  end function upwind

end module advance_module

