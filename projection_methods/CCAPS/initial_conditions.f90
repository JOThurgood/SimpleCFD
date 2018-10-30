module initial_conditions

  use shared_data
  
  implicit none

  private 

  public :: set_ic

  contains

  subroutine set_ic

    ! NB: nx, t_end etc is set in control.f90

    real(num) :: x,y

    do iy = -1,ny+1
    do ix = -1,nx+1
      x = xc(ix)
      y = yc(iy)
      u(ix,iy) = 1.0_num - 2.0_num * cos(2.0_num*pi*x)*sin(2.0_num*pi*y)
      v(ix,iy) = 1.0_num + 2.0_num * sin(2.0_num*pi*x)*cos(2.0_num*pi*y)
    enddo
    enddo

!    p = 1.0_num

  end subroutine set_ic

  ! Helper subroutines may go here! 

end module initial_conditions

