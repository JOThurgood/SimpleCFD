module initial_conditions

  use shared_data
  
  implicit none

  private 

  public :: set_ic

  contains

  subroutine set_ic

    ! NB: nx, t_end etc is set in control.f90

    vx = 0.0_num
    vy = 0.0_num
    p = 1.0_num

  end subroutine set_ic

  ! Helper subroutines may go here! 

end module initial_conditions

