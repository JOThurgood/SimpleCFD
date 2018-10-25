module advance_module

  use shared_data

  implicit none

  private 

  public :: advance

  contains

  subroutine advance

    call set_dt 
    print *,'Step', step,'dt', dt
  
  ! Step 1 - calculate 
  end subroutine advance

  subroutine set_dt

    real(num) :: dtx, dty 

    dtx = CFL * dx / maxval(abs(u))
    dty = CFL * dy / maxval(abs(v))
    dt = MIN(dtx,dty)

  end subroutine set_dt


end module advance_module

