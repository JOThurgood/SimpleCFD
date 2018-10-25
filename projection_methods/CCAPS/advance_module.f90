module advance_module

  use shared_data

  implicit none
  private 

  public :: advance

  contains

  subroutine advance
    print *, 'advance','time',time
  end subroutine advance

end module advance_module

