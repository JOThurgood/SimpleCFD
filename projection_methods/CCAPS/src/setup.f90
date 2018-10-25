module setup

  use control
  use shared_data

  implicit none

  private

  public :: initial_setup

  contains

  subroutine initial_setup

    call user_control ! assign user control parameters

    print *, nx

  end subroutine initial_setup

end module setup

