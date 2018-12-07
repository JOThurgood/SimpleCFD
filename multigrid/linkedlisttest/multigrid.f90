module multigrid

  implicit none

  private
  public :: mg_interface

  type grid
    integer :: level 
!    real(num),dimension(:,:), allocatable :: phi
!    real(num),dimension(:,:), allocatable :: f
    type(grid), pointer :: next, prev
!    contains
!      procedure :: allocate_arrs => allocate_arrs
  end type grid

  type(grid), pointer :: head, tail

contains

  subroutine mg_interface
    call initialise_grids
  end subroutine mg_interface

  subroutine create_grid(new_grid)
    type(grid), pointer :: new_grid
    allocate(new_grid)
    nullify(new_grid%next)
    nullify(new_grid%prev)
    new_grid%level = -1
  end subroutine create_grid

  subroutine add_grid(new_grid)
  ! add grid to a list / create a new list if 1st one
    type(grid), pointer :: new_grid
    if (.not. associated(head)) then
      head=>new_grid
      tail=>new_grid
      return
    endif
    tail%next=>new_grid
    new_grid%prev=>tail
    tail=>new_grid
  end subroutine add_grid

  subroutine initialise_grids

    integer :: lev

    type(grid), pointer :: new
    type(grid), pointer :: current
    nullify(new)
    nullify(current)
    nullify(head)
    nullify(tail)

    do lev = 1, 3
      call create_grid(new)
      call add_grid(new)
    enddo

    current =>head
    lev = 0
    do while(associated(current))
      lev = lev + 1
      current%level = lev 
      current=>current%next
    enddo

    current => head
    do while(associated(current))
      print *,current%level
      current=>current%next
    enddo 
 
  end subroutine initialise_grids

end module multigrid

