module multigrid

  implicit none

  private
  public :: mg_interface

  ! will need own private constants and parameters since entirely modular
  integer, parameter :: num=selected_real_kind(p=15)

  ! data object
  type grid
    integer :: level, nx, ny  
    real(num),dimension(:,:), allocatable :: phi
    real(num),dimension(:,:), allocatable :: f
    type(grid), pointer :: next, prev
!    contains
!      procedure :: allocate_arrs => allocate_arrs
  end type grid

  type(grid), pointer :: head, tail

contains

  subroutine mg_interface(nx,ny,nlevels)
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    integer, intent(in) :: nlevels

    call initialise_grids(nx,ny,nlevels)

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

!  subroutine create_grid(new_grid)
!    type(grid), pointer :: new_grid
!    allocate(new_grid)
!    nullify(new_grid%next)
!    nullify(new_grid%prev)
!    new_grid%level = -1
!  end subroutine create_grid


  subroutine initialise_grids(nx,ny,nlevels)

    integer, intent(in) :: nx
    integer, intent(in) :: ny
    integer, intent(in) :: nlevels

    integer :: lev

    type(grid), pointer :: new
    type(grid), pointer :: current
    nullify(new)
    nullify(current)
    nullify(head)
    nullify(tail)

    do lev = 1, nlevels
      call create_grid(new)
      call add_grid(new)
    enddo

    current =>head
    lev = 0
    do while(associated(current))
      lev = lev + 1
      current%level = lev 
      current%nx = nx / (2**(lev-1)) 
      current%ny = ny / (2**(lev-1)) 
      current=>current%next
    enddo

    current => head
    do while(associated(current))
      print *,'level',current%level
      print *,'nx',current%nx
      print *,'ny',current%ny
      current=>current%next
    enddo 
 
  end subroutine initialise_grids

end module multigrid

