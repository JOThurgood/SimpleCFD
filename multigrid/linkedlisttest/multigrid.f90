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

  subroutine allocate_arrays(new_grid)
    type(grid), pointer :: new_grid
    allocate(new_grid%phi(0:new_grid%nx+1,0:new_grid%ny+1))
    allocate(new_grid%f(1:new_grid%nx,1:new_grid%ny))
  end subroutine allocate_arrays

  subroutine grid_report(this)
    type(grid), pointer :: this
      print *,'******'
      print *,'level',this%level
      print *,'nx',this%nx
      print *,'ny',this%ny
      print *,'lbound phi',lbound(this%phi)
      print *,'ubound phi',ubound(this%phi)
      print *,'lbound f',lbound(this%f)
      print *,'ubound f',ubound(this%f)
      print *,'******'
  end subroutine grid_report


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

    ! create a linked list of grids with blank / unallocated data
    do lev = 1, nlevels
      call create_grid(new)
      call add_grid(new)
    enddo

    ! go through the list and set up the arrays 
    current =>head
    lev = 0
    do while(associated(current))
      lev = lev + 1
      current%level = lev 
      current%nx = nx / (2**(lev-1)) 
      current%ny = ny / (2**(lev-1)) 
      call allocate_arrays(current)
      current=>current%next
    enddo

    ! cycle through the list for a report to check all is set up good
    current => head
    do while(associated(current))
      call grid_report(current)
      current=>current%next
    enddo 
 
  end subroutine initialise_grids

end module multigrid

