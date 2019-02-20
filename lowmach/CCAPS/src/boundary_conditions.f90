module boundary_conditions

  ! all boundary conditions except those for the elliptic solvers
  ! to be kept here

  use shared_data

  implicit none

  private

  public :: velocity_bcs, phi_bcs
  public :: rho_bcs

  contains

  subroutine apply_periodic_cc(arr_cc, boundary)

    real(num), dimension(:,:), allocatable, intent(inout) :: arr_cc
    character(len=5), intent(in) :: boundary
!    if (present(arr_cc)) call bc_sanity_check(arr_cc=arr_cc, di=di, varname='vel_bcs')

    if ((boundary /= 'x_min') .and. (boundary /= 'x_max') .and. &
       (boundary /= 'y_min') .and. (boundary /= 'y_max')) then
      print *,'invalid argument passed to apply_periodic_cc'
      print *,'STOP'
      STOP
    endif

    if (boundary == 'x_min') then
      arr_cc(0,:) = arr_cc(nx,:)
      arr_cc(-1,:) = arr_cc(nx-1,:)
    endif 

    if (boundary == 'x_max') then
      arr_cc(nx+1,:) = arr_cc(1,:)
      arr_cc(nx+2,:) = arr_cc(2,:)
    endif 

    if (boundary == 'y_min') then
      arr_cc(:,0) = arr_cc(:,ny)
      arr_cc(:,-1) = arr_cc(:,ny-1)
    endif 

    if (boundary == 'y_max') then
      arr_cc(:,ny+1) = arr_cc(:,1)
      arr_cc(:,ny+2) = arr_cc(:,2)
    endif 

  end subroutine apply_periodic_cc

  subroutine velocity_bcs(arr_cc, arr_xface, arr_yface, di)

    real(num), dimension(:,:), allocatable, optional, intent(inout) :: arr_cc
    real(num), dimension(:,:), allocatable, optional, intent(inout) :: arr_xface
    real(num), dimension(:,:), allocatable, optional, intent(inout) :: arr_yface

    integer, intent(in) :: di

    real(num) :: const 

    ! Sanity checks

    if (present(arr_cc)) call bc_sanity_check(arr_cc=arr_cc, di=di, varname='vel_bcs')
    if (present(arr_xface)) call bc_sanity_check(arr_xface=arr_xface, di=di, varname='vel_bcs')
    if (present(arr_yface)) call bc_sanity_check(arr_yface=arr_yface, di=di, varname='vel_bcs')

    if (.not.(present(arr_cc)) .and. .not.(present(arr_xface)) &
        .and. .not.(present(arr_yface))) then

      print *,'error: empty call to velocity_bcs'
      STOP
    endif

    ! Periodic 

    if (bc_xmin == periodic) then

      ! at xmin, centered and yface are x coordinates 

      if (present(arr_cc)) then
        call apply_periodic_cc(arr_cc = arr_cc, boundary = 'x_min')
      endif

      if (present(arr_yface)) then
        arr_yface(0,:) = arr_yface(nx,:)
        arr_yface(-1,:) = arr_yface(nx-1,:)
      endif

      ! xface

      if (present(arr_xface)) then
        arr_xface(-1,:) = arr_xface(nx-1,:)
        arr_xface(-2,:) = arr_xface(nx-2,:)
      endif

    endif

    if (bc_xmax == periodic) then

      if (present(arr_cc)) then
        call apply_periodic_cc(arr_cc = arr_cc, boundary = 'x_max')
      endif

      if (present(arr_yface)) then
        arr_yface(nx+1,:) = arr_yface(1,:)
        arr_yface(nx+2,:) = arr_yface(2,:)
      endif

      if (present(arr_xface)) then
        arr_xface(nx+1,:) = arr_xface(1,:)
        arr_xface(nx+2,:) = arr_xface(2,:)
      endif

    endif

    if (bc_ymin == periodic) then

      ! centered and xface are same y coordinates

      if (present(arr_cc)) then
        call apply_periodic_cc(arr_cc = arr_cc, boundary = 'y_min')
      endif

      if (present(arr_xface)) then
        arr_xface(:,0) = arr_xface(:,ny)
        arr_xface(:,-1) = arr_xface(:,ny-1)
      endif

      if (present(arr_yface)) then
        arr_yface(:,-1) = arr_yface(:,ny-1)
        arr_yface(:,-2) = arr_yface(:,ny-2)
      endif

    endif

    if (bc_ymax == periodic) then

      if (present(arr_cc)) then
        call apply_periodic_cc(arr_cc = arr_cc, boundary = 'y_max')
      endif

      if (present(arr_xface)) then
        arr_xface(:,ny+1) = arr_xface(:,1)
        arr_xface(:,ny+2) = arr_xface(:,2)
      endif

      if (present(arr_yface)) then
        arr_yface(:,ny+1) = arr_yface(:,1)
        arr_yface(:,ny+2) = arr_yface(:,2)
      endif

    endif


    ! Encode no-slip and more general dirchlet in one go? 
    ! KISS for now ...

    ! No Slip

    if (bc_xmin == no_slip) then

      ! at xmin, centered and yface are x coordinates 

      if (present(arr_cc)) then !this one passes for v but fails for u? likely a corner issue 
        arr_cc(0,:)  = -arr_cc(1,:) 
        arr_cc(-1,:) = -arr_cc(2,:)
      endif

      if (present(arr_yface)) then
        arr_yface(0,:)  = -arr_yface(1,:)
        arr_yface(-1,:) = -arr_yface(2,:)
      endif

      ! xface

      if (present(arr_xface)) then
        arr_xface(-1,:) = -arr_xface(1,:)
        arr_xface(-2,:) = -arr_xface(2,:)
      endif

    endif

    if (bc_xmax == no_slip) then

      if (present(arr_cc)) then
        arr_cc(nx+1,:) = -arr_cc(nx,:)
        arr_cc(nx+2,:) = -arr_cc(nx-1,:)
      endif

      if (present(arr_yface)) then
        arr_yface(nx+1,:) = - arr_yface(nx,:)
        arr_yface(nx+2,:) = - arr_yface(nx-1,:)
      endif

      if (present(arr_xface)) then
        arr_xface(nx+1,:) = - arr_xface(nx-1,:)
        arr_xface(nx+2,:) = - arr_xface(nx-2,:)
      endif

    endif

    if (bc_ymin == no_slip) then

      ! centered and xface are same y coordinates

      if (present(arr_cc)) then
        arr_cc(:,0) = - arr_cc(:,1)
        arr_cc(:,-1) = -arr_cc(:,2)
      endif

      if (present(arr_xface)) then
        arr_xface(:,0) = -arr_xface(:,1)
        arr_xface(:,-1) = -arr_xface(:,2)
      endif

      if (present(arr_yface)) then
        arr_yface(:,-1) = -arr_yface(:,1)
        arr_yface(:,-2) = -arr_yface(:,2)
      endif

    endif

    if (bc_ymax == no_slip) then

      if (present(arr_cc)) then
        arr_cc(:,ny+1) = -arr_cc(:,ny)
        arr_cc(:,ny+2) = -arr_cc(:,ny-1)
      endif

      if (present(arr_xface)) then
        arr_xface(:,ny+1) = -arr_xface(:,ny)
        arr_xface(:,ny+2) = -arr_xface(:,ny-1)
      endif

      if (present(arr_yface)) then
        arr_yface(:,ny+1) = -arr_yface(:,ny-1)
        arr_yface(:,ny+2) = -arr_yface(:,ny-2)
      endif

    endif


    ! Slip - Only for ymin currently

    if (bc_ymin == slip) then

      if (di == 1) then !apply zero gradient to x velocity (i.e. tangent to this face)

        if (present(arr_cc)) then
          arr_cc(:,0) =  arr_cc(:,1)
          arr_cc(:,-1) = arr_cc(:,2)
        endif
  
        if (present(arr_xface)) then
          arr_xface(:,0) = arr_xface(:,1)
          arr_xface(:,-1) = arr_xface(:,2)
        endif
  
        if (present(arr_yface)) then
          arr_yface(:,-1) = arr_yface(:,1)
          arr_yface(:,-2) = arr_yface(:,2)
        endif

      endif ! di == 1

      if (di == 2) then !apply anti-sym to y velocity (i.e. normal to this face)
  
        ! centered and xface are same y coordinates
  
        if (present(arr_cc)) then
          arr_cc(:,0) = - arr_cc(:,1)
          arr_cc(:,-1) = -arr_cc(:,2)
        endif
  
        if (present(arr_xface)) then
          arr_xface(:,0) = -arr_xface(:,1)
          arr_xface(:,-1) = -arr_xface(:,2)
        endif
  
        if (present(arr_yface)) then
          arr_yface(:,-1) = -arr_yface(:,1)
          arr_yface(:,-2) = -arr_yface(:,2)
        endif
  
      endif ! di == 2

    endif



    ! Constant / Dirichlet

    ! Driven - hardcoded as u = 1 for now
    ! problem is distinguishing between u = 1 and v = 0
    
    ! It would be better to handle a class of 'constant' or 'dirchlet' 
    ! and then have the user specify the constant for each dimension

    if ( bc_ymax == dirichlet) then

      if (di == 1) then
        const = drive_vel
      else if (di == 2) then
        const = 0.0_num
      else 
        STOP
      endif

      if (present(arr_cc)) then
        arr_cc(:,ny+1) = 2.0_num * const - arr_cc(:,ny)
        arr_cc(:,ny+2) = 2.0_num * const - arr_cc(:,ny-1)
      endif

      if (present(arr_xface)) then
        arr_xface(:,ny+1) = 2.0_num * const - arr_xface(:,ny)
        arr_xface(:,ny+2) = 2.0_num * const - arr_xface(:,ny-1)
      endif

      if (present(arr_yface)) then
        arr_yface(:,ny+1) = 2.0_num * const - arr_yface(:,ny-1)
        arr_yface(:,ny+2) = 2.0_num * const - arr_yface(:,ny-2)
      endif

    endif


    ! outflow - zero gradient on v and rho, fixed phi

    if (bc_ymax == outflow) then

      if (present(arr_cc)) then
        arr_cc(:,ny+1) = arr_cc(:,ny)
        arr_cc(:,ny+2) = arr_cc(:,ny-1)
      endif

      if (present(arr_xface)) then
        arr_xface(:,ny+1) = arr_xface(:,ny)
        arr_xface(:,ny+2) = arr_xface(:,ny-1)
      endif

      if (present(arr_yface)) then
        arr_yface(:,ny+1) = arr_yface(:,ny-1)
        arr_yface(:,ny+2) = arr_yface(:,ny-2)
      endif

    endif


  end subroutine velocity_bcs


  subroutine phi_bcs

    !Periodic 

    if (bc_xmin == periodic) then
      call apply_periodic_cc(arr_cc = phi, boundary = 'x_min')
    endif
    if (bc_xmax == periodic) then
      call apply_periodic_cc(arr_cc = phi, boundary = 'x_max')
    endif
    if (bc_ymin == periodic) then
      call apply_periodic_cc(arr_cc = phi, boundary = 'y_min')
    endif
    if (bc_ymax == periodic) then
      call apply_periodic_cc(arr_cc = phi, boundary = 'y_max')
    endif

    ! Zero gradient

    if (bc_xmin == zero_gradient) then
      phi(0,:) = phi(1,:)
      phi(-1,:) = phi(2,:)
    endif
    if (bc_xmax == zero_gradient) then
      phi(nx+1,:) = phi(nx,:)
      phi(nx+2,:) = phi(nx-1,:)
    endif
    if (bc_ymin == zero_gradient) then
      phi(:,0) = phi(:,1)
      phi(:,-1) = phi(:,2)
    endif
    if (bc_ymax == zero_gradient) then
      phi(:,ny+1) = phi(:,ny)
      phi(:,ny+2) = phi(:,ny-1)
    endif

    ! No slip should be Neumann for pressure?

    if (bc_xmin == no_slip) then
      phi(0,:) = phi(1,:)
      phi(-1,:) = phi(2,:)
    endif
    if (bc_xmax == no_slip) then
      phi(nx+1,:) = phi(nx,:)
      phi(nx+2,:) = phi(nx-1,:)
    endif
    if (bc_ymin == no_slip) then
      phi(:,0) = phi(:,1)
      phi(:,-1) = phi(:,2)
    endif
    if (bc_ymax == no_slip) then
      phi(:,ny+1) = phi(:,ny)
      phi(:,ny+2) = phi(:,ny-1)
    endif

    ! want to overhaul BC so you don't have so many different things doing the same thing

    if (bc_ymax == dirichlet) then
      phi(:,ny+1) =  phi(:,ny)
      phi(:,ny+2) =  phi(:,ny-1)
    endif

    ! overwrite with experimental no_slip (for atmospheres)

    if (.not. drivenlid_test) then

    if (bc_ymin == no_slip) then
      phi(:,0) =  phi(:,1) - (phi(:,2)-phi(:,1))
      phi(:,-1) = phi(:,1) - (phi(:,2)-phi(:,1))*2.0_num
    endif

    if (bc_ymax == no_slip) then
      phi(:,ny+1) = phi(:,ny) + (phi(:,ny) - phi(:,ny-1))
      phi(:,ny+2) = phi(:,ny) + (phi(:,ny) - phi(:,ny-1))*2.0_num
    endif

    endif

    if (bc_ymax == outflow) then
      phi(:,ny+1) = phi(:,ny) + (phi(:,ny) - phi(:,ny-1))
      phi(:,ny+2) = phi(:,ny) + (phi(:,ny) - phi(:,ny-1))*2.0_num
    endif

    ! slip BC, only for ymin currently

    if (bc_ymin == slip) then
      phi(:,0) =  phi(:,1) - (phi(:,2)-phi(:,1))
      phi(:,-1) = phi(:,1) - (phi(:,2)-phi(:,1))*2.0_num
    endif


  end subroutine phi_bcs

  subroutine rho_bcs(arr_cc, arr_xface, arr_yface)
 
    real(num), dimension(:,:), allocatable, optional, intent(inout) :: arr_cc
    real(num), dimension(:,:), allocatable, optional, intent(inout) :: arr_xface
    real(num), dimension(:,:), allocatable, optional, intent(inout) :: arr_yface
    integer :: di = 1 ! not meaningful for rho (scalar), but set to something
      ! to stop bc_sanity_check from complaining for now

   ! Sanity checks

    if (present(arr_cc)) call bc_sanity_check(arr_cc=arr_cc, di=di, varname='rho_bcs')
    if (present(arr_xface)) call bc_sanity_check(arr_xface=arr_xface, di=di, varname='rho_bcs')
    if (present(arr_yface)) call bc_sanity_check(arr_yface=arr_yface, di=di, varname='rho_bcs')

    if (.not.(present(arr_cc)) .and. .not.(present(arr_xface)) &
        .and. .not.(present(arr_yface))) then

      print *,'error: empty call to rho_bcs'
      STOP
    endif


    ! periodic

    if (bc_xmin == periodic) then

      if (present(arr_cc)) then
!        arr_cc( 0,:) = arr_cc(nx,:)
!        arr_cc(-1,:) = arr_cc(nx-1,:)
        call apply_periodic_cc(arr_cc = arr_cc, boundary = 'x_min')
      endif 

      if (present(arr_xface)) then
        arr_xface(-1,:) = arr_xface(nx-1,:)
        arr_xface(-2,:) = arr_xface(nx-2,:)
      endif 

      if (present(arr_yface)) then
        arr_yface( 0,:) = arr_yface(nx,:)
        arr_yface(-1,:) = arr_yface(nx-1,:)
      endif 

    endif 


    if (bc_xmax == periodic) then

      if (present(arr_cc)) then
!       arr_cc(nx+1,:) = arr_cc(1,:)
!       arr_cc(nx+2,:) = arr_cc(2,:)
        call apply_periodic_cc(arr_cc = arr_cc, boundary = 'x_max')
      endif 

      if (present(arr_xface)) then
        arr_xface(nx+1,:) = arr_xface(1,:)
        arr_xface(nx+2,:) = arr_xface(2,:)
      endif 

      if (present(arr_yface)) then
        arr_yface(nx+1,:) = arr_yface(1,:)
        arr_yface(nx+2,:) = arr_yface(2,:)
      endif 

    endif 

    if (bc_ymin == periodic) then

      if (present(arr_cc)) then
!        arr_cc(:, 0) = arr_cc(:,ny)
!        arr_cc(:,-1) = arr_cc(:,ny-1)
        call apply_periodic_cc(arr_cc = arr_cc, boundary = 'y_min')
      endif

      if (present(arr_xface)) then
        arr_xface(:, 0) = arr_xface(:,ny)
        arr_xface(:,-1) = arr_xface(:,ny-1)
      endif

      if (present(arr_yface)) then
        arr_yface(:,-1) = arr_yface(:,ny-1)
        arr_yface(:,-2) = arr_yface(:,ny-2)
      endif

    endif 

    if (bc_ymax == periodic) then
      if (present(arr_cc)) then
!        arr_cc(:,ny+1) = arr_cc(:,1)
!        arr_cc(:,ny+2) = arr_cc(:,2)
        call apply_periodic_cc(arr_cc = arr_cc, boundary = 'y_max')
      endif
      if (present(arr_xface)) then
        arr_xface(:,ny+1) = arr_xface(:,1)
        arr_xface(:,ny+2) = arr_xface(:,2)
      endif
      if (present(arr_yface)) then
        arr_yface(:,ny+1) = arr_yface(:,1)
        arr_yface(:,ny+2) = arr_yface(:,2)
      endif
    endif 

    ! No slip is just even sym for rho

    if (bc_xmin == no_slip) then

      if (present(arr_cc)) then
        arr_cc( 0,:) = arr_cc(1,:)
        arr_cc(-1,:) = arr_cc(2,:)
      endif 

      if (present(arr_xface)) then
        arr_xface(-1,:) = arr_xface(1,:)
        arr_xface(-2,:) = arr_xface(2,:)
      endif 

      if (present(arr_yface)) then
        arr_yface( 0,:) = arr_yface(1,:)
        arr_yface(-1,:) = arr_yface(2,:)
      endif 

    endif 


    if (bc_xmax == no_slip) then

      if (present(arr_cc)) then
        arr_cc(nx+1,:) = arr_cc(nx,:)
        arr_cc(nx+2,:) = arr_cc(nx-1,:)
      endif 

      if (present(arr_xface)) then
        arr_xface(nx+1,:) = arr_xface(nx-1,:)
        arr_xface(nx+2,:) = arr_xface(nx-2,:)
      endif 

      if (present(arr_yface)) then
        arr_yface(nx+1,:) = arr_yface(nx,:)
        arr_yface(nx+2,:) = arr_yface(nx-1,:)
      endif 

    endif 

    if (bc_ymin == no_slip) then

      if (present(arr_cc)) then
        arr_cc(:, 0) = arr_cc(:,1)
        arr_cc(:,-1) = arr_cc(:,2)
      endif

      if (present(arr_xface)) then
        arr_xface(:, 0) = arr_xface(:,1)
        arr_xface(:,-1) = arr_xface(:,2)
      endif

      if (present(arr_yface)) then
        arr_yface(:,-1) = arr_yface(:,1)
        arr_yface(:,-2) = arr_yface(:,2)
      endif

    endif 

    if (bc_ymax == no_slip) then
      if (present(arr_cc)) then
        arr_cc(:,ny+1) = arr_cc(:,ny)
        arr_cc(:,ny+2) = arr_cc(:,ny-1)
      endif
      if (present(arr_xface)) then
        arr_xface(:,ny+1) = arr_xface(:,ny)
        arr_xface(:,ny+2) = arr_xface(:,ny-1)
      endif
      if (present(arr_yface)) then
        arr_yface(:,ny+1) = arr_yface(:,ny-1)
        arr_yface(:,ny+2) = arr_yface(:,ny-2)
      endif
    endif 


    ! Zero gradient

    ! Driven / general dirichlet 

    ! outflow - zero gradient on v and rho, fixed phi

    if (bc_ymax == outflow) then

      if (present(arr_cc)) then
        arr_cc(:,ny+1) = arr_cc(:,ny)
        arr_cc(:,ny+2) = arr_cc(:,ny-1)
      endif

      if (present(arr_xface)) then
        arr_xface(:,ny+1) = arr_xface(:,ny)
        arr_xface(:,ny+2) = arr_xface(:,ny-1)
      endif

      if (present(arr_yface)) then
        arr_yface(:,ny+1) = arr_yface(:,ny-1)
        arr_yface(:,ny+2) = arr_yface(:,ny-2)
      endif

    endif
    
    ! slip - currenly only for ymin

    if (bc_ymin == slip) then

      if (present(arr_cc)) then
        arr_cc(:, 0) = arr_cc(:,1)
        arr_cc(:,-1) = arr_cc(:,2)
      endif

      if (present(arr_xface)) then
        arr_xface(:, 0) = arr_xface(:,1)
        arr_xface(:,-1) = arr_xface(:,2)
      endif

      if (present(arr_yface)) then
        arr_yface(:,-1) = arr_yface(:,1)
        arr_yface(:,-2) = arr_yface(:,2)
      endif

    endif 

  end subroutine rho_bcs

  subroutine bc_sanity_check(arr_cc, arr_xface, arr_yface,di,varname)

    real(num), dimension(:,:), allocatable, optional, intent(inout) :: arr_cc
    real(num), dimension(:,:), allocatable, optional, intent(inout) :: arr_xface
    real(num), dimension(:,:), allocatable, optional, intent(inout) :: arr_yface
    integer, intent(in) :: di
    character(len=7) :: varname

    ! check component is specified
    ! you should change this so x is 1, y is 2 consistent with fortran dimensions
    if ( (di /= 1) .and. (di /=2) ) then 
      print *,'di not given in    ', varname
      STOP
    endif
    
    ! check at least valid one argument is present
    if ( (.not. present(arr_cc)) .and. &
      &  (.not. present(arr_xface)) .and. &
      &  (.not. present(arr_yface)) ) &
    then
      print *, 'bad call to    ',varname,'    no valid arguments'
      STOP
    endif

    ! check the shape and bounds of the arguments are as expected

    if ( present(arr_cc) ) then
      ! print *, shape(arr_cc)
      if ( (lbound(arr_cc,1) /= -1) .or. &
           (lbound(arr_cc,2) /= -1) .or. &
           (ubound(arr_cc,1) /= nx+2) .or. &
           (ubound(arr_cc,2) /= ny+2) ) &
      then
        print *, 'bad bounds on arr_cc passed to     ',varname
        print *,'lbound(arr_cc) =',lbound(arr_cc)
        print *,'ubound(arr_cc) =',ubound(arr_cc)
      endif
           

      if ((size(arr_cc,1) /= nx+4) .or. (size(arr_cc,2) /= ny+4 ) )then
        print *, 'bad sized arr_cc to ',varname
        STOP
      endif
    endif

    if ( present(arr_xface) ) then
      ! print *, shape(arr_xface)
      if ( (lbound(arr_xface,1) /= -2) .or. &
           (lbound(arr_xface,2) /= -1) .or. &
           (ubound(arr_xface,1) /= nx+2) .or. &
           (ubound(arr_xface,2) /= ny+2) ) &
      then
        print *, 'bad bounds on arr_xface passed to     ',varname
        print *,'lbound(arr_xface) =',lbound(arr_xface)
        print *,'ubound(arr_xface) =',ubound(arr_xface)
      endif
           

      if ((size(arr_xface,1) /= nx+5) .or. (size(arr_xface,2) /= ny+4 ) )then
        print *, 'bad sized arr_xface to ',varname
        STOP
      endif
    endif

    if ( present(arr_yface) ) then
      ! print *, shape(arr_yface)
      if ( (lbound(arr_yface,1) /= -1) .or. &
           (lbound(arr_yface,2) /= -2) .or. &
           (ubound(arr_yface,1) /= nx+2) .or. &
           (ubound(arr_yface,2) /= ny+2) ) &
      then
        print *, 'bad bounds on arr_yface passed to     ',varname
        print *,'lbound(arr_yface) =',lbound(arr_yface)
        print *,'ubound(arr_yface) =',ubound(arr_yface)
      endif
           

      if ((size(arr_yface,1) /= nx+4) .or. (size(arr_yface,2) /= ny+5 ) )then
        print *, 'bad sized arr_yface to ',varname
        STOP
      endif

    endif


  end subroutine bc_sanity_check

end module boundary_conditions

