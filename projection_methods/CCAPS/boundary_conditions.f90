module boundary_conditions

  use shared_data

  implicit none

  private

  public :: velocity_bcs, phi_bcs, gradp_bcs, velocity_face_bcs

  contains

  subroutine velocity_bcs

    ! Periodic

    if (bc_xmin == periodic) then
      u(0,:) = u(nx,:)
      u(-1,:) = u(nx-1,:)
      v(0,:) = v(nx,:)
      v(-1,:) = v(nx-1,:)
      ustar(0,:) = ustar(nx,:)
      ustar(-1,:) = ustar(nx-1,:)
      vstar(0,:) = vstar(nx,:)
      vstar(-1,:) = vstar(nx-1,:)
    endif 

    if (bc_xmax == periodic) then
      u(nx+1,:) = u(1,:)
      u(nx+2,:) = u(2,:)
      v(nx+1,:) = v(1,:)
      v(nx+2,:) = v(2,:)
      ustar(nx+1,:) = ustar(1,:)
      ustar(nx+2,:) = ustar(2,:)
      vstar(nx+1,:) = vstar(1,:)
      vstar(nx+2,:) = vstar(2,:)
    endif

    if (bc_ymin == periodic) then 
      u(:,0) = u(:,ny)
      u(:,-1) = u(:,ny-1)
      v(:,0) = v(:,ny)
      v(:,-1) = v(:,ny-1)
      ustar(:,0) = ustar(:,ny)
      ustar(:,-1) = ustar(:,ny-1)
      vstar(:,0) = vstar(:,ny)
      vstar(:,-1) = vstar(:,ny-1)
    endif

    if (bc_ymax == periodic) then
      u(:,ny+1) = u(:,1)
      u(:,ny+2) = u(:,2)
      v(:,ny+1) = v(:,1)
      v(:,ny+2) = v(:,2)
      ustar(:,ny+1) = ustar(:,1)
      ustar(:,ny+2) = ustar(:,2)
      vstar(:,ny+1) = vstar(:,1)
      vstar(:,ny+2) = vstar(:,2)
    endif 

    ! Zero gradient
    if (bc_xmin == zero_gradient) then
      u(0,:) = u(1,:) 
      u(-1,:) = u(2,:)
      v(0,:) = v(1,:)
      v(-1,:) = v(2,:)
      ustar(0,:) = ustar(1,:)
      ustar(-1,:) = ustar(2,:)
      vstar(0,:) = vstar(1,:)
      vstar(-1,:) = vstar(2,:)
    endif

    if (bc_xmax == zero_gradient) then
      u(nx+1,:) = u(nx,:)
      u(nx+2,:) = u(nx-1,:)
      v(nx+1,:) = v(nx,:)
      v(nx+2,:) = v(nx-1,:)
      ustar(nx+1,:) = ustar(nx,:)
      ustar(nx+2,:) = ustar(nx-1,:)
      vstar(nx+1,:) = vstar(nx,:)
      vstar(nx+2,:) = vstar(nx-1,:)
    endif

    if (bc_ymin == zero_gradient) then
      u(:,0) = u(:,1)
      u(:,-1) = u(:,2)
      v(:,0) = v(:,1)
      v(:,-1) = v(:,2)
      ustar(:,0) = ustar(:,1)
      ustar(:,-1) = ustar(:,2)
      vstar(:,0) = vstar(:,1)
      vstar(:,-1) = vstar(:,2)
    endif

    if (bc_ymax == zero_gradient) then
      u(:,ny+1) = u(:,ny)
      u(:,ny+2) = u(:,ny-1)
      v(:,ny+1) = v(:,ny)
      v(:,ny+2) = v(:,ny-1)
      ustar(:,ny+1) = ustar(:,ny)
      ustar(:,ny+2) = ustar(:,ny-1)
      vstar(:,ny+1) = vstar(:,ny)
      vstar(:,ny+2) = vstar(:,ny-1)
    endif

    ! No slip (u = v = 0 on the face) 

    if (bc_xmin == no_slip) then
      u(0,:) = -u(1,:) 
      u(-1,:) = -u(2,:)
      v(0,:) = -v(1,:)
      v(-1,:) = -v(2,:)
      ustar(0,:) = -ustar(1,:)
      ustar(-1,:) =- ustar(2,:)
      vstar(0,:) = -vstar(1,:)
      vstar(-1,:) = -vstar(2,:)
    endif

    if (bc_xmax == no_slip) then
      u(nx+1,:) = -u(nx,:)
      u(nx+2,:) = -u(nx-1,:)
      v(nx+1,:) = -v(nx,:)
      v(nx+2,:) = -v(nx-1,:)
      ustar(nx+1,:) = -ustar(nx,:)
      ustar(nx+2,:) = -ustar(nx-1,:)
      vstar(nx+1,:) = -vstar(nx,:)
      vstar(nx+2,:) = -vstar(nx-1,:)
    endif

    if (bc_ymin == no_slip) then
      u(:,0) = -u(:,1)
      u(:,-1) = -u(:,2)
      v(:,0) = -v(:,1)
      v(:,-1) = -v(:,2)
      ustar(:,0) = -ustar(:,1)
      ustar(:,-1) = -ustar(:,2)
      vstar(:,0) = -vstar(:,1)
      vstar(:,-1) = -vstar(:,2)
    endif

    if (bc_ymax == no_slip) then
      u(:,ny+1) = -u(:,ny)
      u(:,ny+2) = -u(:,ny-1)
      v(:,ny+1) = -v(:,ny)
      v(:,ny+2) = -v(:,ny-1)
      ustar(:,ny+1) = -ustar(:,ny)
      ustar(:,ny+2) = -ustar(:,ny-1)
      vstar(:,ny+1) = -vstar(:,ny)
      vstar(:,ny+2) = -vstar(:,ny-1)
    endif

  end subroutine velocity_bcs

  subroutine velocity_face_bcs 

    ! uha is x faced
    ! vha is yfaced
    ! i.e. depending on the centered-ness relative to a 
    ! given boundary, different relationships must be used

    ! Periodic

    if (bc_xmin == periodic) then
      uha(  -2,:) = uha(nx-2,:)
      uha(  -1,:) = uha(nx-1,:)
      vha(   0,:) = vha(nx,:)
      vha(  -1,:) = vha(nx+1,:)
      uhx(  -2,:) = uhx(nx-2,:)
      uhx(  -1,:) = uhx(nx-1,:)
      vhx(  -2,:) = vhx(nx-2,:)
      vhx(  -1,:) = vhx(nx-1,:)
      uhy(   0,:) = uhy(nx,:)
      uhy(  -1,:) = uhy(nx+1,:)
      vhy(   0,:) = vhy(nx,:)
      vhy(  -1,:) = vhy(nx+1,:)
    endif 

    if (bc_xmax == periodic) then
      uha(nx+1,:) = uha(1,:)
      uha(nx+2,:) = uha(2,:)
      vha(nx+1,:) = vha(1,:)
      vha(nx+2,:) = vha(2,:)
      uhx(nx+1,:) = uhx(1,:)
      uhx(nx+2,:) = uhx(2,:)
      vhx(nx+1,:) = vhx(1,:)
      vhx(nx+2,:) = vhx(2,:)
      uhy(nx+1,:) = uhy(1,:)
      uhy(nx+2,:) = uhy(2,:)
      vhy(nx+1,:) = vhy(1,:)
      vhy(nx+2,:) = vhy(2,:)
    endif 

    if (bc_ymin == periodic) then
      uha(:,  -1) = uha(:,nx-1)
      uha(:,   0) = uha(:,nx)
      vha(:,  -1) = vha(:,nx-1)
      vha(:,  -2) = vha(:,nx-2)
      uhx(:,  -1) = uhx(:,nx-1)
      uhx(:,   0) = uhx(:,nx)
      vhx(:,  -1) = vhx(:,nx-1)
      vhx(:,   0) = vhx(:,nx)
      uhy(:,  -1) = uhy(:,nx-1)
      uhy(:,  -2) = uhy(:,nx-2)
      vhy(:,  -1) = vhy(:,nx-1)
      vhy(:,  -2) = vhy(:,nx-2)
    endif

    if (bc_ymax == periodic) then
      uha(:,ny+1) = uha(:,1)
      uha(:,ny+2) = uha(:,2)
      vha(:,ny+1) = vha(:,1)
      vha(:,ny+2) = vha(:,2)
      uhx(:,ny+1) = uhx(:,1)
      uhx(:,ny+2) = uhx(:,2)
      vhx(:,ny+1) = vhx(:,1)
      vhx(:,ny+2) = vhx(:,2)
      uhy(:,ny+1) = uhy(:,1)
      uhy(:,ny+2) = uhy(:,2)
      vhy(:,ny+1) = vhy(:,1)
      vhy(:,ny+2) = vhy(:,2)
    endif

    ! Zero gradient 

    if (bc_xmin == zero_gradient) then
      uha(  -2,:) = uha(2,:)
      uha(  -1,:) = uha(1,:)
      vha(   0,:) = vha(1,:)
      vha(  -1,:) = vha(2,:)
      uhx(  -2,:) = uhx(2,:)
      uhx(  -1,:) = uhx(1,:)
      vhx(  -2,:) = vhx(2,:)
      vhx(  -1,:) = vhx(1,:)
      uhy(   0,:) = uhy(1,:)
      uhy(  -1,:) = uhy(2,:)
      vhy(   0,:) = vhy(1,:)
      vhy(  -1,:) = vhy(2,:)
    endif 

    if (bc_xmax == zero_gradient) then
      uha(nx+1,:) = uha(nx-1,:)
      uha(nx+2,:) = uha(nx-2,:)
      vha(nx+1,:) = vha(nx,:)
      vha(nx+2,:) = vha(nx-1,:)
      uhx(nx+1,:) = uhx(nx-1,:)
      uhx(nx+2,:) = uhx(nx-2,:)
      vhx(nx+1,:) = vhx(nx-1,:)
      vhx(nx+2,:) = vhx(nx-2,:)
      uhy(nx+1,:) = uhy(nx,:)
      uhy(nx+2,:) = uhy(nx-1,:)
      vhy(nx+1,:) = vhy(nx,:)
      vhy(nx+2,:) = vhy(nx-1,:)
    endif 

    if (bc_ymin == zero_gradient) then
      uha(:,  -1) = uha(:,2)
      uha(:,   0) = uha(:,1)
      vha(:,  -1) = vha(:,1)
      vha(:,  -2) = vha(:,2)
      uhx(:,  -1) = uhx(:,2)
      uhx(:,   0) = uhx(:,1)
      vhx(:,  -1) = vhx(:,2)
      vhx(:,   0) = vhx(:,1)
      uhy(:,  -1) = uhy(:,1)
      uhy(:,  -2) = uhy(:,2)
      vhy(:,  -1) = vhy(:,1)
      vhy(:,  -2) = vhy(:,2)
    endif

    if (bc_ymax == zero_gradient) then
      uha(:,ny+1) = uha(:,ny)
      uha(:,ny+2) = uha(:,ny-1)
      vha(:,ny+1) = vha(:,ny-1)
      vha(:,ny+2) = vha(:,ny-2)
      uhx(:,ny+1) = uhx(:,ny)
      uhx(:,ny+2) = uhx(:,ny-1)
      vhx(:,ny+1) = vhx(:,ny)
      vhx(:,ny+2) = vhx(:,ny-1)
      uhy(:,ny+1) = uhy(:,ny-1)
      uhy(:,ny+2) = uhy(:,ny-2)
      vhy(:,ny+1) = vhy(:,ny-1)
      vhy(:,ny+2) = vhy(:,ny-2)
    endif

    ! No slip

    if (bc_xmin == no_slip) then
      uha(  -2,:) = -uha(2,:)
      uha(  -1,:) = -uha(1,:)
      vha(   0,:) = -vha(1,:)
      vha(  -1,:) = -vha(2,:)
      uhx(  -2,:) = -uhx(2,:)
      uhx(  -1,:) = -uhx(1,:)
      vhx(  -2,:) = -vhx(2,:)
      vhx(  -1,:) = -vhx(1,:)
      uhy(   0,:) = -uhy(1,:)
      uhy(  -1,:) = -uhy(2,:)
      vhy(   0,:) = -vhy(1,:)
      vhy(  -1,:) = -vhy(2,:)
    endif 

    if (bc_xmax == no_slip) then
      uha(nx+1,:) = -uha(nx-1,:)
      uha(nx+2,:) = -uha(nx-2,:)
      vha(nx+1,:) = -vha(nx,:)
      vha(nx+2,:) = -vha(nx-1,:)
      uhx(nx+1,:) = -uhx(nx-1,:)
      uhx(nx+2,:) = -uhx(nx-2,:)
      vhx(nx+1,:) = -vhx(nx-1,:)
      vhx(nx+2,:) = -vhx(nx-2,:)
      uhy(nx+1,:) = -uhy(nx,:)
      uhy(nx+2,:) = -uhy(nx-1,:)
      vhy(nx+1,:) = -vhy(nx,:)
      vhy(nx+2,:) = -vhy(nx-1,:)
    endif 

    if (bc_ymin == no_slip) then
      uha(:,  -1) = -uha(:,2)
      uha(:,   0) = -uha(:,1)
      vha(:,  -1) = -vha(:,1)
      vha(:,  -2) = -vha(:,2)
      uhx(:,  -1) = -uhx(:,2)
      uhx(:,   0) = -uhx(:,1)
      vhx(:,  -1) = -vhx(:,2)
      vhx(:,   0) = -vhx(:,1)
      uhy(:,  -1) = -uhy(:,1)
      uhy(:,  -2) = -uhy(:,2)
      vhy(:,  -1) = -vhy(:,1)
      vhy(:,  -2) = -vhy(:,2)
    endif

    if (bc_ymax == no_slip) then
      uha(:,ny+1) = -uha(:,ny)
      uha(:,ny+2) = -uha(:,ny-1)
      vha(:,ny+1) = -vha(:,ny-1)
      vha(:,ny+2) = -vha(:,ny-2)
      uhx(:,ny+1) = -uhx(:,ny)
      uhx(:,ny+2) = -uhx(:,ny-1)
      vhx(:,ny+1) = -vhx(:,ny)
      vhx(:,ny+2) = -vhx(:,ny-1)
      uhy(:,ny+1) = -uhy(:,ny-1)
      uhy(:,ny+2) = -uhy(:,ny-2)
      vhy(:,ny+1) = -vhy(:,ny-1)
      vhy(:,ny+2) = -vhy(:,ny-2)
    endif

  end subroutine velocity_face_bcs 

  subroutine phi_bcs

    !Periodic 

    if (bc_xmin == periodic) then
      phi(0,:) = phi(nx,:)
      phi(-1,:) = phi(nx-1,:)
    endif
    if (bc_xmax == periodic) then
      phi(nx+1,:) = phi(1,:)
      phi(nx+2,:) = phi(2,:)
    endif
    if (bc_ymin == periodic) then
      phi(:,0) = phi(:,ny)
      phi(:,-1) = phi(:,ny-1)
    endif
    if (bc_ymax == periodic) then
      phi(:,ny+1) = phi(:,1)
      phi(:,ny+2) = phi(:,2)
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

    ! No slip (no actually sure what is appropriate for phi at this stage?)

!!!    if (bc_xmin == no_slip) then
!!!      phi(0,:) = -phi(1,:)
!!!      phi(-1,:) = -phi(2,:)
!!!    endif
!!!    if (bc_xmax == no_slip) then
!!!      phi(nx+1,:) = -phi(nx,:)
!!!      phi(nx+2,:) = -phi(nx-1,:)
!!!    endif
!!!    if (bc_ymin == no_slip) then
!!!      phi(:,0) = -phi(:,1)
!!!      phi(:,-1) = -phi(:,2)
!!!    endif
!!!    if (bc_ymax == no_slip) then
!!!      phi(:,ny+1) = -phi(:,ny)
!!!      phi(:,ny+2) = -phi(:,ny-1)
!!!    endif

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

  end subroutine phi_bcs

  subroutine gradp_bcs

    ! Periodic 

    if (bc_xmin == periodic) then
      gradp_x(0,:) = gradp_x(nx,:)
      gradp_x(-1,:) = gradp_x(nx-1,:)
      gradp_y(0,:) = gradp_y(nx,:)
      gradp_y(-1,:) = gradp_y(nx-1,:)
    endif

    if (bc_xmax == periodic) then
      gradp_x(nx+1,:) = gradp_x(1,:)
      gradp_x(nx+2,:) = gradp_x(2,:)
      gradp_y(nx+1,:) = gradp_y(1,:)
      gradp_y(nx+2,:) = gradp_y(2,:)
    endif

    if (bc_ymin == periodic) then
      gradp_x(:,0) = gradp_x(:,ny)
      gradp_x(:,-1) = gradp_x(:,ny-1)
      gradp_y(:,0) = gradp_y(:,ny)
      gradp_y(:,-1) = gradp_y(:,ny-1)
    endif

    if (bc_ymin == periodic) then
      gradp_x(:,ny+1) = gradp_x(:,1)
      gradp_x(:,ny+2) = gradp_x(:,2)
      gradp_y(:,ny+1) = gradp_y(:,1)
      gradp_y(:,ny+2) = gradp_y(:,2)
    endif

    ! Zero Gradient
    if (bc_xmin == zero_gradient) then
      gradp_x(0,:) = gradp_x(1,:)
      gradp_x(-1,:) = gradp_x(2,:)
      gradp_y(0,:) = gradp_y(1,:)
      gradp_y(-1,:) = gradp_y(2,:)
    endif

    if (bc_xmax == zero_gradient) then
      gradp_x(nx+1,:) = gradp_x(nx,:)
      gradp_x(nx+2,:) = gradp_x(nx-1,:)
      gradp_y(nx+1,:) = gradp_y(nx,:)
      gradp_y(nx+2,:) = gradp_y(nx-1,:)
    endif

    if (bc_ymin == zero_gradient) then
      gradp_x(:,0) = gradp_x(:,1)
      gradp_x(:,-1) = gradp_x(:,2)
      gradp_y(:,0) = gradp_y(:,1)
      gradp_y(:,-1) = gradp_y(:,2)
    endif

    if (bc_ymax == zero_gradient) then
      gradp_x(:,ny+1) = gradp_x(:,ny)
      gradp_x(:,ny+2) = gradp_x(:,ny-1)
      gradp_y(:,ny+1) = gradp_y(:,ny)
      gradp_y(:,ny+2) = gradp_y(:,ny-1)
    endif

    ! No slip - similarly unsure of gradp conditions at this stage

!!!    if (bc_xmin == no_slip) then
!!!      gradp_x(0,:) = -gradp_x(1,:)
!!!      gradp_x(-1,:) = -gradp_x(2,:)
!!!      gradp_y(0,:) = -gradp_y(1,:)
!!!      gradp_y(-1,:) = -gradp_y(2,:)
!!!    endif
!!!
!!!    if (bc_xmax == no_slip) then
!!!      gradp_x(nx+1,:) = -gradp_x(nx,:)
!!!      gradp_x(nx+2,:) = -gradp_x(nx-1,:)
!!!      gradp_y(nx+1,:) = -gradp_y(nx,:)
!!!      gradp_y(nx+2,:) = -gradp_y(nx-1,:)
!!!    endif
!!!
!!!    if (bc_ymin == no_slip) then
!!!      gradp_x(:,0) = -gradp_x(:,1)
!!!      gradp_x(:,-1) = -gradp_x(:,2)
!!!      gradp_y(:,0) = -gradp_y(:,1)
!!!      gradp_y(:,-1) = -gradp_y(:,2)
!!!    endif
!!!
!!!    if (bc_ymin == no_slip) then
!!!      gradp_x(:,ny+1) = -gradp_x(:,ny)
!!!      gradp_x(:,ny+2) = -gradp_x(:,ny-1)
!!!      gradp_y(:,ny+1) = -gradp_y(:,ny)
!!!      gradp_y(:,ny+2) = -gradp_y(:,ny-1)
!!!    endif

    if (bc_xmin == no_slip) then
      gradp_x(0,:) = gradp_x(1,:)
      gradp_x(-1,:) = gradp_x(2,:)
      gradp_y(0,:) = gradp_y(1,:)
      gradp_y(-1,:) = gradp_y(2,:)
    endif

    if (bc_xmax == no_slip) then
      gradp_x(nx+1,:) = gradp_x(nx,:)
      gradp_x(nx+2,:) = gradp_x(nx-1,:)
      gradp_y(nx+1,:) = gradp_y(nx,:)
      gradp_y(nx+2,:) = gradp_y(nx-1,:)
    endif

    if (bc_ymin == no_slip) then
      gradp_x(:,0) = gradp_x(:,1)
      gradp_x(:,-1) = gradp_x(:,2)
      gradp_y(:,0) = gradp_y(:,1)
      gradp_y(:,-1) = gradp_y(:,2)
    endif

    if (bc_ymax == no_slip) then
      gradp_x(:,ny+1) = gradp_x(:,ny)
      gradp_x(:,ny+2) = gradp_x(:,ny-1)
      gradp_y(:,ny+1) = gradp_y(:,ny)
      gradp_y(:,ny+2) = gradp_y(:,ny-1)
    endif

  end subroutine gradp_bcs



end module boundary_conditions

