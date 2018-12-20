module advance_module

  use shared_data
  use boundary_conditions
  use diagnostics
  use gauss_seidel
  use multigrid

  implicit none

  private 

  public :: advance_dt 

  type(mg_input) :: input

  contains

  subroutine advance_dt 

    call set_dt 

    call step_1 ! Calculate time-centered normal velocities on the interfaces

    call step_2 ! Step 2 : MAC Projection 

    call step_3 ! Step 3: Reconstruct interface states consistent with constraint

    if (use_vardens) call advect_dens ! consv. update of density based on mac vels

    call step_4 ! Step 4: Provisional update for full dt (star state)

    call step_5 ! Step 5: Project provisional field to constraint
 
    time = time + dt

  end subroutine advance_dt

  subroutine step_1

    real(num) :: du, dv
    real(num) :: transv
    real(num) :: LU_l, LU_r ! vector laplacian of U left/right - for viscous forcing 

    ! 1 - Calculate the advective velocities

    ! 1A -  calculate  time-centered interface states for the normal 
    ! velocities

    ! be careful with indicies (conversion between i+1/2 type notation 
    ! and separate indicies for xc and xb values)

    ! x faced data 

    call velocity_bcs(arr_cc = u,di = 1)
    call velocity_bcs(arr_cc = v,di = 2)

    do iy = 1, ny 
    do ix = 0, nx  !xb counts from 0 to nx, <0 and >nx are ghosts 
  
      if (use_minmod) then
        du = minmod((u(ix,iy)-u(ix-1,iy))/dx,(u(ix+1,iy)-u(ix,iy))/dx)
      else
        du = ( u(ix+1,iy) - u(ix-1,iy) ) /2.0_num / dx
      endif

      uhxl(ix,iy) = u(ix,iy) + &
        0.5_num * (1.0_num - dt * u(ix,iy) / dx ) * du * dx

      if (use_minmod) then
        du = minmod((u(ix+1,iy)-u(ix,iy))/dx,(u(ix+2,iy)-u(ix+1,iy))/dx)
      else
        du= ( u(ix+2,iy)-u(ix,iy) ) / 2.0_num / dx
      endif

      uhxr(ix,iy) = u(ix+1,iy) - &
        0.5_num * (1.0_num + dt * u(ix+1,iy) / dx ) * du * dx

      if (use_minmod) then
        dv = minmod((v(ix,iy)-v(ix-1,iy))/dx, (v(ix+1,iy)-v(ix,iy))/dx) 
      else
        dv = (v(ix+1,iy) - v(ix-1,iy))/2.0_num/dx
      endif
      vhxl(ix,iy) = v(ix,iy) + & 
        & 0.5_num * (1.0_num - dt * u(ix,iy) / dx ) * dv * dx

      if (use_minmod) then
        dv = minmod((v(ix+1,iy)-v(ix,iy))/dx, (v(ix+2,iy)-v(ix+1,iy))/dx) 
      else
        dv = (v(ix+2,iy)-v(ix,iy)) /2.0_num/dx
      endif
      vhxr(ix,iy) = v(ix+1,iy) - & 
        & 0.5_num * (1.0_num + dt * u(ix+1,iy) / dx) * dv * dx 
      
    enddo 
    enddo

    ! y faced data

    do iy = 0, ny 
    do ix = 1, nx

      if (use_minmod) then
        du = minmod( (u(ix,iy)-u(ix,iy-1))/dy, (u(ix,iy+1)-u(ix,iy))/dy)
      else
        du = (u(ix,iy+1) - u(ix,iy-1)) / 2.0_num / dy
      endif

      uhyl(ix,iy) = u(ix,iy) + &
        & 0.5_num * (1.0_num - dt * v(ix,iy) / dy) * du * dy

      if (use_minmod) then
        du = minmod( (u(ix,iy+1)-u(ix,iy))/dy, (u(ix,iy+2)-u(ix,iy+1))/dy)
      else
        du = (u(ix,iy+2)-u(ix,iy)) / 2.0_num/dy
      endif

      uhyr(ix,iy) = u(ix,iy+1) - &
        & 0.5_num * (1.0_num + dt * v(ix,iy+1) / dy) * du * dy

      if (use_minmod) then
        dv = minmod( (v(ix,iy)-v(ix,iy-1))/dy, (v(ix,iy+1)-v(ix,iy))/dy)
      else
        dv = (v(ix,iy+1)-v(ix,iy-1))/2.0_num/dy
      endif
      vhyl(ix,iy) = v(ix,iy) + &
        & 0.5_num * (1.0_num - dt * v(ix,iy) / dy ) * dv * dy

      if (use_minmod) then
        dv = minmod( (v(ix,iy+1)-v(ix,iy))/dy, (v(ix,iy+2)-v(ix,iy+1))/dy)
      else
        dv = (v(ix,iy+2)-v(ix,iy))/2.0_num/dy
      endif
      vhyr(ix,iy) = v(ix,iy+1) - &
        & 0.5_num * (1.0_num + dt * v(ix,iy+1) / dy) * dv * dy

    enddo
    enddo

    ! 1B - Use the normal velocities to calculate the advective vel
    ! by solving a Riemann problem

    do iy = 0, ny
    do ix = 0, nx
      if (iy /= 0) then !can do the xface stuff
        uha(ix,iy) = riemann(uhxl(ix,iy),uhxr(ix,iy))
      endif
      if (ix /= 0) then !can do the yface stuff
        vha(ix,iy) = riemann(vhyl(ix,iy), vhyr(ix,iy)) 
      endif
    enddo
    enddo

    ! 1C - Upwind the hat states (normal velocity predictions)
    ! using the advective vels 

    do iy = 0, ny
    do ix = 0, nx
      if (iy /= 0) then !can do the xface stuff
        uhx(ix,iy) = upwind(uha(ix,iy),uhxl(ix,iy),uhxr(ix,iy))
        vhx(ix,iy) = upwind(uha(ix,iy),vhxl(ix,iy),vhxr(ix,iy)) 
      endif
      if (ix /= 0) then !can do the yface stuff
        uhy(ix,iy) = upwind(vha(ix,iy),uhyl(ix,iy),uhyr(ix,iy))
        vhy(ix,iy) = upwind(vha(ix,iy),vhyl(ix,iy),vhyr(ix,iy)) 
      endif
    enddo
    enddo

    ! 1D construct the full left and right predictions of normal
    ! velocities on the interfaces

    ! (actually get them on all interfaces, as needed later steps)

    call velocity_bcs(arr_xface = uha, di = 1)
    call velocity_bcs(arr_yface = vha, di = 2)
    call velocity_bcs(arr_xface = uhx, arr_yface = uhy, di = 1)
    call velocity_bcs(arr_xface = vhx, arr_yface = vhy, di = 2)

    do iy = 0, ny
    do ix = 0, nx
      if (iy /= 0) then !can do the xface stuff
        ! normal components

        ! left
        transv = -0.5_num * dt * 0.5_num * (vha(ix,iy-1) + vha(ix,iy)) &
          & * (uhy(ix,iy)-uhy(ix,iy-1)) / dy
        uxl(ix,iy) = uhxl(ix,iy)  + transv + 0.5_num * dt * get_force_cc(ix,iy,1)

        ! right
        transv = -0.5_num * dt * 0.5_num *(vha(ix+1,iy-1)+vha(ix+1,iy))&
          & * (uhy(ix+1,iy)-uhy(ix+1,iy-1 )) /dy
        uxr(ix,iy) = uhxr(ix,iy)  + transv + 0.5_num * dt * get_force_cc(ix+1,iy,1)


        ! also calc the tangential vel states for step 3

        ! left 
        transv = -0.5_num * dt * 0.5_num * (vha(ix,iy-1) + vha(ix,iy)) &
          & * (vhy(ix,iy)-vhy(ix,iy-1)) / dy
        vxl(ix,iy) = vhxl(ix,iy) + transv + 0.5_num * dt * get_force_cc(ix,iy,2) 

        ! right 
        transv = -0.5_num * dt * 0.5_num *(vha(ix+1,iy-1)+vha(ix+1,iy))&
          & * (vhy(ix+1,iy)-vhy(ix+1,iy-1 )) /dy
        vxr(ix,iy) = vhxr(ix,iy) + transv + 0.5_num * dt * get_force_cc(ix+1,iy,2) 
      endif
      if (ix /= 0) then !can do the yface stuff
        ! normal components
        transv = -0.5_num * dt * 0.5_num * (uha(ix-1,iy) + uha(ix,iy)) &
          & * (vhx(ix,iy) - vhx(ix-1,iy)) / dx
        vyl(ix,iy) = vhyl(ix,iy) + transv + 0.5_num * dt * get_force_cc(ix,iy,2)

        transv = -0.5_num * dt * 0.5_num *(uha(ix-1,iy+1)+uha(ix,iy+1))&
          & * (vhx(ix,iy+1) - vhx(ix-1,iy+1)) / dx
        vyr(ix,iy) = vhyr(ix,iy) + transv + 0.5_num * dt * get_force_cc(ix,iy+1,2) 

        ! also calc the tangential vel states for step 3
        transv = -0.5_num * dt * 0.5_num * (uha(ix-1,iy) + uha(ix,iy)) &
          & * (uhx(ix,iy) - uhx(ix-1,iy)) / dx
        uyl(ix,iy) = uhyl(ix,iy) + transv + 0.5_num * dt * get_force_cc(ix,iy,1)

        transv = -0.5_num * dt * 0.5_num *(uha(ix-1,iy+1)+uha(ix,iy+1))&
          & * (uhx(ix,iy+1) - uhx(ix-1,iy+1)) / dx
        uyr(ix,iy) = uhyr(ix,iy) + transv + 0.5_num * dt * get_force_cc(ix,iy+1,1)

     endif
    enddo
    enddo

    ! also include viscous forcing if using viscosity

    if (use_viscosity) then

      if (use_vardens) then 
        print *,'warning: using variable density and viscosity. Dont forget to add 1/rho terms to visc forcing'
        print *,'stopping in step 1D'
        STOP
      endif

      do ix = 0, nx
      do iy = 0, ny
        if (iy /=0 ) then ! x faces

          ! vector laplacian in cartesians L(U) = (L u_x, L u_y, L u_z)
          ! x component, in cell centers to left and right of interface

          LU_l = get_LU_cc(ix,iy,1)
          LU_r = get_LU_cc(ix+1,iy,1) 
  
          uxl(ix,iy) = uxl(ix,iy) + 0.5_num * dt * visc * LU_l 
          uxr(ix,iy) = uxr(ix,iy) + 0.5_num * dt * visc * LU_r 
  
          ! vector laplacian of y component
  
          LU_l = get_LU_cc(ix,iy,2)
          LU_r = get_LU_cc(ix+1,iy,2)


          vxl(ix,iy) = vxl(ix,iy) + 0.5_num * dt * visc * LU_l
          vxr(ix,iy) = vxr(ix,iy) + 0.5_num * dt * visc * LU_r
  
  
        endif
        if (ix /=0) then ! yfaces

          !_l and _r is now vector laplacian at cc "below" and "above" y face
  
          LU_l = get_LU_cc(ix,iy,1)
          LU_r = get_LU_cc(ix,iy+1,1) 

          uyl(ix,iy) = uyl(ix,iy) + 0.5_num * dt * visc * LU_l
          uyr(ix,iy) = uyr(ix,iy) + 0.5_num * dt * visc * LU_r

          LU_l = get_LU_cc(ix,iy,2)
          LU_r = get_LU_cc(ix,iy+1,2)

          vyl(ix,iy) = vyl(ix,iy) + 0.5_num * dt * visc * LU_l
          vyr(ix,iy) = vyr(ix,iy) + 0.5_num * dt * visc * LU_r

        endif
      enddo
      enddo
    endif


    ! 1E Final riemann solve + upwinding for full normal velocities 
    ! (sometimes AKA the MAC velocities)

    ! if you find bugs here check 3E also

    do iy = 0, ny
    do ix = 0, nx
      if (iy /= 0) then !can do the xface stuff
        ua(ix,iy) = riemann(uxl(ix,iy),uxr(ix,iy))
      endif
      if (ix /= 0) then !can do the yface stuff
        va(ix,iy) = riemann(vyl(ix,iy), vyr(ix,iy)) 
      endif
    enddo
    enddo

    do iy = 0, ny
    do ix = 0, nx
      if (iy /= 0) then !can do the xface stuff
        macu(ix,iy) = upwind(ua(ix,iy),uxl(ix,iy),uxr(ix,iy))
      endif
      if (ix /= 0) then !can do the yface stuff
        macv(ix,iy) = upwind(va(ix,iy),vyl(ix,iy),vyr(ix,iy)) 
      endif
    enddo
    enddo

    print *, 'Step #1 completed normally'

  end subroutine step_1

  subroutine step_2

    real(num) :: correction

    print *, 'Step #2'
    print *, '*** start'

    ! calc divU at cc using the MAC velocities
    do iy = 1, ny
    do ix = 1, nx
      divu(ix,iy) = (macu(ix,iy) - macu(ix-1,iy) ) /dx &
        & + (macv(ix,iy) - macv(ix,iy-1))/dy
    enddo
    enddo
 
!call plot_divergence_now ! debug
!if (step /=0) call plot_divergence_now ! debug

    if (use_vardens) call rho_bcs ! needed for any OOB in relax and correction 

    if (use_vardens) then
      call solve_variable_elliptic(phigs = phi, f = divu(1:nx,1:ny), &
        & eta= 1.0_num / rho(0:nx+1,0:ny+1), &
        & use_old_phi = .false., tol = 1e-18_num) 
    else
!      call solve_const_Helmholtz(phigs = phi, f = divu(1:nx,1:ny), &
!        & alpha = 0.0_num, beta = -1.0_num, &
!        & use_old_phi = .false., tol = 1e-18_num) 

      input = mg_input(tol = 1e-16_num, nx = nx, ny = ny, dx = dx, dy = dy, &
            & f = divu(1:nx,1:ny), phi=phi, &
            & bc_xmin = mg_bc_xmin, bc_xmax = mg_bc_xmax, &
            & bc_ymin = mg_bc_ymin, bc_ymax = mg_bc_ymax ) 

      call mg_interface(input)
      phi = input%phi

print *,'warning hardcoded perodic for MG, needs general handling'



    endif

    print *, '*** max divu before cleaning',maxval(abs(divu))

    do ix = 0, nx
    do iy = 0, ny
      if (iy /= 0) then !can do the xface stuff
        correction = (phi(ix+1,iy) - phi(ix,iy))/dx
        if (use_vardens) correction = correction / &
            (0.5_num * (rho(ix,iy) + rho(ix+1,iy)))
        macu(ix,iy) = macu(ix,iy) - correction 
      endif
      if (ix /= 0) then !can do the yface stuff
        correction = (phi(ix,iy+1)-phi(ix,iy))/dy
        if (use_vardens) correction = correction / &
            (0.5_num * (rho(ix,iy) + rho(ix,iy+1)))
        macv(ix,iy) = macv(ix,iy) - correction 
      endif
    enddo
    enddo

    ! calculate the new divergence
    ! calc divU at cc using the MAC velocities

    do iy = 1, ny
    do ix = 1, nx
      divu(ix,iy) = (macu(ix,iy) - macu(ix-1,iy) ) /dx &
        & + (macv(ix,iy) - macv(ix,iy-1))/dy
    enddo
    enddo


    print *, '*** max divu after cleaning',maxval(abs(divu))
    print *, '*** complete'

!if (step /=0) call plot_divergence_now ! debug
!call plot_divergence_now ! debug
  end subroutine step_2

  subroutine step_3

    ! reconstruct the interface states using the MAC velocities
    ! for consistency
    ! (redo some of step 1 but use mac velocities  for upwinding)

    ! Because you haven't been overwriting or deallocating 
    ! such arrays, most of it doesnt have to be recalculated

    ! We only need to re-do E

    ! Step 3E Upwind face components based upon MAC vels
    

    do iy = 0, ny
    do ix = 0, nx
      if (iy /= 0) then !can do the xface stuff
        ux(ix,iy) = upwind(macu(ix,iy),uxl(ix,iy),uxr(ix,iy))
        vx(ix,iy) = upwind(macu(ix,iy),vxl(ix,iy),vxr(ix,iy))
      endif
      if (ix /= 0) then !can do the yface stuff
        uy(ix,iy) = upwind(macv(ix,iy),uyl(ix,iy),uyr(ix,iy))
        vy(ix,iy) = upwind(macv(ix,iy),vyl(ix,iy),vyr(ix,iy)) 
      endif
    enddo
    enddo

    print *, 'Step #3 completed normally'
  end subroutine step_3


  subroutine step_4
    real(num) :: Au, Av ! evaluation of advection term

    real(num),dimension(1:nx,1:ny) :: f !viscous only

    ! do inviscid calculation first since the terms
    ! used to evaluate ustar are used in the "f" term
    do iy = 1, ny
    do ix = 1, nx
      Au = get_Au(ix,iy)
      Av = get_Av(ix,iy) 
      ustar(ix,iy) = u(ix,iy) - dt * Au + dt * get_force_cc(ix,iy,1)
      vstar(ix,iy) = v(ix,iy) - dt * Av + dt * get_force_cc(ix,iy,2)
    enddo
    enddo
  
    if (use_viscosity) then

      print *,'Step #4 explicit viscosity on '
      print *,'*** begin implicit solve for ustar'

      do ix = 1, nx
      do iy = 1, ny
        f(ix,iy) = ustar(ix,iy) + 0.5_num * dt * visc * get_LU_cc(ix,iy,1)
      enddo
      enddo    

      if (use_vardens) then
        !stub
        print *,'cannot use viscosity and variable density at same time currently'
        print *,'TERMINATING'
        STOP 
      else

        input = mg_input(tol = 1e-16_num, nx = nx, ny = ny, dx = dx, dy = dy, &
              & f = f, phi=ustar, &
!              & use_as_init_guess = .true., & 
              & bc_xmin = mg_bc_xmin, bc_xmax = mg_bc_xmax, &
              & bc_ymin = mg_bc_ymin, bc_ymax = mg_bc_ymax ,&
              & phi_bc_ymax = 1.0_num, & 
              & const_helmholtz = .true. , ch_alpha = 1.0_num, ch_beta = dt*visc/2.0_num) 
  
        call mg_interface(input)
        ustar(1:nx,1:ny) = input%phi(1:nx,1:ny)

!        call solve_const_Helmholtz(phigs = ustar, f = f, &
!          alpha = 1.0_num, beta = dt*visc/2.0_num, &
!          & use_old_phi = .true., tol = 1e-16_num) 
      endif

      print *,'*** begin implicit solve for vstar' 

      do ix = 1, nx
      do iy = 1, ny
        f(ix,iy) = vstar(ix,iy) + 0.5_num * dt * visc * get_LU_cc(ix,iy,2)
      enddo
      enddo    

      if (use_vardens) then
        !stub
      else
!        call solve_const_Helmholtz(phigs = vstar, f = f, &
!          alpha = 1.0_num, beta = dt*visc/2.0_num, &
!          & use_old_phi = .true., tol = 1e-16_num) 

        input = mg_input(tol = 1e-16_num, nx = nx, ny = ny, dx = dx, dy = dy, &
              & f = f, phi=vstar, &
 !             & use_as_init_guess = .true., & 
              & bc_xmin = mg_bc_xmin, bc_xmax = mg_bc_xmax, &
              & bc_ymin = mg_bc_ymin, bc_ymax = mg_bc_ymax ,&
              & const_helmholtz = .true. , ch_alpha = 1.0_num, ch_beta = dt*visc/2.0_num) 
  
        call mg_interface(input)
        vstar(1:nx,1:ny) = input%phi(1:nx,1:ny)
      endif

    endif

    print *, 'Step #4 completed normally'
  end subroutine step_4

  subroutine step_5

    real(num) :: correction, gpsi

    print *,'Step #5'
    print *, '*** start'

    ! calc divU at cc using the star velocities which themselves are cc
    ! (this differs to step two which uses face vars to get a CC var)

    call velocity_bcs(arr_cc = ustar, di = 1)
    call velocity_bcs(arr_cc = vstar, di = 2)

    do iy = 1, ny
    do ix = 1, nx
      divu(ix,iy) = (ustar(ix+1,iy) - ustar(ix-1,iy))/dx/2.0_num &
        & + (vstar(ix,iy+1) - vstar(ix,iy-1))/dy/2.0_num
    enddo
    enddo

!if (step /=0) call plot_divergence_now
!call plot_divergence_now


    divu = divu/dt

    if (use_vardens) then
      call solve_variable_elliptic(phigs = phi, f = divu(1:nx,1:ny), &
        & eta= 1.0_num / rho(0:nx+1,0:ny+1), &
        & use_old_phi = .true., tol = 1e-18_num) 
    else
!      call solve_const_Helmholtz(phigs = phi, f = divu(1:nx,1:ny), &
!        alpha = 0.0_num, beta = -1.0_num, &
!        & use_old_phi = .true., tol = 1e-16_num) 

      input = mg_input(tol = 1e-16_num, nx = nx, ny = ny, dx = dx, dy = dy, &
            & f = divu(1:nx,1:ny), phi=phi, &
            & bc_xmin = mg_bc_xmin, bc_xmax = mg_bc_xmax, &
            & bc_ymin = mg_bc_ymin, bc_ymax = mg_bc_ymax ) 

      call mg_interface(input)
      phi = input%phi

print *,'warning hardcoded perodic for MG, needs general handling'
    endif

    print *, '*** max divu before cleaning',maxval(abs(divu)*dt)


    call phi_bcs

    do iy = 1, ny
    do ix = 1, nx

      gpsi = (phi(ix+1,iy) - phi(ix-1,iy))/dx/2.0_num
      correction = dt * gpsi
      if (use_vardens) correction = correction/rho(ix,iy) !cc here
      u(ix,iy) = ustar(ix,iy) - correction 

      gpsi = (phi(ix,iy+1)-phi(ix,iy-1))/dy/2.0_num
      correction = dt * gpsi
      if (use_vardens) correction = correction/rho(ix,iy)
      v(ix,iy) = vstar(ix,iy) - correction 

    enddo
    enddo

    ! calculate the divergence of the updated velocity field

    call velocity_bcs(arr_cc = u, di = 1)
    call velocity_bcs(arr_cc = v, di = 2)

    do iy = 1, ny
    do ix = 1, nx
      divu(ix,iy) = (u(ix+1,iy) - u(ix-1,iy))/dx/2.0_num &
        & + (v(ix,iy+1) - v(ix,iy-1))/dy/2.0_num
    enddo
    enddo

    print *, '*** max divu after cleaning',maxval(abs(divu))

!if (step /=0) call plot_divergence_now ! debug
!call plot_divergence_now ! debug
!call plot_vel_now ! debug

    ! update the pressure gradient 
  
    do ix = 0, nx+1 
    do iy = 0, ny+1 

      gpsi = (phi(ix+1,iy) - phi(ix-1,iy))/dx/2.0_num
      gradp_x(ix,iy) = gradp_x(ix,iy) + gpsi

      gpsi = (phi(ix,iy+1)-phi(ix,iy-1))/dy/2.0_num
      gradp_y(ix,iy) = gradp_y(ix,iy) + gpsi
 
    enddo
    enddo

!call plot_gradp_now ! debug

  end subroutine step_5

  subroutine advect_dens

    real(num) :: drho

    call rho_bcs(arr_cc = rho)

    ! calculate rhohat states on faces 

    do ix = 0, nx
    do iy = 0, ny
      if (iy /= 0) then ! xface vars
      
        if (use_minmod) then
          drho = minmod((rho(ix,iy)-rho(ix-1,iy))/dx,(rho(ix+1,iy)-rho(ix,iy))/dx)
        else
          drho = ( rho(ix+1,iy) - rho(ix-1,iy) ) /2.0_num / dx
        endif

        rhohxl(ix,iy) = rho(ix,iy) + &
          & 0.5_num * (1.0_num - dt * macu(ix,iy) / dx ) * drho * dx 

        if (use_minmod) then
          drho = minmod((rho(ix+1,iy)-rho(ix,iy))/dx,(rho(ix+2,iy)-rho(ix+1,iy))/dx)
        else
          drho= ( rho(ix+2,iy)-rho(ix,iy) ) / 2.0_num / dx
        endif

        rhohxr(ix,iy) = rho(ix+1,iy) - &
          & 0.5_num * (1.0_num + dt * macu(ix,iy) / dx) * drho * dx

      endif

      if (ix /= 0) then !y face vars

        if (use_minmod) then
          drho = minmod( (rho(ix,iy)-rho(ix,iy-1))/dy, (rho(ix,iy+1)-rho(ix,iy))/dy)
        else
          drho = (rho(ix,iy+1) - rho(ix,iy-1)) / 2.0_num / dy
        endif

        rhohyl(ix,iy) = rho(ix,iy) + &
          & 0.5_num * (1.0_num - dt * macv(ix,iy) / dy ) * drho * dy 
          
        if (use_minmod) then
          drho = minmod( (rho(ix,iy+1)-rho(ix,iy))/dy, (rho(ix,iy+2)-rho(ix,iy+1))/dy)
        else
          drho = (rho(ix,iy+2)-rho(ix,iy)) / 2.0_num/dy
        endif

        rhohyr(ix,iy) = rho(ix,iy+1) - &
          & 0.5_num * (1.0_num + dt * macv(ix,iy) / dy) * drho * dy

      endif
    enddo
    enddo

    ! upwind using the MAC velocities

    do ix = 0, nx
    do iy = 0, ny
      if (iy /= 0) then ! xface vars
        rhohx(ix,iy) = upwind(macu(ix,iy),rhohxl(ix,iy),rhohxr(ix,iy))
      endif
      if (ix /= 0) then !y face vars
        rhohy(ix,iy) = upwind(macv(ix,iy),rhohyl(ix,iy),rhohyr(ix,iy))
      endif
    enddo
    enddo

    ! calculate full states with transverse terms

    call rho_bcs(arr_xface = rhohx, arr_yface = rhohy) 
    call velocity_bcs(arr_xface = macu, di = 1)
    call velocity_bcs(arr_yface = macv, di = 2)

    do ix = 0, nx
    do iy = 0, ny
      if (iy /= 0) then ! xface vars
        rhoxl(ix,iy) = rhohxl(ix,iy) &
          & - 0.5_num * dt * rho(ix,iy) * (macu(ix,iy) - macu(ix-1,iy)) / dx &
          & - 0.5_num * dt / dy * ( rhohy(ix,iy)*macv(ix,iy) - &
                                      & rhohy(ix,iy-1)*macv(ix,iy-1) )
        rhoxr(ix,iy) = rhohxr(ix,iy) &
          & - 0.5_num * dt * rho(ix+1,iy) * (macu(ix+1,iy) - macu(ix,iy)) / dx &
          & - 0.5_num * dt / dy * ( rhohy(ix+1,iy)*macv(ix+1,iy) - &
                                    & rhohy(ix+1,iy-1)*macv(ix+1,iy-1) )
      endif
      if (ix /= 0) then ! yface vars
         rhoyl(ix,iy) = rhohyl(ix,iy) &
          & - 0.5_num * dt * rho(ix,iy) * (macv(ix,iy) - macv(ix,iy-1)) / dy &
          & - 0.5_num * dt / dx * ( rhohx(ix,iy)*macu(ix,iy) - &
                                      & rhohx(ix-1,iy)*macu(ix-1,iy) )
        rhoyr(ix,iy) = rhohyr(ix,iy) &
          & - 0.5_num * dt * rho(ix,iy+1) * (macv(ix,iy+1) - macv(ix,iy)) / dy &
          & - 0.5_num * dt / dx * ( rhohx(ix,iy+1)*macu(ix,iy+1) - &
                                    & rhohx(ix-1,iy+1)*macu(ix-1,iy+1) )
      endif
    enddo
    enddo

    ! resolve states via upwind

    do ix = 0, nx
    do iy = 0, ny
      if (iy /= 0) then ! xface vars
        rhox(ix,iy) = upwind(macu(ix,iy),rhoxl(ix,iy),rhoxr(ix,iy))
      endif
      if (ix /= 0) then !y face vars
        rhoy(ix,iy) = upwind(macv(ix,iy),rhoyl(ix,iy),rhoyr(ix,iy))
      endif
    enddo
    enddo

    ! simple conservative update to new time level

    do ix = 1, nx
    do iy = 1, ny
      rho(ix,iy) = rho(ix,iy) &
        & - dt / dx * (rhox(ix,iy)*macu(ix,iy) - rhox(ix-1,iy)*macu(ix-1,iy)) &
        & - dt / dy * (rhoy(ix,iy)*macv(ix,iy) - rhoy(ix,iy-1)*macv(ix,iy-1))
    enddo
    enddo

  end subroutine advect_dens
  
  subroutine set_dt

    real(num) :: dtx, dty
    real(num) :: dtf

    ! need to call bcs to capture velocities on driven boundaries
    call velocity_bcs(arr_cc = u, di = 1)  
    call velocity_bcs(arr_cc = v, di = 2)

    dtx = CFL * dx / maxval(abs(u))
    dty = CFL * dy / maxval(abs(v))
    dt = MIN(dtx,dty)

    if (sqrt(grav_x**2 + grav_y**2) > 1e-16_num) then
      if (use_vardens) then
        dtf = CFL * sqrt(2.0_num * dx / maxval(abs(gradp_x-grav_x*rho)))
        dt = MIN(dt,dtf)
        dtf = CFL * sqrt(2.0_num * dy / maxval(abs(gradp_y-grav_y*rho)))
        dt = MIN(dt,dtf)
      else
        dtf = CFL * sqrt(2.0_num * dx / maxval(abs(gradp_x-grav_x)))
        dt = MIN(dt,dtf)
        dtf = CFL * sqrt(2.0_num * dy / maxval(abs(gradp_y-grav_y)))
        dt = MIN(dt,dtf)
      endif
    endif 

    print *, 'hydro dt = ',dt

  end subroutine set_dt

  real(num) function minmod(a,b)  
    real(num), intent(in) :: a, b 
    if ( (abs(a) < abs(b)) .and. (a*b > 0.0_num) ) then
      minmod = a        
    else if ( (abs(a) > abs(b)) .and. (a*b > 0.0_num) ) then
      minmod = b        
    else                
      minmod = 0.0_num
    endif               
  end function minmod

  real(num) function riemann(a,b)  
    real(num), intent(in) :: a, b 
    if ( (a > 0.0_num) .and. ( (a+b)>0.0_num) ) then
      riemann = a
    else if ( (a <= 0.0_num) .and. (b >= 0.0_num) ) then
      riemann = 0.0_num
    else 
      riemann = b
    endif 
  end function riemann

  real(num) function upwind(sadv, a, b)
    real(num), intent(in) :: sadv,a,b
    if (sadv > 0.0_num) then
      upwind = a
    else if (sadv < 0.0_num) then
      upwind = b
    else
      upwind = 0.5_num * (a + b)
    endif
  end function upwind

  real(num) function get_force_cc(ix,iy,di)
    integer,intent(in) :: ix,iy, di
    real(num) :: grav_tmp_x
    real(num) :: grav_tmp_y

    grav_tmp_x = grav_x
    grav_tmp_y = grav_y


! debug: uncomment to turn grav off at the closest cc's to the edges
!    grav_tmp_x = 0.0_num
!    grav_tmp_y = 0.0_num
!    if ( (ix > 1) .and. (ix < nx) ) grav_tmp_x = grav_x
!    if ( (iy > 1) .and. (iy < ny) ) grav_tmp_y = grav_y


    if (di==1) then
      if (use_vardens) then
        get_force_cc = -gradp_x(ix,iy)/rho(ix,iy) + grav_tmp_x
      else
        get_force_cc = -gradp_x(ix,iy) + grav_tmp_x
      endif
    else if (di==2) then
      if (use_vardens) then
        get_force_cc = -gradp_y(ix,iy)/rho(ix,iy) + grav_tmp_y
      else
        get_force_cc = -gradp_y(ix,iy) + grav_tmp_y
      endif
  
    else
      print *,'error get_force_cc given invalid dimension'
      print *,'di = 1(x) or =2(y)'
      print *,'STOP'
      STOP
    endif

  end function get_force_cc

  real(num) function get_LU_cc(ix,iy,di) ! calculate vector laplacian at a coordinate
    integer, intent(in) :: ix, iy, di

    if (di == 1) then

      get_LU_cc = (u(ix+1,iy) - 2.0_num*u(ix,iy) + u(ix-1,iy)) / dx**2 + & 
         (u(ix,iy+1) - 2.0_num*u(ix,iy) + u(ix,iy-1)) / dy**2 

    else if (di == 2) then

      get_LU_cc = (v(ix+1,iy) - 2.0_num*v(ix,iy) + v(ix-1,iy)) / dx**2 + & 
         (v(ix,iy+1) - 2.0_num*v(ix,iy) + v(ix,iy-1)) / dy**2 
    else 

      print *,'error: get_LU_cc (calc vector laplacian) not given valid dimension'
      print *,'di = 1 (x) or = 2 (y)'
      print *,'Terminating early'
      STOP

    endif

  endfunction get_LU_cc

  real(num) function get_Au(ix,iy) 

    integer, intent(in) :: ix,iy 
    get_Au = 0.5_num * (macu(ix-1,iy)+macu(ix,iy)) * (ux(ix,iy)-ux(ix-1,iy))/dx &
       &+ 0.5_num * (macv(ix,iy-1)+macv(ix,iy)) * (uy(ix,iy)-uy(ix,iy-1))/dy 

  endfunction get_Au

  real(num) function get_Av(ix,iy) 

    integer, intent(in) :: ix,iy 
    get_Av = 0.5_num * (macu(ix-1,iy)+macu(ix,iy)) * (vx(ix,iy)-vx(ix-1,iy))/dx &
     &+ 0.5_num * (macv(ix,iy-1)+macv(ix,iy)) * (vy(ix,iy)-vy(ix,iy-1))/dy 

  endfunction get_Av

end module advance_module

