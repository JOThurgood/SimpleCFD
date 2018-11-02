module advance_module

  use shared_data
  use boundary_conditions
  use diagnostics

  implicit none

  private 

  public :: advance_dt 

  contains

  subroutine advance_dt 


    call set_dt 

    call step_1 ! Calculate time-centered normal velocities on the interfaces

    call step_2 ! Step 2 : MAC Projection 

    call step_3 ! Step 3: Reconstruct interface states consistent with MAC-projected
    ! velocities 

    call step_4 ! Step 4: Provisional update for full dt (star state)

    call step_5 ! Step 5: Project provisional field to constraint
 
    time = time + dt

  end subroutine advance_dt

  subroutine step_1

    real(num) :: gp
    real(num) :: du, dv
    real(num) :: transv

    ! 1 - Calculate the advective velocities

    ! 1A -  calculate  time-centered interface states for the normal 
    ! velocities

    ! be careful with indicies (conversion between i+1/2 type notation 
    ! and separate indicies for xc and xb values)

    ! x faced data 

    call velocity_bcs

    do iy = 1, ny 
    do ix = 0, nx  !xb counts from 0 to nx, <0 and >nx are ghosts 
  
      du = minmod((u(ix,iy)-u(ix-1,iy))/dx,(u(ix+1,iy)-u(ix,iy))/dx)
      uhxl(ix,iy) = u(ix,iy) + &
        0.5_num * (1.0_num - dt * u(ix,iy) / dx ) * du * dx

      du = minmod((u(ix+1,iy)-u(ix,iy))/dx,(u(ix+2,iy)-u(ix+1,iy))/dx)
      uhxr(ix,iy) = u(ix+1,iy) - &
        0.5_num * (1.0_num + dt * u(ix+1,iy) / dx ) * du * dx

      dv = minmod((v(ix,iy)-v(ix-1,iy))/dx, (v(ix+1,iy)-v(ix,iy))/dx) 
      vhxl(ix,iy) = v(ix,iy) + & 
        & 0.5_num * (1.0_num - dt * u(ix,iy) / dx ) * dv * dx

      dv = minmod((v(ix+1,iy)-v(ix,iy))/dx, (v(ix+2,iy)-v(ix+1,iy))/dx) 
      vhxr(ix,iy) = v(ix+1,iy) - & 
        & 0.5_num * (1.0_num + dt * u(ix+1,iy) / dx) * dv * dx 
      
    enddo 
    enddo

    ! y faced data

    do iy = 0, ny 
    do ix = 1, nx

      du = minmod( (u(ix,iy)-u(ix,iy-1))/dy, (u(ix,iy+1)-u(ix,iy))/dy)
      uhyl(ix,iy) = u(ix,iy) + &
        & 0.5_num * (1.0_num - dt * u(ix,iy) / dy) * du * dy

      du = minmod( (u(ix,iy+1)-u(ix,iy))/dy, (u(ix,iy+2)-u(ix,iy+1))/dy)
      uhyr(ix,iy) = u(ix,iy+1) - &
        & 0.5_num * (1.0_num + dt * u(ix,iy+1) / dy) * du * dy

      dv = minmod( (v(ix,iy)-v(ix,iy-1))/dy, (v(ix,iy+1)-v(ix,iy))/dy)
      vhyl(ix,iy) = v(ix,iy) + &
        & 0.5_num * (1.0_num - dt * v(ix,iy) / dy ) * dv * dy

      dv = minmod( (v(ix,iy+1)-v(ix,iy))/dy, (v(ix,iy+2)-v(ix,iy+1))/dy)
      vhyr(ix,iy) = v(ix,iy+1) - &
        & 0.5_num * (1.0_num + dt * v(ix,iy+1) / dy) * dv * dy

    enddo
    enddo

    ! 1B - Use the normal velocities to calculate the advective vel
    ! by solving a Riemann problem

    do iy = 0, ny
    do ix = 0, ny
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
    do ix = 0, ny
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

    ! if you find bugs in 1D - check not also in 3D which is similar

    do iy = 0, ny
    do ix = 0, ny
      if (iy /= 0) then !can do the xface stuff
        ! normal components
        ! left
        transv = -0.5_num * dt * 0.5_num * (vha(ix,iy-1) + vha(ix,iy)) &
          & * (uhy(ix,iy)-uhy(ix,iy-1)) / dy
        gp = -0.5_num * dt * gradp_x(ix,iy) 
        uxl(ix,iy) = uhxl(ix,iy)  + transv + gp
        ! right
        transv = -0.5_num * dt * 0.5_num *(vha(ix+1,iy-1)+vha(ix+1,iy))&
          & * (uhy(ix+1,iy)-uhy(ix+1,iy-1 )) /dy
        gp = -0.5_num * dt * gradp_x(ix+1,iy) 
        uxr(ix,iy) = uhxr(ix,iy)  + transv + gp
        ! also calc the tangential vel states for step 3
        ! left 
        transv = -0.5_num * dt * 0.5_num * (vha(ix,iy-1) + vha(ix,iy)) &
          & * (vhy(ix,iy)-vhy(ix,iy-1)) / dy
        gp = -0.5_num * dt * gradp_y(ix,iy) 
        vxl(ix,iy) = vhxl(ix,iy) + transv + gp
        ! right 
        transv = -0.5_num * dt * 0.5_num *(vha(ix+1,iy-1)+vha(ix+1,iy))&
          & * (vhy(ix+1,iy)-vhy(ix+1,iy-1 )) /dy
        gp = -0.5_num * dt * gradp_y(ix+1,iy) 
        vxr(ix,iy) = vhxr(ix,iy) + transv + gp 
      endif
      if (ix /= 0) then !can do the yface stuff
        ! normal components
        transv = -0.5_num * dt * 0.5_num * (uha(ix-1,iy) + uha(ix,iy)) &
          & * (vhx(ix,iy) - vhx(ix-1,iy)) / dx
        gp = -0.5_num * dt * gradp_y(ix,iy) 
        vyl(ix,iy) = vhyl(ix,iy) + transv + gp

        transv = -0.5_num * dt * 0.5_num *(uha(ix-1,iy+1)+uha(ix,iy+1))&
          & * (vhx(ix,iy+1) - vhx(ix-1,iy+1)) / dx
        gp = -0.5_num * dt * gradp_y(ix,iy+1) 
        vyr(ix,iy) = vhyr(ix,iy) + transv + gp

        ! also calc the tangential vel states for step 3
        transv = -0.5_num * dt * 0.5_num * (uha(ix-1,iy) + uha(ix,iy)) &
          & * (uhx(ix,iy) - uhx(ix-1,iy)) / dx
        gp = -0.5_num * dt * gradp_x(ix,iy)
        uyl(ix,iy) = uhyl(ix,iy) + transv + gp

        transv = -0.5_num * dt * 0.5_num *(uha(ix-1,iy+1)+uha(ix,iy+1))&
          & * (uhx(ix,iy+1) - uhx(ix-1,iy+1)) / dx
        gp = -0.5_num * dt * gradp_x(ix,iy+1)
        uyr(ix,iy) = uhyr(ix,iy) + transv + gp

     endif
    enddo
    enddo

    ! 1E Final riemann solve + upwinding for full normal velocities 
    ! (sometimes AKA the MAC velocities)

    ! if you find bugs here check 3E also

    do iy = 0, ny
    do ix = 0, ny
      if (iy /= 0) then !can do the xface stuff
        ua(ix,iy) = riemann(uxl(ix,iy),uxr(ix,iy))
      endif
      if (ix /= 0) then !can do the yface stuff
        va(ix,iy) = riemann(vyl(ix,iy), vyr(ix,iy)) 
      endif
    enddo
    enddo

    do iy = 0, ny
    do ix = 0, ny
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

!    call plot_divergence_now
!    if (step /=0) call plot_divergence_now

    call step2_gauss_seidel

    print *, '*** max divu before cleaning',maxval(abs(divu))

    do ix = 0, nx
    do iy = 0, ny
      if (iy /= 0) then !can do the xface stuff
        correction = (phi(ix+1,iy) - phi(ix,iy))/dx
        macu(ix,iy) = macu(ix,iy) - correction 
      endif
      if (ix /= 0) then !can do the yface stuff
        correction = (phi(ix,iy+1)-phi(ix,iy))/dy
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

!    if (step /=0) call plot_divergence_now
    !call plot_divergence_now
  end subroutine step_2

  subroutine step2_gauss_seidel

    logical :: gsrb=.true. !should move to control eventually

    real(num) :: tol = 1e-18_num
    real(num) :: L2, L2_old !norms
    real(num) :: L_phi !lagrangian of phi
    integer :: maxir = 100000000
    integer :: ir = 0 
    logical :: verbose=.false. 

 
    print *, '*** begining relaxation to solve for phi.'
    print *, '*** this can take a while, use VERBOSE if you want to monitor stepping'

    phi = 0.0_num
    L2_old = 1e6_num

    do
      ir = ir + 1 
   
!      do iy = 1, ny  
!      do ix = 1, nx  
!        !is this the appropriate form of GS for these equations!!!!
!        phi(ix,iy) = 0.25_num * ( & 
!          & phi(ix+1,iy) + phi(ix-1,iy) + phi(ix,iy+1) + phi(ix,iy-1) &
!          - dx**2 * divu(ix,iy) ) 
!   
!      end do
!      end do 

      ! if not using redblack order
      if (.not. gsrb) then 
        do iy = 1, ny  
        do ix = 1, nx  
          phi(ix,iy) = 0.25_num * ( & 
            & phi(ix+1,iy) + phi(ix-1,iy) + phi(ix,iy+1) + phi(ix,iy-1) &
            - dx**2 * divu(ix,iy) ) 
        end do
        end do 
      else !use red black
        ! odd iteration
        do iy = 1, ny  
        do ix = 1, nx  
          if (modulo(ix+iy,2) == 1) then
            phi(ix,iy) = 0.25_num * ( & 
              & phi(ix+1,iy) + phi(ix-1,iy) + phi(ix,iy+1) + phi(ix,iy-1) &
              - dx**2 * divu(ix,iy) ) 
          endif
        end do
        end do 
!        call phi_bcs
        ! even iteration
        do iy = 1, ny  
        do ix = 1, nx  
          if (modulo(ix+iy,2) == 0) then
            phi(ix,iy) = 0.25_num * ( & 
              & phi(ix+1,iy) + phi(ix-1,iy) + phi(ix,iy+1) + phi(ix,iy-1) &
              - dx**2 * divu(ix,iy) ) 
          endif
        end do
        end do 
   
      endif
   
      ! Apply periodic boundary conditions on phi's ghost cells

      call phi_bcs  
 
      L2 = 0.0_num
      do iy = 1, ny  
      do ix = 1, nx  
        L_phi = (phi(ix+1,iy) - 2.0_num*phi(ix,iy) + phi(ix-1,iy)) &
          & / dx**2 + & 
          & (phi(ix,iy+1) - 2.0_num*phi(ix,iy) + phi(ix,iy-1)) / dy**2 
   
        L2 = L2 + abs(divu(ix,iy)-L_phi)**2
      end do
      end do 
      L2 = sqrt( L2 / REAL(nx*ny,num))
   
   
      if (verbose) &
        & print *, 'GS-Step2-iteration',ir,'complete. L2 is',L2,'|L2-L2_old| is',abs(L2-L2_old),'tol is',tol
   
     !exit conditions 
     !if ((step >= nsteps .and. nsteps >= 0) .or. (L2 <= tol)) exit
     ! alt, exit if the difference between L2 and L2 prev is small - might
     ! indicate convergence
      if ((ir >= maxir .and. maxir >= 0) .or. (abs(L2-L2_old) <= tol)) exit
      L2_old = L2
    end do



  print *, '*** Gauss Seidel relaxation completed in',ir,'steps'
  print *, '*** L2 norm on D(velocity_field)- L(phi)',L2
 

  end subroutine step2_gauss_seidel


  subroutine step_3

    ! reconstruct the interface states using the MAC velocities
    ! for consistency
    ! (redo some of step 1 but use mac velocities  for upwinding)

    ! Because you haven't been deallocating arrays etc, you can directly
    ! skip a few of the recalculations

    ! I think we only need to re-do D (extra terms) and E 

    ! Step 1D calculated the full tangential velocity states 
    ! in anticipation of this - only need E which differs 

    ! Step 3E Upwind face components based upon MAC vels
    

    do iy = 0, ny
    do ix = 0, ny
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

!!!    ! If you uncommment the following, it suggests that step_3 does nothing for
!!!    ! the test IC...
!!!    ! might just be that you're divergence cleaning adjustment was consistent 
!!!    ! with the initial upwinding.... Not 100% happy with this yet
!!! 
!!!
!!!   ! check the equivalent to the macu and macv (ux and vy) still has
!!!   ! small divergence?
!!!
!!!     do iy = 1, ny
!!!     do ix = 1, nx
!!!       divu(ix,iy) = (ux(ix,iy) - ux(ix-1,iy) ) /dx &
!!!         & + (vy(ix,iy) - vy(ix,iy-1))/dy
!!!     enddo
!!!     enddo
!!!  
!!!     print *, 'max divu after 3E',maxval(abs(divu))
!!!  

    print *, 'Step #3 completed normally'
  end subroutine step_3

  subroutine step_4

    real(num) :: Au, Av ! evaluation of advection term

    do iy = 1, ny
    do ix = 1, nx
      Au = 0.5_num * (macu(ix-1,iy)+macu(ix,iy)) * (ux(ix,iy)-ux(ix-1,iy))/dx &
        &+ 0.5_num * (macv(ix,iy-1)+macv(ix,iy)) * (uy(ix,iy)-uy(ix,iy-1))/dy 
      Av = 0.5_num * (macu(ix-1,iy)+macu(ix,iy)) * (vx(ix,iy)-vx(ix-1,iy))/dx &
        &+ 0.5_num * (macv(ix,iy-1)+macv(ix,iy)) * (vy(ix,iy)-vy(ix,iy-1))/dy 
      ustar(ix,iy) = u(ix,iy) - dt * Au - dt * gradp_x(ix,iy)
      vstar(ix,iy) = v(ix,iy) - dt * Av - dt * gradp_y(ix,iy)
    enddo
    enddo

    print *, 'Step #4 completed normally'
  end subroutine step_4

  subroutine step_5

    real(num) :: correction, gpsi

    print *,'Step #5'

    ! calc divU at cc using the star velocities which themselves are cc
    ! (this differs to step two which uses face vars to get a CC var)

    call velocity_bcs

    do iy = 1, ny
    do ix = 1, nx
      divu(ix,iy) = (ustar(ix+1,iy) - ustar(ix-1,iy))/dx/2.0_num &
        & + (vstar(ix,iy+1) - vstar(ix,iy-1))/dy/2.0_num
    enddo
    enddo

!    if (step /=0) call plot_divergence_now

    divu = divu/dt

    call step5_gauss_seidel

    print *, '*** max divu before cleaning',maxval(abs(divu)*dt)
!    print *, '*** max divu/dt before cleaning',maxval(abs(divu))

    do iy = 1, ny
    do ix = 1, nx

      gpsi = (phi(ix+1,iy) - phi(ix-1,iy))/dx/2.0_num
      correction = dt * gpsi
      u(ix,iy) = ustar(ix,iy) - correction 

      gpsi = (phi(ix,iy+1)-phi(ix,iy-1))/dy/2.0_num
      correction = dt * gpsi
      v(ix,iy) = vstar(ix,iy) - correction 

    enddo
    enddo

    ! calculate the divergence of the updated velocity field
    call velocity_bcs

    do iy = 1, ny
    do ix = 1, nx
      divu(ix,iy) = (u(ix+1,iy) - u(ix-1,iy))/dx/2.0_num &
        & + (v(ix,iy+1) - v(ix,iy-1))/dy/2.0_num
    enddo
    enddo

!    if (step /=0) call plot_divergence_now

    print *, '*** max divu after cleaning',maxval(abs(divu))
   !   print *, '*** max divu/dt after cleaning',maxval(abs(divu/dt))


    ! update the pressure gradient 
  
    do ix = 1, nx
    do iy = 1, ny

      gpsi = (phi(ix+1,iy) - phi(ix-1,iy))/dx/2.0_num
      gradp_x(ix,iy) = gradp_x(ix,iy) + gpsi

      gpsi = (phi(ix,iy+1)-phi(ix,iy-1))/dy/2.0_num
      gradp_y(ix,iy) = gradp_y(ix,iy) + gpsi
 
    enddo
    enddo

    call gradp_bcs

  end subroutine step_5

  subroutine step5_gauss_seidel

    logical :: gsrb=.true. !should move to control eventually

    real(num) :: tol = 1e-16_num
    real(num) :: L2, L2_old !norms
    real(num) :: L_phi !lagrangian of phi
    integer :: maxir = 10000000
    integer :: ir = 0 
    logical :: verbose=.false. 

 
    print *, '*** begining relaxation to solve for phi.'
    print *, '*** this can take a while, use VERBOSE if you want to monitor stepping'

    call phi_bcs !use old phi as initial guess
    L2_old = 1e6_num

    do
      ir = ir + 1 
   
!      do iy = 1, ny  
!      do ix = 1, nx  
!        !is this the appropriate form of GS for these equations!!!!
!        phi(ix,iy) = 0.25_num * ( & 
!          & phi(ix+1,iy) + phi(ix-1,iy) + phi(ix,iy+1) + phi(ix,iy-1) &
!          - dx**2 * divu(ix,iy) ) 
!   
!      end do
!      end do 


      ! if not using redblack order
      if (.not. gsrb) then 
        do iy = 1, ny  
        do ix = 1, nx  
          phi(ix,iy) = 0.25_num * ( & 
            & phi(ix+1,iy) + phi(ix-1,iy) + phi(ix,iy+1) + phi(ix,iy-1) &
            - dx**2 * divu(ix,iy) ) 
        end do
        end do 
     
      else !use red black
   
        ! odd iteration
        do iy = 1, ny  
        do ix = 1, nx  
          if (modulo(ix+iy,2) == 1) then
            phi(ix,iy) = 0.25_num * ( & 
              & phi(ix+1,iy) + phi(ix-1,iy) + phi(ix,iy+1) + phi(ix,iy-1) &
              - dx**2 * divu(ix,iy) ) 
          endif
        end do
        end do 
   
        ! even iteration
   
        do iy = 1, ny  
        do ix = 1, nx  
          if (modulo(ix+iy,2) == 0) then
            phi(ix,iy) = 0.25_num * ( & 
              & phi(ix+1,iy) + phi(ix-1,iy) + phi(ix,iy+1) + phi(ix,iy-1) &
              - dx**2 * divu(ix,iy) ) 
          endif
        end do
        end do 
   
      endif

      ! Apply periodic boundary conditions on phi's ghost cells

      call phi_bcs  
 
      L2 = 0.0_num
      do iy = 1, ny  
      do ix = 1, nx  
        L_phi = (phi(ix+1,iy) - 2.0_num*phi(ix,iy) + phi(ix-1,iy)) &
          & / dx**2 + & 
          & (phi(ix,iy+1) - 2.0_num*phi(ix,iy) + phi(ix,iy-1)) / dy**2 
   
        L2 = L2 + abs(divu(ix,iy)-L_phi)**2
      end do
      end do 
      L2 = sqrt( L2 / REAL(nx*ny,num))
   
   
      if (verbose) &
        & print *, 'GS-Step5-iteration',ir,'complete. L2 is',L2,'|L2-L2_old| is',abs(L2-L2_old),'tol is',tol
   
     !exit conditions 
     !if ((step >= nsteps .and. nsteps >= 0) .or. (L2 <= tol)) exit
     ! alt, exit if the difference between L2 and L2 prev is small - might
     ! indicate convergence
      if ((ir >= maxir .and. maxir >= 0) .or. (abs(L2-L2_old) <= tol)) exit
      L2_old = L2
    end do



    print *, '*** Gauss Seidel relaxation completed in',ir,'steps'
    print *, '*** L2 norm on D(velocity_field)- L(phi)',L2
 

  end subroutine step5_gauss_seidel

  subroutine set_dt

    real(num) :: dtx, dty 

    dtx = CFL * dx / maxval(abs(u))
    dty = CFL * dy / maxval(abs(v))
    dt = MIN(dtx,dty)

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

end module advance_module

