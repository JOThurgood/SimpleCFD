module advance_module

  use shared_data
  use boundary_conditions

  implicit none

  private 

  public :: advance_dt 

  contains

  subroutine advance_dt 


    call set_dt 
    print *,'Step', step,'dt', dt

    call velocity_bcs

    call step_1 ! Calculate time-centered normal velocities on the interfaces

    call step_2 ! Step 2 : MAC Projection 

    call step_3 ! Step 3: Reconstruct interface states consistent with MAC-projected
    ! velocities 

    call step_4 ! Step 4: Provisional update for full dt (star state)

    ! Step 5: Project provisional field to constraint
 
  end subroutine advance_dt

  subroutine step_1

    real(num) :: du, dv
    real(num) :: gp
    real(num) :: transv

    ! 1 - Calculate the advective velocities

    ! 1A -  calculate  time-centered interface states for the normal 
    ! velocities

    ! be careful with indicies (conversion between i+1/2 type notation 
    ! and separate indicies for xc and xb values)

    ! x faced data 

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
        gp = -0.0_num 
        uxl(ix,iy) = uhxl(ix,iy)  + transv + gp
        ! right
        transv = -0.5_num * dt * 0.5_num *(vha(ix+1,iy-1)+vha(ix+1,iy))&
          & * (uhy(ix+1,iy)-uhy(ix+1,iy-1 )) /dy
        gp = -0.0_num 
        uxr(ix,iy) = uhxr(ix,iy)  + transv + gp
        ! also calc the tangential vel states for step 3
        ! left 
        transv = -0.5_num * dt * 0.5_num * (vha(ix,iy-1) + vha(ix,iy)) &
          & * (vhy(ix,iy)-vhy(ix,iy-1)) / dy
        gp = -0.0_num 
        vxl(ix,iy) = vhxl(ix,iy) + transv + gp
        ! right 
        transv = -0.5_num * dt * 0.5_num *(vha(ix+1,iy-1)+vha(ix+1,iy))&
          & * (vhy(ix+1,iy)-vhy(ix+1,iy-1 )) /dy
        gp = -0.0_num 
        vxr(ix,iy) = vhxr(ix,iy) + transv + gp 
      endif
      if (ix /= 0) then !can do the yface stuff
        ! normal components
        transv = -0.5_num * dt * 0.5_num * (uha(ix-1,iy) + uha(ix,iy)) &
          & * (vhx(ix,iy) - vhx(ix-1,iy)) / dx
        gp = -0.0_num
        vyl(ix,iy) = vhyl(ix,iy) + transv + gp

        transv = -0.5_num * dt * 0.5_num *(uha(ix-1,iy+1)+uha(ix,iy+1))&
          & * (vhx(ix,iy+1) - vhx(ix-1,iy+1)) / dx
        gp = -0.0_num
        vyr(ix,iy) = vhyr(ix,iy) + transv + gp

        ! also calc the tangential vel states for step 3
        transv = -0.5_num * dt * 0.5_num * (uha(ix-1,iy) + uha(ix,iy)) &
          & * (uhx(ix,iy) - uhx(ix-1,iy)) / dx
        gp = -0.0_num
        uyl(ix,iy) = uhyl(ix,iy) + transv + gp

        transv = -0.5_num * dt * 0.5_num *(uha(ix-1,iy+1)+uha(ix,iy+1))&
          & * (uhx(ix,iy+1) - uhx(ix-1,iy+1)) / dx
        gp = -0.0_num
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


  end subroutine step_1
  
  subroutine step_2

    real(num) :: residual, L_phi
    real(num) :: tol = 5e-3_num
    real(num) :: correction
    integer :: relaxation_steps = 0
    integer :: max_relaxation_steps = 100000

    ! MAC Projection via Gauss-Seidel

    print *, 'Step2'


    ! calc divU at cc using the MAC velocities
    do iy = 1, ny
    do ix = 1, nx
      divu(ix,iy) = (macu(ix,iy) - macu(ix-1,iy) ) /dx &
        & + (macv(ix,iy) - macv(ix,iy-1))/dy
    enddo
    enddo

    print *, 'max divu before cleaning',maxval(abs(divu))


    do 

      relaxation_steps = relaxation_steps + 1

      call phi_bcs 

      do iy = 1, ny
      do ix = 1, nx
        phi(ix,iy) = 0.25_num * ( & 
          & phi(ix+1,iy) + phi(ix-1,iy) + phi(ix,iy+1) + phi(ix,iy-1) &
          - dx*dy * divu(ix,iy) ) 
      enddo
      enddo

      ! reset the maxmimum val of residual
      residual = -1e6_num 
   
      do iy = 1, ny  
      do ix = 1, nx  
        L_phi = (phi(ix+1,iy) - 2.0_num*phi(ix,iy) + phi(ix-1,iy)) &
          & / dx**2 + & 
          & (phi(ix,iy+1) - 2.0_num*phi(ix,iy) + phi(ix,iy-1)) / dy**2 
   
        residual = max(residual, abs(divu(ix,iy) - L_phi))
!        print *,abs(divu(ix,iy)-L_phi)  
    end do
      end do 

      if (residual <= tol) exit
      if (relaxation_steps > max_relaxation_steps) then
        print *, 'failed to converge quickly, STOP'
        STOP 
      endif
    enddo

    print *, 'Step', step,'relaxation finished in',relaxation_steps, &
      & 'residual',residual

    ! perform the divergence cleaning

    call phi_bcs

    do iy = 0, ny
    do ix = 0, nx
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

    print *, 'max divu after cleaning',maxval(abs(divu))

  end subroutine step_2

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

  end subroutine step_3

  subroutine step_4

    real(num) :: Au, Av, gp ! evaluation of advection term

    do iy = 1, ny
    do ix = 1, nx
      Au = 0.5_num * (macu(ix-1,iy)+macu(ix,iy)) * (ux(ix,iy)-ux(ix-1,iy))/dx &
        &+ 0.5_num * (macv(ix,iy-1)+macv(ix,iy)) * (uy(ix,iy)-uy(ix,iy-1))/dy 
      Av = 0.5_num * (macu(ix-1,iy)+macu(ix,iy)) * (vx(ix,iy)-vx(ix-1,iy))/dx &
        &+ 0.5_num * (macv(ix,iy-1)+macv(ix,iy)) * (vy(ix,iy)-vy(ix,iy-1))/dy 
      gp = 0.0_num 
      ustar(ix,iy) = u(ix,iy) - dt * Au - dt * gp
      vstar(ix,iy) = v(ix,iy) - dt * Av - dt * gp
    enddo
    enddo

  end subroutine step_4

  subroutine step_5 ! Approximate projection 

  end subroutine step_5


  subroutine set_dt

    real(num) :: dtx, dty 

    dtx = CFL * dx / maxval(abs(u))
    dty = CFL * dy / maxval(abs(v))
    dt = MIN(dtx,dty)

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

