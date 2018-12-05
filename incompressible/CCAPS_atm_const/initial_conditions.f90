module initial_conditions

  use shared_data
  
  implicit none

  contains

  subroutine set_ic ! NB: nx, t_end etc is set in control.f90

    real(num) :: rho_lo
    real(num) :: rho_hi
    real(num) :: p0_lo
    real(num) :: fwtm, d
    real(num) :: y

    ! flow 

    u = 0.0_num
    v = 0.0_num

    ! free to either set rho0, then get rho from it via
    !   rho(ix,iy = rho0(iy) + perturbation
    ! or set rho(ix,iy) then calc 
    !   rho0(iy) =  1/Nx Sum_x (rho(ix,iy))
    !
    ! (i.e. user choice, all must be done in initial_conditions though)

    ! Simple case of unperturbed heavy on top of light (i.e. rho0(iy) = rho(:,iy) )

    rho_lo = 1.0_num
    rho_hi = 7.0_num
    d = x_max - x_min
    fwtm = 0.1_num *d

    do iy = -1, ny+2
      y = yc(iy)-0.5_num
      y = y * 1.818_num / fwtm
      rho0(iy) = rho_lo + 0.5_num * (rho_hi-rho_lo) * (1.0_num + tanh(y))
      rho(-1:nx+2,iy) = rho0(iy)
    enddo

    ! Calculate a HSE background pressure p0 

    p0_lo = 10.0_num ! just sets a Gauge? 
        ! However, due to p0^(<1) terms, cant tolerate a negative number in the profile
        ! so needs to be sufficiently high to not NAN the betas 

    p0(-1) = p0_lo
    do iy = 0, ny+2
      p0(iy) = p0(iy-1) + 0.5_num * dy * (rho0(iy-1) + rho0(iy)) * grav_y
    enddo



  end subroutine set_ic

  subroutine shear_test_ic

    real(num) :: x,y, rho_s, delta_s

    rho_s = 42.0_num
    delta_s = 0.05_num

    do iy = -1,ny+1
    do ix = -1,nx+1
      x = xc(ix)
      y = yc(iy)
      if (y <= 0.5_num) then
        u(ix,iy) =  tanh(rho_s * (y-0.25_num))
      else
        u(ix,iy) = tanh(rho_s * (0.75_num -y))
      endif 
      v(ix,iy) = delta_s * sin(2.0_num * pi * x)
    enddo
    enddo

    ! for homogeneous tests, should just be able to set background to uniform with 
    ! an arbitary pressure 

    rho = 1.0_num
    rho0 = 1.0_num
    p0 = 1.0_num

  end subroutine shear_test_ic

  subroutine minion_test_ic

    real(num) :: x,y

    do iy = -1,ny+1
    do ix = -1,nx+1
      x = xc(ix)
      y = yc(iy)
      u(ix,iy) = 1.0_num - 2.0_num * cos(2.0_num*pi*x)*sin(2.0_num*pi*y)
      v(ix,iy) = 1.0_num + 2.0_num * sin(2.0_num*pi*x)*cos(2.0_num*pi*y)
    enddo
    enddo

    do iy = 1, ny
    do ix = 1, nx
      divu(ix,iy) = (u(ix+1,iy) - u(ix-1,iy))/dx/2.0_num &
        & + (v(ix,iy+1) - v(ix,iy-1))/dy/2.0_num
    enddo
    enddo

    ! for homogeneous tests, should just be able to set background to uniform with 
    ! an arbitary pressure 

    rho = 1.0_num
    rho0 = 1.0_num
    p0 = 1.0_num


  end subroutine minion_test_ic

  subroutine vardens_adv_test_ic

    real(num) :: x,y,r

    u = 0.0_num
    v = 1.0_num

    rho = 1.0_num
    do iy = -1,ny+1
    do ix = -1,nx+1
      x = xc(ix)-0.5_num
      y = yc(iy)-0.5_num
      r = sqrt(x**2 + y**2)
      if (r <= 0.1_num) rho(ix,iy) = 2.0_num
    enddo
    enddo
    print *, 'rho on grid',sum(rho(1:nx,1:ny)*dx*dy)

    print *,'max numerical divu in ICs',maxval(abs(divu))
  end subroutine vardens_adv_test_ic

! rho0 now clashes with array for background rho - rewrite
  subroutine rti1_ic

    real(num) :: x, y 
    real(num) :: d, rho0, rho1
    real(num) :: fwtm , amp

    amp = 0.1_num

    d = x_max - x_min
    fwtm = 0.1_num *d
    rho0 = 1.0_num
    rho1 = 7.0_num 

    u = 0.0_num
    v = 0.0_num 

    rho = rho0

    do ix = -1, nx+1
    do iy = -1, ny+2
      x = xc(ix)

!      if (step /= 0) then ! trust post step0 to have HSE built in . dont call on bootstrap?
!          y = yc(iy) !+ amp * d * (cos(2.0_num * pi * x / d))
!          y = y * 1.818_num / fwtm
          rho(ix,iy) = rho0 + 0.5_num * (rho1-rho0) * (1.0_num + tanh(y))
!          gradp_x(ix,iy) = rho(ix,iy) * grav_x
!          gradp_y(ix,iy) = rho(ix,iy) * grav_y
!      endif

      !now overwrite rho with the perturbed + desingularised profile
      y = yc(iy) + amp * d * (cos(2.0_num * pi * x / d))
      y = y * 1.818_num / fwtm
      rho(ix,iy) = rho0 + 0.5_num * (rho1-rho0) * (1.0_num + tanh(y))

    enddo
    enddo

  end subroutine rti1_ic

  subroutine blob1_ic

    real(num) :: x,y,r

    ! this should probabally be de-singularised with a tanh profile on rho=rho(r)

    u = 0.0_num
    v = 0.0_num

    grav_y = -1.0_num

    rho = 1.0_num

    do iy = -1,ny+1
    do ix = -1,nx+1
      x = xc(ix)-0.5_num
      y = yc(iy)-0.5_num
      r = sqrt(x**2 + y**2)
      if (r <= 0.1_num) rho(ix,iy) = 2.0_num
    enddo
    enddo

    print *, 'rho on grid',sum(rho(1:nx,1:ny)*dx*dy)

  end subroutine blob1_ic

end module initial_conditions

