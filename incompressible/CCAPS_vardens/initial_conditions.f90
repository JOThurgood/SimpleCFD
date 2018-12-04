module initial_conditions

  use shared_data
  
  implicit none

  contains

  subroutine set_ic ! NB: nx, t_end etc is set in control.f90

    u = 0.0_num
    v = 0.0_num

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

    rho = 1.0_num

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

    rho = 1.0_num

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

  subroutine rti1_ic

    real(num) :: x, y 
    real(num) :: d, rho0, rho1
    real(num) :: fwtm , amp

    amp = 0.1_num

    d = x_max - x_min
    fwtm = 0.1_num *d
!    fwtm = 0.01_num * d
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

    grav_x = 0.0_num
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

