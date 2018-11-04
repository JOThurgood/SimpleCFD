module initial_conditions

  use shared_data
  
  implicit none

  private 

  public :: set_ic

  contains

  subroutine set_ic ! NB: nx, t_end etc is set in control.f90

    !call shear_problem
    call minion_convergence_test

  end subroutine set_ic

  subroutine shear_problem

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

    shear_test = .true.
  end subroutine shear_problem


  subroutine minion_convergence_test
    real(num) :: x,y

    do iy = -1,ny+1
    do ix = -1,nx+1
      x = xc(ix)
      y = yc(iy)
      u(ix,iy) = 1.0_num - 2.0_num * cos(2.0_num*pi*x)*sin(2.0_num*pi*y)
      v(ix,iy) = 1.0_num + 2.0_num * sin(2.0_num*pi*x)*cos(2.0_num*pi*y)
    enddo
    enddo

    minion_test = .true.

  end subroutine minion_convergence_test


  ! Helper subroutines may go here! 

end module initial_conditions

