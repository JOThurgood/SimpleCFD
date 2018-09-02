! ers_euler_ideal_1d.f90
!
! Exact Riemann Solver for ideal 1D euler equations, producing a 
! converged exact solution.
! 
! <link to wiki> 

module shared_data ! a common block essentially

  implicit none

  integer, parameter :: num=selected_real_kind(p=15) 

  real(num) :: rhol, ul, pl, rhor, ur, pr !W_l and W_r
  real(num) :: gamma, gm, gp, g1, g2, g3, g4, g5, g6, g7 !gamma + const
  real(num) :: cl, cr !sound speeds
  real(num) :: ps !star region vars

  public 

  contains 

  subroutine constants
    !gamma 
    gm = gamma - 1.0_num
    gp = gamma + 1.0_num 
    g1 = gm / 2.0_num / gamma
    g2 = gp / 2.0_num / gamma
    g3 = 2.0_num * gamma / gm
    g4 = 2.0_num / gm
    g5 = 2.0_num / gp
    g6 = gm / gp
    g7 = gm / 2.0_num
    !sound speeds
    cl = sqrt(gamma * pl /rhol)
    cr = sqrt(gamma * pr / rhor)
  end subroutine constants

  !predefined initial conditions also kept here 

  subroutine test_1 !Sod's test
    rhol = 1.0_num
    rhor = 0.125_num !1.0_num / 8.0_num 
    ul = 0.0_num
    ur = 0.0_num
    pl = 1.0_num
    pr = 0.10_num 
  end subroutine test_1 

end module shared_data

module user

  use shared_data 

  implicit none 

  public 

  contains 

  subroutine initial_conditions

    ! set these for your problem
    rhol = 1.0_num
    rhor = 1.0_num
    ul = 0.0_num
    ur = 0.0_num
    pl = 1.0_num
    pr = 1.0_num
    gamma = 1.4_num

    ! call one of these to overwrite with one of the test cases
     call test_1 ! sods problem
    ! call test_2 ! 1,2,3 problem 
    ! call test_3 ! left-half of Woodward + Colella ()'s blast-wave 
    ! call test_4 ! right-half of ^ 
    ! call test_5 ! full blast(?) 

  end subroutine initial_conditions
  
  subroutine control
    !stub for setting sampling and output options
  end subroutine 

end module user 

module riemann !subroutines related to calculating star states

  use shared_data

  implicit none

  private 

  public :: check_positivity, pstar 

  !if when all is said and done there are any vars only used at 
  ! this stage it would be sensible to move here from shared_data

  contains 


  subroutine check_positivity
    !stub
  end subroutine check_positivity

  subroutine pstar
    call guess_ps 

    print *, ps
  end subroutine pstar

  subroutine guess_ps !initial guess on pstar 
    ! in the interests of speed you, can base your guess on the value of
    ! pstar as the  given by an appropriate choice of approximate
    ! R solver. This should mean that the Newton-Raphson iterates in as
    ! few steps as possible. 
    ps = 0.5_num * (pl + pr) ! For now just always use arithmetic mean
  end subroutine guess_ps
 
end module riemann

program ers_euler_ideal_1d

  use shared_data
  use user 
  use riemann 
  implicit none   

  call initial_conditions
!  call control 
  call constants 
! call check_positivity

  call pstar

end program ers_euler_ideal_1d

