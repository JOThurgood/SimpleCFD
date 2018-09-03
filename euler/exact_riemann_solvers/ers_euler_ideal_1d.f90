! ers_euler_ideal_1d.f90
!
! Exact Riemann Solver for ideal 1D euler equations, producing a 
! converged exact solution.
! 
! <link to wiki> 

module shared_data ! a common block essentially

  implicit none

  integer, parameter :: num=selected_real_kind(p=6) 
!  integer, parameter :: num=selected_real_kind(p=15) 

  real(num) :: rhol, ul, pl, rhor, ur, pr !W_l and W_r
  real(num) :: gamma, gm, gp, g1, g2, g3, g4, g5, g6, g7 !gamma + const
  real(num) :: al,ar, bl, br !more constants
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
    !data dependant constants for finding the root of the "f" function
    al = g5 / rhol
    ar = g5 / rhor
    bl = g6 * pl
    br = g6 * pr
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
    gamma = 1.4_num
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

    print *,'rhol,ul,pl',rhol,ul,pl
    print *,'rhor,ur,pr',rhor,ur,pr
    print *,'p0',ps
    print *,'f(p0)', f(ps)
    print *,'fprime(p0)', fprime(ps)
    call newton_raphson
    print *, ps
  end subroutine pstar

  subroutine guess_ps !initial guess on pstar 
    ! in the interests of speed you, can base your guess on the value of
    ! pstar as the  given by an appropriate choice of approximate
    ! R solver. This should mean that the Newton-Raphson iterates in as
    ! few steps as possible. 
    ps = 0.5_num * (pl + pr) ! For now just always use arithmetic mean
  end subroutine guess_ps
 
  subroutine newton_raphson
    integer :: i = 0
    real(num) :: pold
    real(num) :: tol = 1e-6_num
    real(num) :: rpc = 1e15_num ! init as huge for exit logix
    do
      if (rpc < tol) exit !condition on tolerance
      pold = ps 
      ps = ps - f(ps)/fprime(ps)
      rpc = 2.0_num * abs((ps - pold) / (ps + pold)) 
!      rpc = 2.0_num * abs(ps - pold) / abs(ps + pold)
      i = i + 1
      print *,'i',i,'pold',pold,'psnew',ps, 'rpc',rpc
    enddo 
    print *, 'Newton-Raphson converged in',i,'iterations'
  end subroutine newton_raphson

  real(num) function f(p) !root function
    real(num), intent(in) :: p
    real(num) :: ak,bk, pk, prat, ck
    integer :: i 

    f = 0.0_num
    do i = 0,1 
      if (i == 0) then 
        ak = al
        bk = bl
        ck = cl
        pk =  pl
        prat = p / pk
      else 
        ak = ar
        bk = br
        ck = cr
        pk = pr
        prat = p / pk
      endif 

      if (p > pk) then
        f = f + (p - pk) * sqrt(ak / (p + bk))   
      else 
        f = f + g4 * ck * (prat**g1 - 1.0_num)
      endif

    enddo 

    f = f + (ur-ul)

    return
  end function

  real(num) function fprime(p) ! 1st derivative
    real(num), intent(in) :: p
    real(num) :: ak,bk,ck, pk, rhok, prat
    integer :: i 

    fprime = 0.0_num

    do i = 0,1 

      if (i == 0) then 
        ak = al
        bk = bl
        ck = cl
        pk =pl
        rhok = rhol
        prat = p / pk
      else 
        ak = ar
        bk = br
        ck = cr
        pk = pr
        rhok = rhor
        prat = p / pk
      endif 

      if (p > pk) then
        fprime = fprime + &
          & sqrt(ak / bk + p) * (1.0_num - 0.5 * (p-pk)/(bk+p)) 
      else 
        fprime = fprime + prat**(-g2) / (rhok * ck)
      endif

    enddo 

    return
  end function

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

