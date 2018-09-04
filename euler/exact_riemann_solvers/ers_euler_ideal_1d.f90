! ers_euler_ideal_1d.f90
!
! Exact Riemann Solver for ideal 1D euler equations, producing a 
! converged exact solution.
! 
! <link to wiki> 

module shared_data ! a common block essentially

  implicit none

  integer, parameter :: num=selected_real_kind(p=6) !p=6 or 15
  integer :: niter !debug

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
    al = g5 / rhol !only used in pressure function - move out of shared?
    ar = g5 / rhor
    bl = g6 * pl
    br = g6 * pr
    !sound speeds
    cl = sqrt(gamma * pl /rhol)
    cr = sqrt(gamma * pr / rhor)
  end subroutine constants


end module shared_data

module user

  use shared_data 

  implicit none 

  public 

  contains 

  subroutine initial_conditions

    ! set these for your problem - default is Sod
    rhol = 1.0_num
    ul = 0.0_num
    pl = 1.0_num
    rhor = 0.125_num  
    ur = 0.0_num
    pr = 0.10_num 
    gamma = 1.4_num

  end subroutine initial_conditions
  
  subroutine control
    !stub for setting sampling and output options
  end subroutine 

end module user 

module riemann !subroutines related to calculating star states

  use shared_data

  implicit none

  private 

  public :: check_positivity, pstar, newton_raphson 

  !if when all is said and done there are any vars only used at 
  ! this stage it would be sensible to move here from shared_data

  contains 


  subroutine check_positivity
    !stub
  end subroutine check_positivity

 subroutine pstar
    call guess_ps
    call newton_raphson
    print *, ps
  end subroutine pstar

  subroutine guess_ps !initial guess on pstar 
    ! in the interests of speed you, can base your guess on the value of
    ! pstar as the  given by an appropriate choice of approximate
    ! R solver. This should mean that the Newton-Raphson iterates in as
    ! few steps as possible. 

    ps = 0.5_num * (pl + pr) ! For now just always use arithmetic mean
          !and see if you get same convergence as in Table 4.2 of Toro
  end subroutine guess_ps
 
  subroutine newton_raphson

    integer :: i 
    real(num) :: pold, delta
    real(num) :: tol = 1.0e-6_num
    real(num) :: rpc  

    i = 0 !define these here instead of header so can call subroutine
    rpc = 1e15_num !repeatedly during testing

    do
      if (rpc <= tol) exit !condition on tolerance
      pold = ps 
      delta = f(pold)/fprime(pold)
      ps = pold - delta !f(pold)/fprime(pold)
      rpc = 2.0_num * abs((ps - pold) / (ps + pold)) 
      i = i + 1
      if (ps < 0.0_num) ps = tol !to correct -ve p guesses 
!      print *,'i',i,'pold',pold,'psnew',ps, 'rpc',rpc,'delta', delta
      if (i > 20) print *,'nonconvergence'
      if (i > 20) exit
    enddo 

!    print *, 'Newton-Raphson converged in',i,'iterations'
    niter = i !so test module can keep track
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

    f = f + (ur - ul)

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
        pk = pl
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
          & sqrt(ak / (bk + p)) * (1.0_num - 0.5_num * (p-pk)/(bk+p)) 
      else 
        fprime = fprime + prat**(-g2) / (rhok * ck)
      endif

    enddo 

    return
  end function

end module riemann

module tests !subroutines for automatic testing

  use shared_data
  use riemann 

  implicit none


  contains

  subroutine test_pstar !check pstar and iteration info
    logical :: test1 = .true.
    logical :: verbose = .false.
    real(num) :: tol = 1e-6_num 
    real(num), dimension(5) :: p_toro
    real(num), dimension(4,5) :: p0_toro
    integer, dimension(4,5) :: i_toro
    integer :: i

    !toro
    p_toro = (/0.30313_num, 0.00189_num, 460.894_num, 46.0950_num, &
      & 1691.64_num /)

    i_toro = reshape( (/3,5,3,5, &
                      1, 8, 8 ,9, &
                      5, 4, 3, 4, &
                      5, 4, 3, 4, &
                      4, 5, 4, 6/), shape(i_toro))!, order = (/2,1/) )
 
    p0_toro = reshape((/0.30677_num, 0.55_num, 0.31527_num, 0.55_num, &
      & 0.00189387344_num, tol, tol, 0.4_num, &
      & 912.449_num, 500.005_num, 464.108_num, 500.005_num, &  
      & 82.9831_num, 50.005_num, 46.4162_num, 50.005_num, &
      & 2322.65_num, 781.353_num, 1241.21_num, 253.494_num/), &
      & shape(p0_toro) )

!    do i=1,5 
!      print *,'i_toro',i_toro(:,i) 
!      print *,'p0_toro',p0_toro(:,i) 
!    enddo

    call test_1 ! initialise  Sods' test
    do i = 1, 4 
      ps = p0_toro(i,1)
      call newton_raphson 
      if (verbose) & 
        & print *, 'test 1: ','p0', p0_toro(i,1), 'ps',ps,'niter',niter
      if (niter .ne. i_toro(i,1)) then
        print *,'fail'
        test1 = .false.
      endif
    enddo     

    call test_2
    do i = 1, 4 
      ps = p0_toro(i,2)
      call newton_raphson 
      if (verbose) & 
        & print *, 'test 2: ','p0', p0_toro(i,2), 'ps',ps,'niter',niter
      if (niter .ne. i_toro(i,2)) print *,'fail'
    enddo     

    call test_3
    do i = 1, 4 
      ps = p0_toro(i,3)
      call newton_raphson 
      if (verbose) & 
        & print *, 'test 3: ','p0', p0_toro(i,3), 'ps',ps,'niter',niter
      if (niter .ne. i_toro(i,3)) print *,'fail'
    enddo     

    call test_4
    do i = 1, 4 
      ps = p0_toro(i,4)
      call newton_raphson 
      if (verbose) & 
        & print *, 'test 4: ','p0', p0_toro(i,4), 'ps',ps,'niter',niter
      if (niter .ne. i_toro(i,4)) print *,'fail'
    enddo     

    call test_5
    do i = 1, 4 
      ps = p0_toro(i,5)
      call newton_raphson 
      if (verbose) & 
        & print *, 'test 5: ','p0', p0_toro(i,5), 'ps',ps,'niter',niter
      if (niter .ne. i_toro(i,5)) print *,'fail'
    enddo     

    if (test1) then 
      print *,'passed pstar test'
    else 
      print *,'failed pstar test'
    endif

  end subroutine test_pstar


  ! predefined initial conditions

  subroutine test_1 !Sod's test
    rhol = 1.0_num
    ul = 0.0_num
    pl = 1.0_num
    rhor = 0.125_num  
    ur = 0.0_num
    pr = 0.10_num 
    gamma = 1.4_num
    call constants !reset constants
  end subroutine test_1 

  subroutine test_1_r !Reverse Sod's test for symmetry test in p*
    rhor = 1.0_num
    ur = 0.0_num
    pr = 1.0_num
    rhol = 0.125_num !1.0_num / 8.0_num 
    ul = 0.0_num
    pl = 0.10_num 
    gamma = 1.4_num
    call constants !reset constants
  end subroutine test_1_r

  subroutine test_2 !1,2,3 problem
    rhol = 1.0_num
    ul = -2.0_num
    pl = 0.4_num
    rhor = 1.0_num
    ur = 2.0_num
    pr = 0.4_num 
    gamma = 1.4_num
    call constants !reset constants
  end subroutine test_2 

  subroutine test_3 !left half blast
    rhol = 1.0_num
    ul = 0.0_num
    pl = 1000.0_num
    rhor = 1.0_num
    ur = 0.0_num
    pr = 0.01_num
    gamma = 1.4_num
    call constants !reset constants
  end subroutine test_3

  subroutine test_4 !right half
    rhol = 1.0_num
    ul = 0.0_num
    pl = 0.01_num
    rhor = 1.0_num
    ur = 0.0_num
    pr = 100.0_num
    gamma = 1.4_num
    call constants !reset constants
  end subroutine test_4

  subroutine test_5 !both halfs of the W+C blast 
    rhol = 5.99924_num
    ul = 19.5975_num
    pl = 460.894_num
    rhor= 5.99242_num
    ur = -6.19633_num
    pr = 46.0950_num
    gamma = 1.4_num
    call constants !reset constants
  end subroutine 

end module tests



program ers_euler_ideal_1d

  use shared_data
  use riemann 
  use user 
  use tests
  implicit none   

  !check passes numerical tests
  call test_pstar 

  !do calculations for users setup 
  call initial_conditions
!  call control 
  call constants 
! call check_positivity

  call pstar


end program ers_euler_ideal_1d
