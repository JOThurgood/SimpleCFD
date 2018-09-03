! minimal working example of the pstar calculation following the 
! algorithm of toro - shows his table is indeed correct and 
! I was just being rather cavalier trying to 
! embedd it within a more complicated program structure from the get-go!

program pstar

  implicit none

  integer, parameter :: num=selected_real_kind(p=6)

  integer :: i, k

  real(num) :: rhol, rhor, rhok
  real(num) :: ul, ur
  real(num) :: pl, pr, pk
  real(num) :: cl, cr
  real(num) :: ak, bk, ck
 
  real(num) :: ps, pold, delta, change
  real(num) :: f, fprime
  real(num) :: gamma
  real(num) :: tol

 !initial conditions

  !test 1 sod
  rhol = 1.0_num 
  ul = 0.0_num
  pl = 1.0_num
  rhor = 0.125_num
  ur = 0.0_num
  pr = 0.1_num
  gamma = 1.4_num

!  !test 5 - both halfs of the W+C blast 
!  rhol = 5.99924_num
!  ul = 19.5975_num
!  pl = 460.894_num
!  rhor= 5.99242_num
!  ur = -6.19633_num
!  pr = 46.0950_num
!  gamma = 1.4_num

  !initial guess on pstar

  ps = 0.5_num * (pl + pr) 

  !iterate NR for converged pstar according to tolerances
  tol = 1e-6_num 

  i = 0  
  change = 1e15_num
  do
    if (change <= tol) exit !condition on tolerance

    pold = ps 

    !evaluate the function f and fprime at pold
    f = ur - ul !speed contribution to f(p)
    fprime = 0.0_num
    do k = 0, 1 !left and right wave contributions

      if (k == 0) then
        ak = 2.0_num / (gamma + 1.0_num) / rhol
        bk = (gamma-1.0_num)/(gamma+1.0_num) * pl
        ck = sqrt(gamma * pl / rhol)
        pk = pl
        rhok = rhol
      else 
        ak = 2.0_num / (gamma + 1.0_num) / rhor
        bk = (gamma-1.0_num)/(gamma+1.0_num) * pr
        ck = sqrt(gamma * pl / rhor)
        pk = pr
        rhok = rhor
      endif 

      if (ps > pk) then !shock
        f = f + (ps - pk) * sqrt(ak / (ps+bk))
        fprime = fprime + sqrt(ak / (ps+bk)) * &
          & (1.0_num - 0.5_num*(ps-pk)/(bk+ps))
      else !rarefaction
        f = f + 2.0_num / (gamma-1.0_num) * ck * & 
          & ((ps/pk)**((gamma-1.0_num)/2.0_num/gamma) - 1.0_num)
        fprime = fprime + (ps/pk)**(-(gamma+1.0_num)/2.0_num/gamma) / &
          & (rhok * ck)
      endif 

    enddo
    
    !advance the newton raphson and calculate a change to compare to tol
    delta = f /fprime
    ps = pold - delta
    change = 2.0_num * abs((ps - pold) / (ps + pold)) 

    i = i + 1
    if (ps < 0.0_num) ps = tol !to correct -ve p guesses 

    print *,'i',i,'pold',pold,'psnew',ps, 'change',change,'delta', delta

    !if (i > ) exit
  enddo 
  
print *, 'Newton-Raphson converged in',i,'iterations'


end program pstar
