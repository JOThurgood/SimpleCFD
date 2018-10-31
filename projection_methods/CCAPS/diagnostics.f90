module diagnostics

  use shared_data

  implicit none

  private

  public :: test_analytic_sln 

  contains

  subroutine test_analytic_sln

    real(num) :: x,y,t
    real(num) :: utest, vtest
    real(num) :: L2

    L2 = 0.0_num
    do ix = 1, nx
    do iy = 1, ny
      t = time
      x = xc(ix)
      y = yc(iy)
      utest =1.0_num - 2.0_num * cos(2.0_num*pi*(x-t)) * sin(2.0_num*pi*(y-t))
      vtest =1.0_num + 2.0_num * sin(2.0_num*pi*(x-t)) * cos(2.0_num*pi*(y-t))

      L2 = L2 + (u(ix,iy)-utest)**2 + (v(ix,iy)-vtest)**2
    enddo
    enddo
    L2 = sqrt( L2 / real(2*nx*ny,num))

    print *,'Cycle: ',step,'time: ',time,'L2 of velocity vs analytic',L2

  end subroutine test_analytic_sln

end module diagnostics
