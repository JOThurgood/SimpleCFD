module diagnostics

  use shared_data

  implicit none

  real(num), allocatable, dimension(:,:) :: utestarr, vtestarr

  private

  public :: test_analytic_sln, analytic_sln_plots 

  contains

  subroutine test_analytic_sln

    real(num) :: x,y,t
    real(num) :: utest, vtest
    real(num) :: L2


    ! only allocate on first call 
    if (.not. allocated(utestarr)) allocate(utestarr(1:nx,1:ny))
    if (.not. allocated(vtestarr)) allocate(vtestarr(1:nx,1:ny))

    t = time
    L2 = 0.0_num
    do ix = 1, nx
    do iy = 1, ny
      x = xc(ix)
      y = yc(iy)
      utest =1.0_num - 2.0_num * cos(2.0_num*pi*(x-t)) * sin(2.0_num*pi*(y-t))
      vtest =1.0_num + 2.0_num * sin(2.0_num*pi*(x-t)) * cos(2.0_num*pi*(y-t))

      utestarr(ix,iy) = utest
      vtestarr(ix,iy) = vtest

      L2 = L2 + abs(u(ix,iy)-utest)**2 + abs(v(ix,iy)-vtest)**2
    enddo
    enddo
    L2 = sqrt( L2 / real(2*nx*ny,num))

    print *,'Cycle: ',step,'time: ',time,'L2 of velocity vs analytic',L2

  end subroutine test_analytic_sln

  subroutine analytic_sln_plots

    print *,'******************************************************************'
    print *, 'Drawing plots'
    print *,'******************************************************************'

    call execute_command_line("rm -rf *.dat *.png")

    open(out_unit, file="xc.dat", access="stream")
    write(out_unit) xc(1:nx)
    close(out_unit)
   
    open(out_unit, file="yc.dat", access="stream")
    write(out_unit) yc(1:ny)
    close(out_unit)
  
    open(out_unit, file="u.dat", access="stream")
    write(out_unit) u(1:nx,1:ny)
    close(out_unit)
  
    open(out_unit, file="v.dat", access="stream")
    write(out_unit) v(1:nx,1:ny)
    close(out_unit)

    open(out_unit, file="utest.dat", access="stream")
    write(out_unit) utestarr(1:nx,1:ny)
    close(out_unit)
  
    open(out_unit, file="vtest.dat", access="stream")
    write(out_unit) vtestarr(1:nx,1:ny)
    close(out_unit)

    open(out_unit, file="udiff.dat", access="stream")
    write(out_unit) u(1:nx,1:ny) - utestarr(1:nx,1:ny)
    close(out_unit)
  
    open(out_unit, file="vdiff.dat", access="stream")
    write(out_unit) v(1:nx,1:ny) - vtestarr(1:nx,1:ny)
    close(out_unit)



    call execute_command_line("python analytic_sln_plots.py")
    call execute_command_line("rm -rf xc.dat yc.dat u.dat v.dat utest.dat vtest.dat udiff.dat vdiff.dat")

  end subroutine analytic_sln_plots

end module diagnostics
