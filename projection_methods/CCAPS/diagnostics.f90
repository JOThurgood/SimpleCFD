module diagnostics

  use shared_data

  implicit none

  real(num), allocatable, dimension(:,:) :: utestarr, vtestarr
  real(num), allocatable, dimension(:,:) :: curlu

  private

  public :: test_minion, minion_plots, sln_plots, &
    & plot_divergence_now, test_steady

  contains

  subroutine test_steady
    ! try to assess if solution has converged to a steady state.
    ! Involves keeping a copy of the previous solution
    ! which would be quite wasteful at scale

    logical :: first_call = .true.
    real(num), dimension(1:nx,1:ny) :: uold, vold
    real(num) :: L2 

    if (first_call) then
      uold = u(1:nx,1:ny)
      vold = v(1:nx,1:ny)
      first_call = .false.
    else
    
      L2 = 0.0_num
      do ix = 1, nx
      do iy = 1, ny
        L2 = L2 + abs(u(ix,iy)-uold(ix,iy))**2 + abs(v(ix,iy)-vold(ix,iy))**2
      enddo
      enddo
      L2 = sqrt( L2 / real(2*nx*ny,num))

      uold = u(1:nx,1:ny)
      vold = v(1:nx,1:ny)

    endif

    print *,''
    print *,' L2 of velocity vs previous step:',L2
    print *,''

  end subroutine test_steady


  subroutine test_minion

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

    print *,''
    print *,' Morrisons Convergence test: - L2 of velocity vs analytic sln:',L2
    print *,''

  end subroutine test_minion

  subroutine sln_plots

    print *,'******************************************************************'
    print *, 'Drawing plots'
    print *,'******************************************************************'

    call execute_command_line("rm -rf *.dat *.png")

    open(out_unit, file="time.dat", access="stream")
    write(out_unit) time
    close(out_unit)

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

    open(out_unit, file="divu.dat", access="stream")
    write(out_unit) divu(1:nx,1:ny)
    close(out_unit)

    call get_vorticity

    open(out_unit, file="curlu.dat", access="stream")
    write(out_unit) curlu(1:nx,1:ny)
    close(out_unit)


    !open(out_unit, file="phi.dat", access="stream")
    !write(out_unit) phi(1:nx,1:ny)
    !close(out_unit)

    call execute_command_line("python python_plotting/sln_plots.py")

    call execute_command_line("rm -rf *.dat")

!    call execute_command_line("rm -rf xc.dat yc.dat u.dat v.dat utest.dat vtest.dat udiff.dat vdiff.dat")

  end subroutine sln_plots

  subroutine minion_plots

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

    open(out_unit, file="divu.dat", access="stream")
    write(out_unit) divu(1:nx,1:ny)
    close(out_unit)

    open(out_unit, file="phi.dat", access="stream")
    write(out_unit) phi(1:nx,1:ny)
    close(out_unit)

    call execute_command_line("python python_plotting/minion_plots.py")
    call execute_command_line("rm -rf *.dat")

  end subroutine minion_plots
  
  subroutine get_vorticity !calc z component of curl U

    real(num) :: dv_dx, du_dy

    if (.not. allocated(curlu)) allocate( curlu(1:nx,1:ny))

    !call velocity_bcs ! not currently necessary

    do ix = 1, nx
    do iy = 1, ny
      dv_dx = (v(ix+1,iy)-v(ix-1,iy))/dx/2.0_num
      du_dy = (u(ix,iy+1)-u(ix,iy-1))/dy/2.0_num
      curlu(ix,iy) = dv_dx - du_dy 
    enddo
    enddo
  end subroutine get_vorticity


  subroutine plot_divergence_now
    !for runtime debugging

    print *,'******************************************************************'
    print *, 'call plot_divercence_now'
    print *,'******************************************************************'

    call execute_command_line("rm -rf divu_rt.dat xc.dat yc.dat")

    open(out_unit, file="divu_rt.dat", access="stream")
    write(out_unit) divu(1:nx,1:ny)
    close(out_unit)

    open(out_unit, file="xc.dat", access="stream")
    write(out_unit) xc(1:nx)
    close(out_unit)
   
    open(out_unit, file="yc.dat", access="stream")
    write(out_unit) yc(1:ny)
    close(out_unit)

    call execute_command_line("python python_plotting/plot_divercence_now.py")
    call execute_command_line("rm -rf divu_rt.dat")

    print *,'******************************************************************'
    print *, 'Early stop'
    print *,'******************************************************************'
    STOP

  end subroutine plot_divergence_now


end module diagnostics
