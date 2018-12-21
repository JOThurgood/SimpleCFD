! can test compiling with and without the openmp flag to see how it goes
!
! Seems to automagically work on my mac
!
!C17775472:OpenMP_test tflv6$ gfortran -o out omp_par_do.f90
!C17775472:OpenMP_test tflv6$ time ./out
! hi
!   313.292236
!
!real	5m13.676s
!user	5m12.808s
!sys	0m0.488s

!
!C17775472:OpenMP tflv6$ gfortran -fopenmp -o out omp_par_do.f90
!C17775472:OpenMP_test tflv6$ time ./out
! hi
! hi
! hi
! hi
! hi
! hi
! hi
! hi
!
!real	1m52.110s
!user	14m12.182s

program omp_par_do

  use omp_lib

  implicit none

  integer, parameter :: n = 5e3
  real, dimension(n) :: dat, result
  real :: start, finish
  integer :: i,j,k
  integer :: threads


  !$OMP PARALLEL
!    threads = omp_get_num_threads() !this will throw an error if not compiled with the ompi flag, otherwise you find it reports faithfully the number of threads
    print *,'hi'
  !$OMP END PARALLEL
!    print *, 'threads', threads

!  call cpu_time(start) !this will give you the total cpu time not a wall or real time
  !$OMP PARALLEL DO
  do i = 1, n
  do j = 1, n
  do k = 1, n
     result(i) = my_function(dat(i))
  end do
  end do
  end do
  !$OMP END PARALLEL DO
!  call cpu_time(finish)

!  print *,finish-start
!  print *,(finish-start)/threads
contains

  function my_function(d) result(y)
    real, intent(in) :: d
    real :: y

    ! do something complex with data to calculate y
  end function my_function

end program omp_par_do
