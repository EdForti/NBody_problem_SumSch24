module nbody_module
  implicit none
  real, parameter :: SOFTENING = 1.0e-9

  type :: Body
     real :: x, y, z, vx, vy, vz
  end type Body
contains

  subroutine randomizeBodies(data, n)
    real, dimension(:), intent(out) :: data
    integer, intent(in) :: n
    integer :: i
    real :: r

    call random_seed()
    do i = 1, n
       call random_number(r)
       data(i) = 2.0 * (r - 0.5)
    end do
  end subroutine randomizeBodies

  subroutine bodyForce(p, dt, n, Fx, Fy, Fz)
    type(Body), dimension(:), intent(inout) :: p
    real, intent(in) :: dt
    integer, intent(in) :: n
    real, dimension(:), intent(inout) :: Fx, Fy, Fz
    integer :: i, j
    real :: dx, dy, dz, distSqr, invDist, invDist3

    do i = 1, n
       Fx(i) = 0.0
       Fy(i) = 0.0
       Fz(i) = 0.0
       do j = 1, n
          if (i /= j) then
             dx = p(j)%x - p(i)%x
             dy = p(j)%y - p(i)%y
             dz = p(j)%z - p(i)%z
             distSqr = dx * dx + dy * dy + dz * dz + SOFTENING
             invDist = 1.0 / sqrt(distSqr)
             invDist3 = invDist * invDist * invDist
             Fx(i) = Fx(i) + dx * invDist3
             Fy(i) = Fy(i) + dy * invDist3
             Fz(i) = Fz(i) + dz * invDist3
          end if
       end do
       p(i)%vx = p(i)%vx + dt * Fx(i)
       p(i)%vy = p(i)%vy + dt * Fy(i)
       p(i)%vz = p(i)%vz + dt * Fz(i)
    end do
  end subroutine bodyForce

  subroutine saveForcesToFile(filename, nBodies, p, Fx, Fy, Fz)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: nBodies
    type(Body), dimension(:), intent(in) :: p
    real, dimension(:), intent(in) :: Fx, Fy, Fz
    integer :: i
    open(unit=10, file=filename, status='unknown')
    do i = 1, nBodies
       write(10, '(A, I0, A, F6.3, A, F6.3, A, F6.3, A, F6.3, A, F6.3, A, F6.3)') &
            'Body ', i, ': x = ', p(i)%x, ', y = ', p(i)%y, ', z = ', p(i)%z, &
            ', Fx = ', Fx(i), ', Fy = ', Fy(i), ', Fz = ', Fz(i)
    end do
    close(10)
  end subroutine saveForcesToFile

end module nbody_module

program nbody
  use nbody_module
  implicit none
  integer :: nBodies, nIters, i, iter
  real :: dt, totalTime, avgTime, rate
  type(Body), allocatable :: p(:)
  real, allocatable :: buf(:), Fx(:), Fy(:), Fz(:)
  real :: tStart, tEnd, tIter

  nBodies = 30
  nIters = 5
  dt = 0.01

  allocate(p(nBodies))
  allocate(buf(6 * nBodies))
  allocate(Fx(nBodies), Fy(nBodies), Fz(nBodies))

  call randomizeBodies(buf, 6 * nBodies)
  do i = 1, nBodies
     p(i)%x = buf(6 * (i - 1) + 1)
     p(i)%y = buf(6 * (i - 1) + 2)
     p(i)%z = buf(6 * (i - 1) + 3)
     p(i)%vx = buf(6 * (i - 1) + 4)
     p(i)%vy = buf(6 * (i - 1) + 5)
     p(i)%vz = buf(6 * (i - 1) + 6)
  end do

  totalTime = 0.0

  do iter = 1, nIters
     call cpu_time(tStart)
     call bodyForce(p, dt, nBodies, Fx, Fy, Fz)

     do i = 1, nBodies
        p(i)%x = p(i)%x + p(i)%vx * dt
        p(i)%y = p(i)%y + p(i)%vy * dt
        p(i)%z = p(i)%z + p(i)%vz * dt
     end do

     call cpu_time(tEnd)

        if (iter > 1) then
        tIter = tEnd - tStart
        totalTime = totalTime + tIter
        print 100, 'Iteration ', iter, ': ', tIter, ' seconds, Total time: ', totalTime, " seconds"
        end if

  100 format(A,I5,A,F14.8,A,F14.8,A)
  end do

  call saveForcesToFile('forces.txt', nBodies, p, Fx, Fy, Fz)

  avgTime = totalTime / real(nIters - 1)
  rate = real(nBodies) / avgTime

  print 101, 'Average time for iteration is: ', avgTime, ' sec'
  print 102, 'Average rate for iterations 2 through ', nIters, ': ', rate, ' steps per second.'
  print 103, '#Bodies: ',nBodies,', Average ', 1.0e-9 * nBodies * nBodies / avgTime, ' Billion Interactions/second'

  101 format(/,A,F14.8,A)
  102 format(A,I5,A,F12.2,A)
  103 format(A,I5,A,F10.6,A,/)

  deallocate(p)
  deallocate(Fx)
  deallocate(Fy)
  deallocate(Fz)
end program nbody

