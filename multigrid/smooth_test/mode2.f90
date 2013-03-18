! Do Gauss-Seidel relaxation on u" = 0 with initial guess of cos(nx), with 3 
! different wave numbers to demonstrate that the shorter wavelengths vanish 
! fastest.

! M. Zingale (2005-03-29)

program fmgfv

  implicit none

  integer, parameter :: nx = 128          ! number of interior zones
  integer, parameter :: ng = 1           ! number of guardcells

  integer, parameter :: nsmooth = 1000   ! number of smoothing blocks

  double precision, parameter :: pi = 3.141592653d0
  
  ! set the left and right boundary conditions
  double precision, parameter :: lbc = 0.d0
  double precision, parameter :: rbc = 0.d0

  ! imin and imax will always point to the starting and ending index of 
  ! the interior zones on the current level
  integer :: imin, imax                   

  double precision :: xmin, xmax, dx, x
  double precision :: source_norm

  character (len=2) imode
  character (len=5) ifile

  integer :: i, j, n
  
  double precision, dimension(:), allocatable :: v, f
  double precision, dimension(:), allocatable :: w

  ! the function g holds the RHS, the function true holds the analytic 
  ! solution.
  double precision :: g, true, error
  double precision :: temp, temp2

  
  ! initialize the solution and rhs arrays
  allocate(f(nx+2*ng))
  allocate(v(nx+2*ng))
  allocate(w(nx+2*ng))

  f(:) = 0.d0
  v(:) = 0.d0

  
  ! setup the grid
  xmin = 0.d0
  xmax = 1.d0
  dx = (xmax - xmin)/dble(nx)

  imin = ng + 1           ! we are using 1 based indexing
  imax = ng + nx 


  ! fill the RHS of the finest level the true RHS
  do i = imin, imax
     x = dble(i - ng - 1 + 0.5d0)*dx + xmin
     v(i) = (sin(2.d0*pi*x)+sin(2.d0*pi*8*x)+sin(2.d0*pi*16*x))/3.d0
     f(i) = 0.0
  enddo


  ! relax
  do j = 1, nsmooth
     call smooth(nx, ng, dx, lbc, rbc, v, f, 1)


     ! compare to the true solution, u = 0
     w(:) = 0.d0
     do i = imin, imax
        w(i) = - v(i)
     enddo

     ! print out every 10th iteration
     if (log10(dble(j)) == int(log10(dble(j))) ) then

        write (ifile,'(i5.5)') j
        open (unit=10, file='relax.mixed.'//ifile, status='unknown')

        do i = imin, imax
           x = dble(i - ng - 1 + 0.5d0)*dx + xmin
           write (10,*) x, v(i)
        enddo

        close (unit=10)

        
     endif

  enddo

  
end program fmgfv
 


!=============================================================================
! error
!=============================================================================

function error(nx, ng, dx, v)

  ! compute the L2 norm 
  
  implicit none

  integer :: nx, ng
  double precision :: dx
  double precision, dimension(nx + 2*ng) :: v

  integer :: i, imin, imax
  double precision error

  imin = ng + 1
  imax = ng + nx

  error = 0.d0

  do i = imin, imax
     error = error + v(i)**2
  enddo
  
  error = dx*error    ! make it grid invariant
  error = sqrt(error)

  return
end function error



!=============================================================================
! smooth
!=============================================================================

subroutine smooth(nx, ng, dx, lbc, rbc, v, f, nsmooth)

  ! given a solution vector, v, and a RHS vector, f,
  ! smooth v to better satisfy the equation.  This is
  ! done in place, using Red-Black Gauss-Seidel

  ! lbc and rbc are the left and right Dirichlet boundary conditions.
  ! Because we are finite-volume, and therefore, cell-centered, we 
  ! need to extrapolate to match the desired Dirichlet BC.

  implicit none

  integer :: nx, ng, nsmooth
  double precision :: dx
  double precision :: lbc, rbc
  double precision, dimension(nx + 2*ng) :: v, f

  integer :: i, m, ioff, color
  integer :: imin, imax

  imin = ng + 1
  imax = ng + nx

  ! do some smoothing -- Red-Black Gauss-Seidel
  do m = 1, nsmooth

    ioff = 0
    do color = 0, 1

       ! set the guardcells to give the proper boundary condition, using
       ! extrapolation
       v(ng) = 2*lbc - v(imin)
       v(ng+nx+1) = 2*rbc - v(imax)

       do i = imin+ioff, imax, 2
          v(i) = 0.5d0*(v(i-1) + v(i+1) - dx*dx*f(i))
       enddo
       
       ioff = 1 - ioff
    enddo

 enddo

  return
end subroutine smooth



