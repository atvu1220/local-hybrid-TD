module boundary
      implicit none
      contains

      subroutine boundary_vector(b)
            use dimensions
            use inputs, only: boundx
            implicit none
            real, intent(inout):: b(nx,ny,nz,3)

            if (boundx .eq. 1) then
                  call periodic(b)
            elseif (boundx .eq. 2) then
                  call periodic_xy(b)
            elseif (boundx .eq. 4) then
                  call periodic_y(b)
            elseif (boundx .eq. 5) then
                  call periodic_yz(b)    
            else
                  write(*,*) 'Boundary conditions are unspecified'
                  stop
            endif
            
      end subroutine boundary_vector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      subroutine boundary_scalar(b)
            use dimensions
            use inputs, only: boundx
            implicit none
            real, intent(inout):: b(nx,ny,nz)

            if (boundx .eq. 1) then
                  call periodic_scalar(b)
            elseif (boundx .eq. 2) then
                  call periodic_scalar_xy(b)
            elseif (boundx .eq. 4) then
                  call periodic_scalar_y(b)
            elseif (boundx .eq. 5) then
                  call periodic_scalar_yz(b)     
            else
                  write(*,*) 'Boundary conditions are unspecified'
                  stop
            endif
            
      end subroutine boundary_scalar
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      subroutine add_boundary_vector(b)
            use dimensions
            use inputs, only: boundx
            implicit none
            real, intent(inout):: b(nx,ny,nz,3)
            if (boundx .eq. 1) then
                  b(nx-1,:,:,:) = b(nx-1,:,:,:)+b(1,:,:,:)
                  b(:,ny-1,:,:) = b(:,ny-1,:,:)+b(:,1,:,:)
                  b(:,:,nz-1,:) = b(:,:,nz-1,:)+b(:,:,1,:)
            endif
            if (boundx .eq. 2) then
                  b(nx-1,:,:,:) = b(nx-1,:,:,:)+b(1,:,:,:)
                  b(:,ny-1,:,:) = b(:,ny-1,:,:)+b(:,1,:,:)
            endif
            if (boundx .eq. 4) then
                  b(:,ny-1,:,:) = b(:,ny-1,:,:)+b(:,1,:,:)
            endif
      end subroutine add_boundary_vector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      subroutine add_boundary_scalar(b)
            use dimensions
            use inputs, only: boundx
            implicit none
            real, intent(inout):: b(nx,ny,nz)
            if (boundx .eq. 1) then
                  b(nx-1,:,:) = b(nx-1,:,:)+b(1,:,:)
                  b(:,ny-1,:) = b(:,ny-1,:)+b(:,1,:)
                  b(:,:,nz-1) = b(:,:,nz-1)+b(:,:,1)
            endif
            if (boundx .eq. 2) then
                  b(nx-1,:,:) = b(nx-1,:,:)+b(1,:,:)
                  b(:,ny-1,:) = b(:,ny-1,:)+b(:,1,:)
            endif
            if (boundx .eq. 4) then
                  b(:,ny-1,:) = b(:,ny-1,:)+b(:,1,:)
            endif
      end subroutine add_boundary_scalar            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
      subroutine periodic(b)
            use dimensions
            implicit none
            real, intent(inout):: b(nx,ny,nz,3)
            integer:: i,j,k,m
            
!       X direction

            do j=1,ny
                  do k=1,nz
                        do m=1,3
                              b(1,j,k,m) = b(nx-1,j,k,m)
                              b(nx,j,k,m) = b(2,j,k,m)
                        enddo
                  enddo
            enddo

!       Y direction

            do i=1,nx
                  do k=1,nz
                        do m=1,3
                              b(i,1,k,m) = b(i,ny-1,k,m)
                              b(i,ny,k,m) = b(i,2,k,m)
                        enddo
                  enddo
            enddo

!       Z direction

            do i=1,nx
                  do j=1,ny
                        do m=1,3
                              b(i,j,1,m) = b(i,j,nz-1,m)
                              b(i,j,nz,m) = b(i,j,2,m)
                        enddo
                  enddo
            enddo

      end subroutine periodic
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine periodic_xy(b)
            use dimensions
            implicit none
            real, intent(inout):: b(nx,ny,nz,3)
            integer:: i,j,k,m
            
!       X direction

            do j=1,ny
                  do k=1,nz
                        do m=1,3
                              b(1,j,k,m) = b(nx-1,j,k,m)
                              b(nx,j,k,m) = b(2,j,k,m)
                        enddo
                  enddo
            enddo

!       Y direction

            do i=1,nx
                  do k=1,nz
                        do m=1,3
                              b(i,1,k,m) = b(i,ny-1,k,m)
                              b(i,ny,k,m) = b(i,2,k,m)
                        enddo
                  enddo
            enddo
            
!       Z direction is not periodic
            do i=1,nx
                  do j=1,ny
                        do m=1,3
                              b(i,j,nz,m) = b(i,j,nz-1,m)
                              b(i,j,1,m) = b(i,j,2,m)
                        enddo
                  enddo
            enddo
            
      end subroutine periodic_xy
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine periodic_y(b)
                  use dimensions
                  implicit none
                  real, intent(inout):: b(nx,ny,nz,3)
                  integer:: i,j,k,m
                  
      !       X direction is not periodic

                  do j=1,ny
                        do k=1,nz
                              do m=1,3
                                    b(nx,j,k,m) = b(nx-1,j,k,m)
                              enddo
                        enddo
                  enddo

      !       Y direction

                  do i=1,nx
                        do k=1,nz
                              do m=1,3
                                    b(i,1,k,m) = b(i,ny-1,k,m)
                                    b(i,ny,k,m) = b(i,2,k,m)
                              enddo
                        enddo
                  enddo
                  
      !       Z direction is not periodic
                  do i=1,nx
                        do j=1,ny
                              do m=1,3
                                    b(i,j,nz,m) = b(i,j,nz-1,m)
                                    b(i,j,1,m) = b(i,j,2,m)
                              enddo
                        enddo
                  enddo
                  
            end subroutine periodic_y
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine periodic_yz(b)
                  use dimensions
                  implicit none
                  real, intent(inout):: b(nx,ny,nz,3)
                  integer:: i,j,k,m
                  
      !       X direction is not periodic
                  do j=1,ny
                        do k=1,nz
                              do m=1,3
                                    !b(1,j,k,m) = b(2,j,k,m)
                                    b(nx,j,k,m) = b(nx-1,j,k,m)
                              enddo
                        enddo
                  enddo

      !       Y direction is periodic
                  do i=1,nx
                        do k=1,nz
                              do m=1,3
                                    b(i,1,k,m) = b(i,ny-1,k,m)
                                    b(i,ny,k,m) = b(i,2,k,m)
                              enddo
                        enddo
                  enddo
                  
      !       Z direction is periodic
                  do i=1,nx
                        do j=1,ny
                              do m=1,3
                                    b(i,j,nz,m) = b(i,j,2,m)
                                    b(i,j,1,m) = b(i,j,nz-1,m)
                              enddo
                        enddo
                  enddo
                  
end subroutine periodic_yz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine periodic_scalar(b)
            use dimensions
            implicit none
            real, intent(inout):: b(nx,ny,nz)
            integer:: i,j,k
            
!       X surfaces
            do j=1,ny
                  do k=1,nz
                        b(1,j,k) = b(nx-1,j,k)
                        b(nx,j,k) = b(2,j,k)
                  enddo
            enddo

!       Y surfaces
            do i=1,nx
                  do k=1,nz
                        b(i,1,k) = b(i,ny-1,k)
                        b(i,ny,k) = b(i,2,k)
                  enddo
            enddo
            
!       Z surfaces
            do i=1,nx
                  do j=1,ny
                        b(i,j,1) = b(i,j,nz-1)
                        b(i,j,nz) = b(i,j,2)
                  enddo
            enddo    
            
      end subroutine periodic_scalar
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine periodic_scalar_xy(b)
            use dimensions
            implicit none
            real, intent(inout):: b(nx,ny,nz)
            integer:: i,j,k
            
!       X surfaces
            do j=1,ny
                  do k=1,nz
                        b(1,j,k) = b(nx-1,j,k)
                        b(nx,j,k) = b(2,j,k)
                  enddo
            enddo

!       Y surfaces
            do i=1,nx
                  do k=1,nz
                        b(i,1,k) = b(i,ny-1,k)
                        b(i,ny,k) = b(i,2,k)
                  enddo
            enddo
            
!       Z surfaces are not periodic
            do i=1,nx
                  do j=1,ny
                        b(i,j,1) = b(i,j,2)
                        b(i,j,nz) = b(i,j,nz-1)
                  enddo
            enddo    
            
      end subroutine periodic_scalar_xy
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine periodic_scalar_y(b)
            use dimensions
            implicit none
            real, intent(inout):: b(nx,ny,nz)
            integer:: i,j,k
            
!       X surfaces are not periodic
            do j=1,ny
                  do k=1,nz
                        b(nx,j,k) = b(nx-1,j,k)
                  enddo
            enddo

!       Y surfaces
            do i=1,nx
                  do k=1,nz
                        b(i,1,k) = b(i,ny-1,k)
                        b(i,ny,k) = b(i,2,k)
                  enddo
            enddo
            
!       Z surfaces are not periodic
            do i=1,nx
                  do j=1,ny
                        b(i,j,1) = b(i,j,2)
                        b(i,j,nz) = b(i,j,nz-1)
                  enddo
            enddo
            
      end subroutine periodic_scalar_y
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine periodic_scalar_yz(b)
            use dimensions
            implicit none
            real, intent(inout):: b(nx,ny,nz)
            integer:: i,j,k
            
!       X surfaces are not periodic
            do j=1,ny
                  do k=1,nz
                        !b(1,j,k) = b(2,j,k)
                        b(nx,j,k) = b(nx-1,j,k)
                  enddo
            enddo

!       Y surfaces is periodic
            do i=1,nx
                  do k=1,nz
                        b(i,1,k) = b(i,ny-1,k)
                        b(i,ny,k) = b(i,2,k)
                  enddo
            enddo
            
!       Z surfaces are periodic
            do i=1,nx
                  do j=1,ny
                        b(i,j,1) = b(i,j,nz-1)
                        b(i,j,nz) = b(i,j,2)
                  enddo
            enddo
            
      end subroutine periodic_scalar_yz      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
!!!!!!!!!RANDOM NUMBER GENERATOR!!!!!!!!!!!!!!

      real function pad_ranf()
            implicit none
            call random_number(pad_ranf)
      end function pad_ranf
      
      
      subroutine random_initialize(seed_input)
!**********************************************************************
!
!! RANDOM_INITIALIZE initializes the FORTRAN90 random number seed.
!
!  Discussion:
!    If you don't initialize the FORTRAN90 random number generator
!    routine RANDOM_NUMBER, which is used by calls like
!
!      call random_number ( harvest = x )
!
!    then its behavior is not specified.  That may be OK for you.  But
!    you may want to be able to force the same sequence to be generated
!    each time, or to force a different sequence.
!
!    To control the sequence of random numbers, you need to set the seed.
!    In FORTRAN90, this is done by calling the RANDOM+SEED routine.
!    You can call it with no arguments, in fact.  But if you call
!    it with no arguments:
!
!      call random_seed ( )
!
!    then its behavior (or more particularly, the behavior of RANDOM_NUMBER)
!    is still not specified.  You might hope that the system will go out
!    and pick a nice random seed for you, but there's no guarantee.
!
!
!    For example, on the DEC ALPHA, if you compile a program that calls
!    RANDOM_NUMBER, then every time you run it, you get the same sequence
!    of "random" values.  If you compile a program that calls RANDOM_SEED
!    with no arguments, and then calls RANDOM_NUMBER, you still get the
!    same sequence each time you run the program.
!
!    In order to actually try to scramble up the random number generator
!    a bit, this routine goes through the tedious process of getting the
!    size of the random number seed, making up values based on the current
!    time, and setting the random number seed.
!
!    Unfortunately, the RANDOM_SEED routine has a very elastic definition.
!    It does not use a single scalar integer SEED.  Instead, it communicates
!    with the user through an integer vector whose size is not specified.
!    You actually have to "ask" the routine to tell you the size of this
!    vector.  Then you can fill up the vector with values to be used to
!    influence the seeding of the random number routine.  The details of
!    how the seed affects the sequence are also unspecified, but the only
!    thing we are reasonably confident about is that specifying the same
!    seed should result in the same sequence, and specifying different
!    seeds should result in different sequences!
!
!    I assume this is far more than you wanted to know.  (It's certainly
!    more than I wanted to know!)
!
!    The argument SEED is an input quantity only, so it is legal
!    to type
!
!      call random_initialize ( 0 )
!
!    or
!
!      call random_initialize ( 18867 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!     Input, integer ( kind = 4 ) SEED_INPUT, the user's "suggestion" for a seed
!     However, if the input value is 0, the routine will come up with
!     its own "suggestion", based on the system clock.   

            implicit none
            integer, intent(in):: seed_input
            integer:: seed, count, count_max, count_rate, seed_size
            integer, allocatable:: seed_vector(:)
            logical, parameter:: debug = .false.
            real:: t
            integer, parameter:: warm_up = 100
            integer:: i
            
            seed = seed_input
            
            !initialize the random seed routine
            call random_seed()
            !Determine the size of the random number seed vector
            call random_seed(size = seed_size)
            !Allocate a vector of the right size to be used as a random seed
            allocate (seed_vector(seed_size))
            !If user supplied a SEED value, use that
            
            !Otherwise, use the system clock value to make up a value that is likely to change based on when the routine is called.
            
            if (seed /= 0 ) then
                  if (debug) then
                        write(*,*) ' '
                        write(*,*) 'RANDOM_INTITIALIZE'
                        write(*,*) 'Initialize RANDOM_NUMBER, user SEED = ', seed
                  endif
            else
                  call system_clock(count,count_rate,count_max)
                  
                  seed = count
                  
                  if (debug) then
                        write(*,*) ' '
                        write(*,*) 'RANDOM_INITIALIZE'
                        write(*,*) 'Initialize RANDOM_NUMBER, arbitrary SEED = ',seed
                  endif
            endif
            
            !Set the seed vector.  Set all entries to seed
            
            seed_vector(1:seed_size) = seed
            
            !Free up the seed space
            
            deallocate(seed_vector)
            
            !Call the random number routine a bunch of time to "warm up"
            
            do i = 1, warm_up
                  call random_number(harvest = t)
            enddo
            
      end subroutine
                  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module boundary

