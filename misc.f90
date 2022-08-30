module misc
      implicit none
      contains
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine seed_mpi(my_rank)
            use iso_fortran_env, only:int64
            implicit none
            integer, intent(IN) :: my_rank
            integer,  allocatable :: seed(:)
            integer :: i, n, un, istat, dt(8)
            integer(int64) :: t

            call random_seed(size=n)
            allocate(seed(n))
            ! First check to see if OS provides a random number generator
            open(newunit=un, file="/dev/urandom", access="stream", &
                  form="unformatted", action="read", status="old", iostat=istat)
            if(istat ==0) then
                  read(un) seed
                  close(un)
            else ! Fallback to using time and rank.
                  call system_clock(t)
                  if(t == 0) then
                        call date_and_time(values=dt)
                        t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
                              + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
                              + dt(3) * 24_int64 * 60 * 60 * 1000 &
                              + dt(5) * 60 * 60 * 1000 &
                              + dt(6) * 60 * 1000 + dt(7) * 1000 &
                              + dt(8)
                  end if

                  ! Use rank for low bits and time for high bits uless rank has a lot of bits
                  if(bit_size(my_rank) <= bit_size(t)) then
                        t = my_rank + ishft(t, bit_size(my_rank))
                  else
                        t = ieor(t, int(my_rank,kind(t)))
                  end if

                  ! Here we're using a crappy RNG to seed the better one.
                  do i=1,n
                  seed(i) = lcg(t)
                  end do

            end if
            call random_seed(put=seed)
      contains
            function lcg(s)
                  implicit none
                  integer :: lcg
                  integer(int64) :: s

                  if(s == 0) then
                        s = 104729
                  else
                        s = mod(s, 4294967296_int64)
                  end if
                  s = mod(s*279470273_int64, 4294967291_int64)
                  lcg = int(mod(s, int(huge(0), int64)), kind(0))
            end function lcg
      end subroutine seed_mpi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_beta(Ni_tot_sys,beta)
            use dimensions
            use grid, only: qx,qy,qz
            use inputs, only: nf_init
            implicit none
            integer(4), intent(in):: Ni_tot_sys
            real, intent(out):: beta
            real:: vol
            
            
            vol = ((qx(nx-1)-qx(1))*(qy(ny-1)-qy(1))*(qz(nz-1)-qz(1)))
            beta = (Ni_tot_sys/vol)/nf_init
            
            !write(*,*) 'beta....',beta
            
      end subroutine get_beta
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_gradP()
            use dimensions
            use grid_interp
            use var_arrays, only: np, gradP
            use grid, only: dx_grid, dy_grid, dz_grid
            use inputs, only: etemp0, mion, kboltz
            implicit none
            real:: np1,gdnp,a0,etemp, gnpf(nx,ny,nz,3)
            integer:: i,j,k
            !real, allocatable :: gnpf(:,:,:,:)
            !allocate(gnpf(nx,ny,nz,3))
            etemp = etemp0*11604.505  !eV to Kelvin
            do i=2,nx-1
                  do j=2,ny-1
                        do k=2,nz-1
                              np1 =  0.5*(np(i+1,j,k)+np(i,j,k))
                              gdnp = (np(i+1,j,k)-np(i,j,k))/dx_grid(i)
                              a0 = kboltz*etemp/(mion*np1)
!                             a(i,j,k,1) = a0*gdnp
                              gnpf(i,j,k,1) = a0*gdnp

                              np1 =  0.5*(np(i,j+1,k)+np(i,j,k))
                              gdnp = (np(i,j+1,k)-np(i,j,k))/dy_grid(j)
                              !write(*,*) 'np1',np1,np(i,1,k),np(i,2,k),np(i,3,k),np(i,4,k)
                              a0 = kboltz*etemp/(mion*np1)
!                              a(i,j,k,2) = a0*gdnp
                              gnpf(i,j,k,2) = a0*gdnp

                              np1 =  0.5*(np(i,j,k+1)+np(i,j,k))
                              gdnp = (np(i,j,k+1)-np(i,j,k))/dz_grid(k)
                              a0 = kboltz*etemp/(mion*np1)
!                              a(i,j,k,3) = a0*gdnp
                              gnpf(i,j,k,3) = a0*gdnp
                        enddo
                  enddo
            enddo
            
            call face_to_center(gnpf,gradP)
            !gradP(:,:,:,:) = 0
            
      end subroutine get_gradP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module misc
