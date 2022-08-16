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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_bndry_Eflux(b1,E,bndry_Eflux)
            use dimensions
            use inputs, only: mO,q,dt,dx,dy,km_to_m,mu0
            implicit none
            real, intent(in):: b1(nx,ny,nz,3), E(nx,ny,nz,3)
            real, intent(inout):: bndry_Eflux
            integer:: i,j,k
            real:: exb_flux, mO_q
            
            mO_q = mO/q
            !k=2 face
            do i = 2, nx
                  do j= 2, ny
!                        m=3
                        k=2
                        exb_flux = (mO_q)**2*(1.0/mu0)*dt*dx*dy* &
                              (E(i,j,k,1)*b1(i,j,k,2) - E(i,j,k,2)*b1(i,j,k,1))* &
                              km_to_m**3
                        bndry_Eflux = bndry_Eflux + exb_flux
                        
                  enddo
            enddo
            
            !k=nx face
            
            do i =2,nx
                  do j=2,ny
                        k = nz-1
                        exb_flux = (mO_q)**2*(1.0/mu0)*dt*dx*dy* &
                              (E(i,j,k,1)*b1(i,j,k,2) - E(i,j,k,2)*b1(i,j,k,1))* &
                              km_to_m**3
                        bndry_Eflux = bndry_Eflux - exb_flux
                        
                  enddo
            enddo
      
      end subroutine get_bndry_Eflux
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_beta(Ni_tot_sys,beta)
            use dimensions
            use grid, only: qx,qy,qz
            use inputs, only: np_top
            implicit none
            integer(4), intent(in):: Ni_tot_sys
            real, intent(out):: beta
            real:: vol
            
            
            vol = ((qx(nx-1)-qx(1))*(qy(ny-1)-qy(1))*(qz(nz-1)-qz(1)))
            beta = (Ni_tot_sys/vol)/np_top
            
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
            real:: np1,gdnp,a0,etemp,gnpf(nx,ny,nz,3)
            integer:: i,j,k
            
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
      subroutine get_v_dist()
            use MPI
            use dimensions
            use var_arrays, only: vdist_init,vdist_add, vpp_init,vpp_add,Ni_init, Ni_tot,vp
            implicit none
            integer:: i,j,k,m,l,vxb,vxe,vyb,vye,vzb,vze,ierr,count
            integer, allocatable:: recvbuf(:)
            
!            integer:: vdist_init(-80:80),vdist_add(-80:80), Ni_init, Ni_tot
!            Ni_init=10
!            Ni_tot = 100
            
            vxb=-80
            vxe=80
            vyb=-80
            vye=80
            vzb=-80
            vze=80
            count = (-vxb+vxe+1)*(-vyb+vye+1)
            allocate(recvbuf(count))
            
            vdist_init(:,:) = 0
            vdist_add(:,:) = 0
            vpp_init(:,:) = 0
            vpp_add(:,:) = 0
            
            do l=1,Ni_init
                  i=floor(vp(l,1)-57.0)
                  j=floor(vp(l,2))
                  if ( (i .lt. vxb) .or. (i .gt. vxe) ) then
                        cycle
                  endif
                  if ( (j .lt. vyb) .or. (j .gt. vye) ) then
                        cycle
                  endif
                  vdist_init(i,j) = vdist_init(i,j) + 1
            enddo
            do l= Ni_init+1, Ni_tot
                  i=floor(vp(l,1)-57.0)
                  j=floor(vp(l,2))
                  if ( (i .lt. vxb) .or. (i .gt. vxe) ) then
                        cycle
                  endif
                  if ( (j .lt. vyb) .or. (j .gt. vye) ) then
                        cycle
                  endif
                  vdist_add(i,j) = vdist_add(i,j) + 1
            enddo
            do l= 1, Ni_init
                  m=floor(sqrt((vp(l,1)-57.0)**2+vp(l,2)**2))  ! -57 inside sqrt
                  k=floor(vp(l,3))
                  if ( (m .lt. vxb) .or. (i .gt. vxe) ) then
                        cycle
                  endif
                  if ( (k .lt. vyb) .or. (j .gt. vye) ) then
                        cycle
                  endif
                  vpp_init(m,k) = vpp_init(m,k) + 1
            enddo
            do l= Ni_init+1, Ni_tot
                  m=floor(sqrt((vp(l,1)-57.0)**2+vp(l,2)**2)) ! -57 inside sqrt
                  k=floor(vp(l,3))
                  if ( (m .lt. vxb) .or. (i .gt. vxe) ) then
                        cycle
                  endif
                  if ( (k .lt. vyb) .or. (j .gt. vye) ) then
                        cycle
                  endif
                  vpp_add(m,k) = vpp_add(m,k) + 1
            enddo
            
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(vdist_init(:,:),recvbuf,count,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
            vdist_init(:,:) = reshape(recvbuf,(/(-vxb+vxe+1),(-vyb+vye+1)/))
            
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(vdist_add(:,:),recvbuf,count,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
            vdist_add(:,:) = reshape(recvbuf,(/(-vxb+vxe+1),(-vyb+vye+1)/))
            
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(vpp_init(:,:),recvbuf,count,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
            vpp_init(:,:) = reshape(recvbuf,(/(-vxb+vxe+1),(-vyb+vye+1)/))
            
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(vpp_add(:,:),recvbuf,count,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
            vpp_add(:,:) = reshape(recvbuf,(/(-vxb+vxe+1),(-vyb+vye+1)/))
            deallocate(recvbuf)
            
      end subroutine get_v_dist
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      subroutine get_np3(nfp,np3)
!            use dimensions
!            use boundary
!            implicit none
!            real, intent(in):: nfp(nx,ny,nz)
!            real, intent(out):: np3(nx,ny,nz,3)
!            integer:: i,j,k
            
!            do i = 2, nx-1
!                  do j= 2, ny-1
!                        do k = 2, nz-1
!                              np3(i,j,k,1) = 0.5*(nfp(i,j,k)+nfp(i+1,j,k))
!                              np3(i,j,k,2) = 0.5*(nfp(i,j,k)+nfp(i,j+1,k))
!                              np3(i,j,k,3) = 0.5*(nfp(i,j,k)+nfp(i,j,k+1))
!                        enddo
!                  enddo
!            enddo
            
!            call boundary_vector(np3)
!            call periodic(np3)
            
!      end subroutine get_np3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            
end module misc
