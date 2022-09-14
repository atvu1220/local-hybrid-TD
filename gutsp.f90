module gutsp
      implicit none
      contains

      subroutine remove_ion(ion_l)
! Removes particles from the simulation that have gone out of bounds
            use dimensions
            use inputs, only: km_to_m
            use var_arrays, only: xp,vp,vp1,Ni_tot,input_E,ijkp,beta,beta_p,m_arr,wght,mix_ind,vplus,vminus
            implicit none
            integer, intent(in):: ion_l
            integer:: l,m,tl
            
            do m=1,3    !remove ion energy from total input energy
                  input_E = input_E - 0.5*m_arr(ion_l)*(vp(ion_l,m)*km_to_m)**2 &
                        / (beta * beta_p(ion_l))
            enddo

            do l = ion_l, Ni_tot-1
                  do m = 1,3
                        xp(l,m) = xp(l+1,m)
                        vp(l,m) = vp(l+1,m)
                        vplus(l,m) = vplus(l+1,m)
                        vminus(l,m) = vminus(l+1,m)
                        vp1(l,m) = vp1(l+1,m)
                        ijkp(l,m) = ijkp(l+1,m)
                  enddo
                  beta_p(l) = beta_p(l+1)
                  m_arr(l) = m_arr(l+1)
                  mix_ind(l) = mix_ind(l+1)
                
            enddo

            do m=1,8
                  do l= ion_l,Ni_tot-1
                        wght(l,m) = wght(l+1,m)
                  enddo
            enddo

            Ni_tot = Ni_tot - 1

      end subroutine remove_ion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine extrapol_up()
! This subroutine does the provisional extrapolation of the particle
! bulk flow velocity to time level n, and replaces up_n-3/2 
! with up_n-1/2
            use dimensions
            use var_arrays, only: up,vp,vp1,np,Ni_tot,beta,beta_p,wght
            implicit none
            real:: v_at_n(Ni_max,3)
            integer:: l,m
            
            do l=1,Ni_tot
                  do m=1,3
                        v_at_n(l,m) = 1.5*vp(l,m) -0.5*vp1(l,m)
                  enddo
            enddo
            
            call update_up(v_at_n)
            
      end subroutine extrapol_up
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_Ep()
            use dimensions
            use grid_interp
            use var_arrays, only: Ep,aj,up,btc,Ni_tot,ijkp,wght, gradP,xp,vp,np
            use inputs, only: mion
            implicit none
            real:: ajc(nx,ny,nz,3), &     !aj at cell center
                   aa(3),bb(3),cc(3),aj3(3),up3(3),btc3(3), gradP3(3)    !dummy variables
            integer:: l,i,j,k,m,ip,jp,kp
            
            
            call face_to_center(aj,ajc)
!            call face_to_center(up,upc)
            
            do l=1, Ni_tot
                  i=ijkp(l,1)
                  j=ijkp(l,2)
                  k=ijkp(l,3)
                  
                  ip=i+1
                  jp=j+1
                  kp=k+1            
                  
                  do m=1,3
                        aj3(m) = ajc(i,j,k,m)*wght(l,1) + ajc(ip,j,k,m)*wght(l,2) &
                              + ajc(i,j,kp,m)*wght(l,3) + ajc(ip,j,kp,m)*wght(l,4) &
                              + ajc(i,jp,k,m)*wght(l,5) + ajc(ip,jp,k,m)*wght(l,6) &
                              + ajc(i,jp,kp,m)*wght(l,7) + ajc(ip,jp,kp,m)*wght(l,8)

                        up3(m) = up(i,j,k,m)*wght(l,1) + up(ip,j,k,m)*wght(l,2) &
                              + up(i,j,kp,m)*wght(l,3) + up(ip,j,kp,m)*wght(l,4) &
                              + up(i,jp,k,m)*wght(l,5) + up(ip,jp,k,m)*wght(l,6) &
                              + up(i,jp,kp,m)*wght(l,7) + up(ip,jp,kp,m)*wght(l,8)

!                        uf3(m) = ufc(i,j,k,m)*wght(l,1) + ufc(ip,j,k,m)*wght(l,2) 
!                              + ufc(i,j,kp,m)*wght(l,3) + ufc(ip,j,kp,m)*wght(l,4) &
!                              + ufc(i,jp,k,m)*wght(l,5) + ufc(ip,jp,k,m)*wght(l,6) &
!                              + ufc(i,jp,kp,m)*wght(l,7) + ufc(ip,jp,kp,m)*wght(l,8)

                        btc3(m) = btc(i,j,k,m)*wght(l,1) & 
                              + btc(ip,j,k,m)*wght(l,2) &
                              + btc(i,j,kp,m)*wght(l,3) &
                              + btc(ip,j,kp,m)*wght(l,4) &
                              + btc(i,jp,k,m)*wght(l,5) &
                              + btc(ip,jp,k,m)*wght(l,6) &
                              + btc(i,jp,kp,m)*wght(l,7) &
                              + btc(ip,jp,kp,m)*wght(l,8)
                               
                               
                       !electron pressure term
                       gradP3(m) = gradP(i,j,k,m)*wght(l,1) &
                               + gradP(ip,j,k,m)*wght(l,2) &
                               + gradP(i,j,kp,m)*wght(l,3) &
                               + gradP(ip,j,kp,m)*wght(l,4) &
                               + gradP(i,jp,k,m)*wght(l,5) &
                               + gradP(ip,jp,k,m)*wght(l,6) &
                               + gradP(i,jp,kp,m)*wght(l,7) &
                               + gradP(ip,jp,kp,m)*wght(l,8) 
                   
                  enddo 
                  do m=1,3 
                        aa(m) = aj3(m) - up3(m)
                        bb(m) = btc3(m)
                  enddo
                  
                  !Cross product
                  cc(1) = aa(2)*bb(3) - aa(3)*bb(2)
                  cc(2) = aa(3)*bb(1) - aa(1)*bb(3)
                  cc(3) = aa(1)*bb(2) - aa(2)*bb(1)
                  
                  
                  do m=1,3
                        Ep(l,m) = cc(m) - gradP3(m) !add in electron pressure term
                  enddo
                 
                  
            enddo
            
      end subroutine get_Ep
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_vplus_vminus()
            use dimensions
            use inputs, only: dt
            use var_arrays, only: Ep,btc,vp,vplus,vminus,Ni_tot,ijkp,wght
            implicit none
            real:: a1,a2,a3,a_d,B2,dt2,Bx,By,Bz,vminus_x_B(3),vminus_dot_B,btc3(3)
            integer:: l,i,j,k,ip,jp,kp,m
            
            do l=1, Ni_tot
                  do m=1,3
                        vminus(l,m) = vp(l,m) + 0.5*dt*Ep(l,m)
                  enddo
            enddo
            
            do l = 1, Ni_tot
                  i=ijkp(l,1)
                  j=ijkp(l,2)
                  k=ijkp(l,3)
                  
                  ip=i+1
                  jp=j+1
                  kp=k+1
                  
                  do m=1,3
                        btc3(m) = btc(i,j,k,m)*wght(l,1) &
                              + btc(ip,j,k,m)*wght(l,2) & 
                              + btc(i,j,kp,m)*wght(l,3) &
                              + btc(ip,j,kp,m)*wght(l,4) &
                              + btc(i,jp,k,m)*wght(l,5) &
                              + btc(ip,jp,k,m)*wght(l,6) &
                              + btc(i,jp,kp,m)*wght(l,7) &
                              + btc(ip,jp,kp,m)*wght(l,8) 
                  enddo
                  
                  vminus_x_B(1) = vminus(l,2)*btc3(3) - &
                        vminus(l,3)*btc3(2) 
                  vminus_x_B(2) = vminus(l,3)*btc3(1) - &
                        vminus(l,1)*btc3(3)   
                  vminus_x_B(3) = vminus(l,1)*btc3(2) - &
                        vminus(l,2)*btc3(1)  

                  vminus_dot_B = vminus(l,1)*btc3(1) + &
                        vminus(l,2)*btc3(2) + &
                        vminus(l,3)*btc3(3) 

                  Bx = btc3(1)
                  By = btc3(2)
                  Bz = btc3(3)
      
                  B2 = Bx*Bx + By*By + Bz*Bz
                  dt2 = dt*dt

                  a_d = 1.0/(1.0 + (0.25*B2*dt2))
                  a1 = (1 - (0.25*B2*dt2))*a_d
                  a2 = dt*a_d
                  a3 = 0.5*dt2*a_d
                  
                  do m=1,3
                        vplus(l,m) = a1*vminus(l,m) + a2*vminus_x_B(m) + &
                              a3*vminus_dot_B*btc3(m)
                  enddo
            enddo
            
      end subroutine get_vplus_vminus
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine improve_up()
! The routine calculates v at time level n, and the associated bulk
! flow velocity up using the v+, v- technique.  The new up at
! time level n replaces the provisional extrapolation for up.
            use dimensions
            use var_arrays, only:vp1,vplus,vminus,up,np,Ni_tot,beta,beta_p,wght
            implicit none
            integer:: l,m
            
            do l=1, Ni_tot
                  do m=1,3
                        vp1(l,m) = 0.5*(vplus(l,m) + vminus(l,m))
                  enddo
            enddo
            
            call update_up(vp1)
            
      end subroutine improve_up
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_vp_final()
            use dimensions
            use inputs, only: dt
            use var_arrays, only: Ep,vp,vp1,vplus,Ni_tot
            implicit none

            integer:: l,m,ierr
            
            do l=1, Ni_tot
                  do m=1,3
                        vp1(l,m) = vp(l,m)      !to be used in extrapol_up for n-3/2
                        vp(l,m) = vplus(l,m) + 0.5*dt*Ep(l,m)
                  enddo
            enddo
            
      end subroutine get_vp_final
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine move_ion_half()
            use dimensions
            use boundary
            use inputs, only: dt,boundx
            use grid, only: qx,qy,qz
            use var_arrays, only: xp,vp,Ni_tot
            implicit none
            real:: dth
            integer:: l
            
            dth = dt/2.0
            do l=1,Ni_tot
                  xp(l,1) = xp(l,1) + dth*vp(l,1)
                  xp(l,2) = xp(l,2) + dth*vp(l,2)
                  xp(l,3) = xp(l,3) + dth*vp(l,3)
            enddo      
            call particle_boundary()      
      
      end subroutine move_ion_half

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_interp_weights()
! Weights are used for trilinear interpolation to/from main cell
! centers to particle positions.  For each particle there are 8
! grid points associated with the interpolation.  These 8 points
! are determined by the location of the particle within the main
! cell.  There are 8 sets of 8 grid points for each cell.
            use dimensions
            use grid, only: qx,qy,qz
            use inputs, only: dx,dy
            use var_arrays, only: xp,vp,Ni_tot,ijkp,wght,Ep
            implicit none
            real:: vol,x1,x2,y1,y2,z1,z2
            integer:: i,j,k,l

            do l=1, Ni_tot

                  call get_pindex(l)
                  i=ijkp(l,1)
                  j=ijkp(l,2)
                  k=ijkp(l,3)

                  vol = 1.0/((qx(i+1)-qx(i))*(qy(j+1)-qy(j))*(qz(k+1)-qz(k)))

                  x1=abs(xp(l,1)-qx(i))
                  x2=abs(xp(l,1)-qx(i+1))
                  y1=abs(xp(l,2)-qy(j))
                  y2=abs(xp(l,2)-qy(j+1))
                  z1=abs(xp(l,3)-qz(k))
                  z2=abs(xp(l,3)-qz(k+1))
                 
                  wght(l,1) = x2*y2*z2*vol
                  wght(l,2) = x1*y2*z2*vol
                  wght(l,3) = x2*y2*z1*vol
                  wght(l,4) = x1*y2*z1*vol
                  wght(l,5) = x2*y1*z2*vol
                  wght(l,6) = x1*y1*z2*vol
                  wght(l,7) = x2*y1*z1*vol
                  wght(l,8) = x1*y1*z1*vol
                
            enddo
            
            
      end subroutine get_interp_weights
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine update_np()
! Weight density to eight nearest grid points
            use dimensions
            use MPI
            use grid, only: dx_grid,dy_grid,dz_grid
            use boundary
            use var_arrays, only: np,Ni_tot,ijkp,beta,beta_p,wght
            implicit none
            real:: volb, recvbuf(nx*ny*nz)
            integer:: i,j,k,l,ip,jp,kp,ierr,count
            
            count = nx*ny*nz
            
            do i=1,nx
                  do j=1,ny
                        do k=1,nz
                              np(i,j,k) = 0.0
                        enddo
                  enddo
            enddo
            
            do l=1, Ni_tot
                  i=ijkp(l,1)
                  j=ijkp(l,2)
                  k=ijkp(l,3)
                  
                  ip=i+1
                  jp=j+1
                  kp=k+1
                  
                  volb = 1.0/(dx_grid(i)*dy_grid(j)*dz_grid(k)*beta*beta_p(l))
                  
                  np(i,j,k) = np(i,j,k) + wght(l,1)*volb
                  np(ip,j,k) = np(ip,j,k) + wght(l,2)*volb
                  np(i,j,kp) = np(i,j,kp) + wght(l,3)*volb
                  np(ip,j,kp) = np(ip,j,kp) + wght(l,4)*volb
                  np(i,jp,k) = np(i,jp,k) + wght(l,5)*volb
                  np(ip,jp,k) = np(ip,jp,k) + wght(l,6)*volb
                  np(i,jp,kp) = np(i,jp,kp) + wght(l,7)*volb
                  np(ip,jp,kp) = np(ip,jp,kp) + wght(l,8)*volb
                  
                  
            enddo
                  
            !Used for periodic boundary conditions
            call add_boundary_scalar(np)
            
            
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(np(:,:,:),recvbuf,count,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
            
            np(:,:,:) = reshape(recvbuf,(/nx,ny,nz/))
            
            call boundary_scalar(np)
            
      end subroutine update_np
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      subroutine update_rho()
! Weight density to eight neares grid points
            use dimensions
            use MPI
            use boundary
            use inputs, only: mion
            use grid, only: qx,qy,qz
            use var_arrays, only: mnp,Ni_tot,ijkp,beta,beta_p,wght
            implicit none
            real:: volb, recvbuf(nx*ny*nz)
            integer:: i,j,k,l,ip,jp,kp,count,ierr
            
            count = nx*ny*nz
            
            do i=1,nx
                  do j=1,ny
                        do k=1,nz
                              mnp(i,j,k) = 0.0
                        enddo
                  enddo
            enddo
            
            do l = 1, Ni_tot
                  i=ijkp(l,1)
                  j=ijkp(l,2)
                  k=ijkp(l,3)
                  
                  ip=i+1
                  jp=j+1
                  kp=k+1
                  
                  volb = (qx(ip)-qx(i))*(qy(jp)-qy(j))*(qz(kp)-qz(k))*beta*beta_p(l)
                  
                  mnp(i,j,k) = mnp(i,j,k) + (wght(l,1))/volb
                  mnp(ip,j,k) = mnp(ip,j,k) + (wght(l,2))/volb
                  mnp(i,j,kp) = mnp(i,j,kp) + (wght(l,3))/volb
                  mnp(ip,j,kp) = mnp(ip,j,kp) + (wght(l,4))/volb
                  mnp(i,jp,k) = mnp(i,jp,k) + (wght(l,5))/volb
                  mnp(ip,jp,k) = mnp(ip,jp,k) + (wght(l,6))/volb
                  mnp(i,jp,kp) = mnp(i,jp,kp) + (wght(l,7))/volb
                  mnp(ip,jp,kp) = mnp(ip,jp,kp) + (wght(l,8))/volb
                  
            enddo
            
            mnp(:,:,:) = mion*mnp(:,:,:) !mass density
            
            !Used for periodic boundary conditions
            call add_boundary_scalar(mnp)
          
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(mnp(:,:,:),recvbuf,count,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
            
            mnp(:,:,:) = reshape(recvbuf,(/nx,ny,nz/))

            call boundary_scalar(mnp)
            
      end subroutine update_rho
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine update_mixed()
! Weight density to eight nearest grid points
            use dimensions
            use MPI
            use var_arrays, only: Ni_tot,ijkp,mix_ind,mixed
            use boundary
            implicit none
            real:: recvbuf(nx*ny*nz)
            real:: mix_cnt(nx,ny,nz)
            integer:: i,j,k,l,count,ierr
            
            count = nx*ny*nz
            
            do i=1,nx
                  do j=1,ny
                        do k=1,nz
                              mixed(i,j,k) = 0.0
                              mix_cnt(i,j,k) = 0.0
                        enddo
                  enddo
            enddo
            
            do l = 1, Ni_tot
                  i=ijkp(l,1)
                  j=ijkp(l,2)
                  k=ijkp(l,3)
                  
                  mixed(i,j,k) = mixed(i,j,k) + mix_ind(l)
                  mix_cnt(i,j,k) = mix_cnt(i,j,k) + 1.0
                  
            enddo
            
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(mixed(:,:,:),recvbuf,count,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
            
            mixed(:,:,:) = reshape(recvbuf,(/nx,ny,nz/))
            
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(mix_cnt(:,:,:),recvbuf,count,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
            
            mix_cnt(:,:,:) = reshape(recvbuf,(/nx,ny,nz/))
            
            where(mix_cnt(:,:,:) .gt. 0.0)
            	mixed(:,:,:) = mixed(:,:,:)/mix_cnt(:,:,:)
            endwhere
            
            call add_boundary_scalar(mixed)
            
            call boundary_scalar(mixed)
            
      end subroutine update_mixed
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine update_up(vp)
            use dimensions
            use MPI
            use boundary
            use grid, only: qx,qy,qz
            use var_arrays, only: np,up,Ni_tot,ijkp,beta,beta_p,wght
            implicit none
            real, intent(in):: vp(Ni_max,3)
            real:: ct(nx,ny,nz,3), recvbuf(nx*ny*nz*3),volb,nvolb
            integer:: i,j,k,m,l,ip,jp,kp,count,ierr
            
            count=nx*ny*nz*3
            
            do i=1,nx
                  do j=1,ny
                        do k=1,nz
                              do m=1,3
                                    up(i,j,k,m) = 0.0
                                    ct(i,j,k,m) = 0.0
                              enddo
                        enddo
                  enddo
            enddo
            
            do l=1,Ni_tot
                  i=ijkp(l,1)
                  j=ijkp(l,2)
                  k=ijkp(l,3)
                  
                  ip=i+1
                  jp=j+1
                  kp=k+1
                  
                  volb = (qx(ip)-qx(i))*(qy(jp)-qy(j))*(qz(kp)-qz(k))*beta*beta_p(l)
                  
!                  if(np(i,j,k) .gt. 0.0) then                  
                  nvolb = 1.0/(np(i,j,k)*volb)
                  ct(i,j,k,1) = ct(i,j,k,1) + vp(l,1)*wght(l,1)*nvolb
                  ct(i,j,k,2) = ct(i,j,k,2) + vp(l,2)*wght(l,1)*nvolb
                  ct(i,j,k,3) = ct(i,j,k,3) + vp(l,3)*wght(l,1)*nvolb
!                  endif

!                  if (np(ip,j,k) .gt. 0.0) then
                  nvolb = 1.0/(np(ip,j,k)*volb)
                  ct(ip,j,k,1) = ct(ip,j,k,1) + vp(l,1)*wght(l,2)*nvolb
                  ct(ip,j,k,2) = ct(ip,j,k,2) + vp(l,2)*wght(l,2)*nvolb
                  ct(ip,j,k,3) = ct(ip,j,k,3) + vp(l,3)*wght(l,2)*nvolb
!                  endif

!                  if (np(i,j,kp) .gt. 0.0) then
                  nvolb = 1.0/(np(i,j,kp)*volb)
                  ct(i,j,kp,1) = ct(i,j,kp,1) + vp(l,1)*wght(l,3)*nvolb
                  ct(i,j,kp,2) = ct(i,j,kp,2) + vp(l,2)*wght(l,3)*nvolb
                  ct(i,j,kp,3) = ct(i,j,kp,3) + vp(l,3)*wght(l,3)*nvolb
!                  endif

!                  if (np(ip,j,kp) .gt. 0.0) then
                  nvolb = 1.0/(np(ip,j,kp)*volb)
                  ct(ip,j,kp,1) = ct(ip,j,kp,1) + vp(l,1)*wght(l,4)*nvolb
                  ct(ip,j,kp,2) = ct(ip,j,kp,2) + vp(l,2)*wght(l,4)*nvolb
                  ct(ip,j,kp,3) = ct(ip,j,kp,3) + vp(l,3)*wght(l,4)*nvolb
!                  endif

!                  if (np(i,jp,k) .gt. 0.0) then
                  nvolb = 1.0/(np(i,jp,k)*volb)
                  ct(i,jp,k,1) = ct(i,jp,k,1) + vp(l,1)*wght(l,5)*nvolb
                  ct(i,jp,k,2) = ct(i,jp,k,2) + vp(l,2)*wght(l,5)*nvolb
                  ct(i,jp,k,3) = ct(i,jp,k,3) + vp(l,3)*wght(l,5)*nvolb
!                  endif

!                  if (np(ip,jp,k) .gt. 0.0) then
                  nvolb = 1.0/(np(ip,jp,k)*volb)
                  ct(ip,jp,k,1) = ct(ip,jp,k,1) + vp(l,1)*wght(l,6)*nvolb
                  ct(ip,jp,k,2) = ct(ip,jp,k,2) + vp(l,2)*wght(l,6)*nvolb
                  ct(ip,jp,k,3) = ct(ip,jp,k,3) + vp(l,3)*wght(l,6)*nvolb
!                  endif

!                  if (np(i,jp,kp) .gt. 0.0) then
                  nvolb = 1.0/(np(i,jp,kp)*volb)
                  ct(i,jp,kp,1) = ct(i,jp,kp,1) + vp(l,1)*wght(l,7)*nvolb
                  ct(i,jp,kp,2) = ct(i,jp,kp,2) + vp(l,2)*wght(l,7)*nvolb
                  ct(i,jp,kp,3) = ct(i,jp,kp,3) + vp(l,3)*wght(l,7)*nvolb
!                  endif

!                  if (np(ip,jp,kp) .gt. 0.0) then
                  nvolb = 1.0/(np(ip,jp,kp)*volb)
                  ct(ip,jp,kp,1) = ct(ip,jp,kp,1) + vp(l,1)*wght(l,8)*nvolb
                  ct(ip,jp,kp,2) = ct(ip,jp,kp,2) + vp(l,2)*wght(l,8)*nvolb
                  ct(ip,jp,kp,3) = ct(ip,jp,kp,3) + vp(l,3)*wght(l,8)*nvolb
!                  endif

            enddo
            
            !Used for periodic boundary conditions
            call add_boundary_vector(ct)
            
            call boundary_vector(ct)
            
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(ct(:,:,:,:),recvbuf,count,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
            
            ct(:,:,:,:) = reshape(recvbuf,(/nx,ny,nz,3/))
            up=ct
            
            call boundary_vector(up)      
            
      end subroutine update_up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine update_up_cold(vp)
            use dimensions
            use MPI
            use boundary
            use grid, only: qx,qy,qz
            use var_arrays, only: np_cold,up_cold,Ni_tot,ijkp,beta,beta_p,wght,mix_ind
            implicit none
            real, intent(in):: vp(Ni_max,3)
            real:: ct(nx,ny,nz,3), recvbuf(nx*ny*nz*3),volb,nvolb
            integer:: i,j,k,m,l,ip,jp,kp,count,ierr
            
            count=nx*ny*nz*3

            do i=1,nx
                  do j=1,ny
                        do k=1,nz
                              do m=1,3
                                    up_cold(i,j,k,m) = 0.0
                                    ct(i,j,k,m) = 0.0
                              enddo
                        enddo
                  enddo
            enddo
            
            l=1
            do while (l .le. Ni_tot) 
            	do while ( (mix_ind(l) .ne. 0) .and. (l .lt. Ni_tot) )
			      l=l+1
	            end do
		      if  (l .ge. Ni_tot) then
			      EXIT
		      endif
            
                  !do l=1,Ni_tot
                  i=ijkp(l,1)
                  j=ijkp(l,2)
                  k=ijkp(l,3)
                  
                  ip=i+1
                  jp=j+1
                  kp=k+1
                  
                  volb = (qx(ip)-qx(i))*(qy(jp)-qy(j))*(qz(kp)-qz(k))*beta*beta_p(l)
                  
                  if(np_cold(i,j,k) .gt. 0.0) then                  
                        nvolb = 1.0/(np_cold(i,j,k)*volb)
                        ct(i,j,k,1) = ct(i,j,k,1) + vp(l,1)*wght(l,1)*nvolb
                        ct(i,j,k,2) = ct(i,j,k,2) + vp(l,2)*wght(l,1)*nvolb
                        ct(i,j,k,3) = ct(i,j,k,3) + vp(l,3)*wght(l,1)*nvolb
                  endif

                  if (np_cold(ip,j,k) .gt. 0.0) then
                        nvolb = 1.0/(np_cold(ip,j,k)*volb)
                        ct(ip,j,k,1) = ct(ip,j,k,1) + vp(l,1)*wght(l,2)*nvolb
                        ct(ip,j,k,2) = ct(ip,j,k,2) + vp(l,2)*wght(l,2)*nvolb
                        ct(ip,j,k,3) = ct(ip,j,k,3) + vp(l,3)*wght(l,2)*nvolb
                  endif

                  if (np_cold(i,j,kp) .gt. 0.0) then
                        nvolb = 1.0/(np_cold(i,j,kp)*volb)
                        ct(i,j,kp,1) = ct(i,j,kp,1) + vp(l,1)*wght(l,3)*nvolb
                        ct(i,j,kp,2) = ct(i,j,kp,2) + vp(l,2)*wght(l,3)*nvolb
                        ct(i,j,kp,3) = ct(i,j,kp,3) + vp(l,3)*wght(l,3)*nvolb
                  endif

                  if (np_cold(ip,j,kp) .gt. 0.0) then
                        nvolb = 1.0/(np_cold(ip,j,kp)*volb)
                        ct(ip,j,kp,1) = ct(ip,j,kp,1) + vp(l,1)*wght(l,4)*nvolb
                        ct(ip,j,kp,2) = ct(ip,j,kp,2) + vp(l,2)*wght(l,4)*nvolb
                        ct(ip,j,kp,3) = ct(ip,j,kp,3) + vp(l,3)*wght(l,4)*nvolb
                  endif

                  if (np_cold(i,jp,k) .gt. 0.0) then
                        nvolb = 1.0/(np_cold(i,jp,k)*volb)
                        ct(i,jp,k,1) = ct(i,jp,k,1) + vp(l,1)*wght(l,5)*nvolb
                        ct(i,jp,k,2) = ct(i,jp,k,2) + vp(l,2)*wght(l,5)*nvolb
                        ct(i,jp,k,3) = ct(i,jp,k,3) + vp(l,3)*wght(l,5)*nvolb
                  endif

                  if (np_cold(ip,jp,k) .gt. 0.0) then
                        nvolb = 1.0/(np_cold(ip,jp,k)*volb)
                        ct(ip,jp,k,1) = ct(ip,jp,k,1) + vp(l,1)*wght(l,6)*nvolb
                        ct(ip,jp,k,2) = ct(ip,jp,k,2) + vp(l,2)*wght(l,6)*nvolb
                        ct(ip,jp,k,3) = ct(ip,jp,k,3) + vp(l,3)*wght(l,6)*nvolb
                  endif

                  if (np_cold(i,jp,kp) .gt. 0.0) then
                        nvolb = 1.0/(np_cold(i,jp,kp)*volb)
                        ct(i,jp,kp,1) = ct(i,jp,kp,1) + vp(l,1)*wght(l,7)*nvolb
                        ct(i,jp,kp,2) = ct(i,jp,kp,2) + vp(l,2)*wght(l,7)*nvolb
                        ct(i,jp,kp,3) = ct(i,jp,kp,3) + vp(l,3)*wght(l,7)*nvolb
                  endif

                  if (np_cold(ip,jp,kp) .gt. 0.0) then
                        nvolb = 1.0/(np_cold(ip,jp,kp)*volb)
                        ct(ip,jp,kp,1) = ct(ip,jp,kp,1) + vp(l,1)*wght(l,8)*nvolb
                        ct(ip,jp,kp,2) = ct(ip,jp,kp,2) + vp(l,2)*wght(l,8)*nvolb
                        ct(ip,jp,kp,3) = ct(ip,jp,kp,3) + vp(l,3)*wght(l,8)*nvolb
                  endif

	    	      l=l+1
            enddo

            !Used for periodic boundary conditions
            call add_boundary_vector(ct)
            
            call boundary_vector(ct)
            
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(ct(:,:,:,:),recvbuf,count,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
            
            ct(:,:,:,:) = reshape(recvbuf,(/nx,ny,nz,3/))
            up_cold=ct

            call boundary_vector(up_cold)      
            
      end subroutine update_up_cold
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine update_up_mixed(vp)
            use dimensions
            use MPI
            use boundary
            use grid, only: qx,qy,qz
            use var_arrays, only: np_mixed,up_mixed,Ni_tot,ijkp,beta,beta_p,wght,mix_ind
            implicit none
            real, intent(in):: vp(Ni_max,3)
            real:: ct(nx,ny,nz,3), recvbuf(nx*ny*nz*3),volb,nvolb
            integer:: i,j,k,m,l,ip,jp,kp,count,ierr
            
            count=nx*ny*nz*3
            
            do i=1,nx
                  do j=1,ny
                        do k=1,nz
                              do m=1,3
                                    up_mixed(i,j,k,m) = 0.0
                                    ct(i,j,k,m) = 0.0
                              enddo
                        enddo
                  enddo
            enddo
            
            l=1
            do while (l .le. Ni_tot) 
            	do while ( (mix_ind(l) .eq. 0) .and. (l .lt. Ni_tot) )
			      l=l+1
		      end do
		      if  (l .ge. Ni_tot) then
			      EXIT
		      endif
            
            !do l=1,Ni_tot
                  i=ijkp(l,1)
                  j=ijkp(l,2)
                  k=ijkp(l,3)
                  
                  ip=i+1
                  jp=j+1
                  kp=k+1
                  
                  volb = (qx(ip)-qx(i))*(qy(jp)-qy(j))*(qz(kp)-qz(k))*beta*beta_p(l)
                  
                  if(np_mixed(i,j,k) .gt. 0.0) then                  
                        nvolb = 1.0/(np_mixed(i,j,k)*volb)
                        ct(i,j,k,1) = ct(i,j,k,1) + vp(l,1)*wght(l,1)*nvolb
                        ct(i,j,k,2) = ct(i,j,k,2) + vp(l,2)*wght(l,1)*nvolb
                        ct(i,j,k,3) = ct(i,j,k,3) + vp(l,3)*wght(l,1)*nvolb
                  endif

                  if (np_mixed(ip,j,k) .gt. 0.0) then
                        nvolb = 1.0/(np_mixed(ip,j,k)*volb)
                        ct(ip,j,k,1) = ct(ip,j,k,1) + vp(l,1)*wght(l,2)*nvolb
                        ct(ip,j,k,2) = ct(ip,j,k,2) + vp(l,2)*wght(l,2)*nvolb
                        ct(ip,j,k,3) = ct(ip,j,k,3) + vp(l,3)*wght(l,2)*nvolb
                  endif

                 if (np_mixed(i,j,kp) .gt. 0.0) then
                        nvolb = 1.0/(np_mixed(i,j,kp)*volb)
                        ct(i,j,kp,1) = ct(i,j,kp,1) + vp(l,1)*wght(l,3)*nvolb
                        ct(i,j,kp,2) = ct(i,j,kp,2) + vp(l,2)*wght(l,3)*nvolb
                        ct(i,j,kp,3) = ct(i,j,kp,3) + vp(l,3)*wght(l,3)*nvolb
                  endif

                  if (np_mixed(ip,j,kp) .gt. 0.0) then
                        nvolb = 1.0/(np_mixed(ip,j,kp)*volb)
                        ct(ip,j,kp,1) = ct(ip,j,kp,1) + vp(l,1)*wght(l,4)*nvolb
                        ct(ip,j,kp,2) = ct(ip,j,kp,2) + vp(l,2)*wght(l,4)*nvolb
                        ct(ip,j,kp,3) = ct(ip,j,kp,3) + vp(l,3)*wght(l,4)*nvolb
                  endif

                  if (np_mixed(i,jp,k) .gt. 0.0) then
                        nvolb = 1.0/(np_mixed(i,jp,k)*volb)
                        ct(i,jp,k,1) = ct(i,jp,k,1) + vp(l,1)*wght(l,5)*nvolb
                        ct(i,jp,k,2) = ct(i,jp,k,2) + vp(l,2)*wght(l,5)*nvolb
                        ct(i,jp,k,3) = ct(i,jp,k,3) + vp(l,3)*wght(l,5)*nvolb
                  endif

                  if (np_mixed(ip,jp,k) .gt. 0.0) then
                        nvolb = 1.0/(np_mixed(ip,jp,k)*volb)
                        ct(ip,jp,k,1) = ct(ip,jp,k,1) + vp(l,1)*wght(l,6)*nvolb
                        ct(ip,jp,k,2) = ct(ip,jp,k,2) + vp(l,2)*wght(l,6)*nvolb
                        ct(ip,jp,k,3) = ct(ip,jp,k,3) + vp(l,3)*wght(l,6)*nvolb
                  endif

                  if (np_mixed(i,jp,kp) .gt. 0.0) then
                        nvolb = 1.0/(np_mixed(i,jp,kp)*volb)
                        ct(i,jp,kp,1) = ct(i,jp,kp,1) + vp(l,1)*wght(l,7)*nvolb
                        ct(i,jp,kp,2) = ct(i,jp,kp,2) + vp(l,2)*wght(l,7)*nvolb
                        ct(i,jp,kp,3) = ct(i,jp,kp,3) + vp(l,3)*wght(l,7)*nvolb
                  endif

                  if (np_mixed(ip,jp,kp) .gt. 0.0) then
                        nvolb = 1.0/(np_mixed(ip,jp,kp)*volb)
                        ct(ip,jp,kp,1) = ct(ip,jp,kp,1) + vp(l,1)*wght(l,8)*nvolb
                        ct(ip,jp,kp,2) = ct(ip,jp,kp,2) + vp(l,2)*wght(l,8)*nvolb
                        ct(ip,jp,kp,3) = ct(ip,jp,kp,3) + vp(l,3)*wght(l,8)*nvolb
                  endif

		      l=l+1

            enddo
            
            !Used for periodic boundary conditions
            call add_boundary_vector(ct)
              
            call boundary_vector(ct)
            
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(ct(:,:,:,:),recvbuf,count,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
            
            ct(:,:,:,:) = reshape(recvbuf,(/nx,ny,nz,3/))
            up_mixed=ct

            call boundary_vector(up_mixed)      
            
      end subroutine update_up_mixed
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_temperature()
            use dimensions
            use MPI
            use boundary
            use inputs, only: mion
            use grid, only: qx,qy,qz
            use var_arrays, only: vp,np,temp_p,Ni_tot,ijkp,beta,beta_p,wght,tp
            implicit none
            real:: recvbuf(nx*ny*nz*3),up2(nx,ny,nz,3),up_ave(nx,ny,nz,3),ct(nx,ny,nz,3),volb,nvolb,mvp(Ni_max,3)
            integer:: i,j,k,l,m,ip,jp,kp,count,ierr
            
            count = nx*ny*nz*3
            
            up2(:,:,:,:) = 0.0
            up_ave(:,:,:,:) = 0.0
            ct(:,:,:,:) = 0.0
            
            do m=1,3
                  mvp(1:Ni_tot,m) = vp(1:Ni_tot,m)
            enddo
            
            do l=1,Ni_tot
                  i=ijkp(l,1)
                  j=ijkp(l,2)
                  k=ijkp(l,3)
                  
                  ip=i+1
                  jp=j+1
                  kp=k+1
                  
                  volb = (qx(ip)-qx(i))*(qy(jp)-qy(j))*(qz(kp)-qz(k))*beta*beta_p(l)
                  
                  if (np(i,j,k) .gt. 0.0) then
                        nvolb = np(i,j,k)*volb
                        ct(i,j,k,1) = ct(i,j,k,1) + mvp(l,1)**2*wght(l,1)/nvolb  
                        ct(i,j,k,2) = ct(i,j,k,2) + mvp(l,2)**2*wght(l,1)/nvolb
                        ct(i,j,k,3) = ct(i,j,k,3) + mvp(l,3)**2*wght(l,1)/nvolb
                  endif

                  if (np(ip,j,k) .gt. 0.0) then
                        nvolb = np(ip,j,k)*volb
                        ct(ip,j,k,1) = ct(ip,j,k,1) + mvp(l,1)**2*wght(l,2)/nvolb
                        ct(ip,j,k,2) = ct(ip,j,k,2) + mvp(l,2)**2*wght(l,2)/nvolb
                        ct(ip,j,k,3) = ct(ip,j,k,3) + mvp(l,3)**2*wght(l,2)/nvolb
                  endif

                  if (np(i,j,kp) .gt. 0.0) then
                        nvolb = np(i,j,kp)*volb
                        ct(i,j,kp,1) = ct(i,j,kp,1) + mvp(l,1)**2*wght(l,3)/nvolb
                        ct(i,j,kp,2) = ct(i,j,kp,2) + mvp(l,2)**2*wght(l,3)/nvolb
                        ct(i,j,kp,3) = ct(i,j,kp,3) + mvp(l,3)**2*wght(l,3)/nvolb
                  endif

                  if (np(ip,j,kp) .gt. 0.0) then
                        nvolb = np(ip,j,kp)*volb
                        ct(ip,j,kp,1) = ct(ip,j,kp,1) + mvp(l,1)**2*wght(l,4)/nvolb
                        ct(ip,j,kp,2) = ct(ip,j,kp,2) + mvp(l,2)**2*wght(l,4)/nvolb
                        ct(ip,j,kp,3) = ct(ip,j,kp,3) + mvp(l,3)**2*wght(l,4)/nvolb
                  endif

                  if (np(i,jp,k) .gt. 0.0) then
                        nvolb = np(i,jp,k)*volb
                        ct(i,jp,k,1) = ct(i,jp,k,1) + mvp(l,1)**2*wght(l,5)/nvolb
                        ct(i,jp,k,2) = ct(i,jp,k,2) + mvp(l,2)**2*wght(l,5)/nvolb
                        ct(i,jp,k,3) = ct(i,jp,k,3) + mvp(l,3)**2*wght(l,5)/nvolb
                  endif

                  if (np(ip,jp,k) .gt. 0.0) then
                        nvolb = np(ip,jp,k)*volb
                        ct(ip,jp,k,1) = ct(ip,jp,k,1) + mvp(l,1)**2*wght(l,6)/nvolb
                        ct(ip,jp,k,2) = ct(ip,jp,k,2) + mvp(l,2)**2*wght(l,6)/nvolb
                        ct(ip,jp,k,3) = ct(ip,jp,k,3) + mvp(l,3)**2*wght(l,6)/nvolb
                  endif

                  if (np(i,jp,kp) .gt. 0.0) then
                        nvolb = np(i,jp,kp)*volb
                        ct(i,jp,kp,1) = ct(i,jp,kp,1) + mvp(l,1)**2*wght(l,7)/nvolb
                        ct(i,jp,kp,2) = ct(i,jp,kp,2) + mvp(l,2)**2*wght(l,7)/nvolb
                        ct(i,jp,kp,3) = ct(i,jp,kp,3) + mvp(l,3)**2*wght(l,7)/nvolb
                  endif

                  if (np(ip,jp,kp) .gt. 0.0) then
                        nvolb = np(ip,jp,kp)*volb
                        ct(ip,jp,kp,1) = ct(ip,jp,kp,1) + mvp(l,1)**2*wght(l,8)/nvolb
                        ct(ip,jp,kp,2) = ct(ip,jp,kp,2) + mvp(l,2)**2*wght(l,8)/nvolb
                        ct(ip,jp,kp,3) = ct(ip,jp,kp,3) + mvp(l,3)**2*wght(l,8)/nvolb
                  endif
         
            enddo
            
            !Used for periodic boundary conditions
            call add_boundary_vector(ct)

            call boundary_vector(ct)
                        
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(ct(:,:,:,:),recvbuf,count,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
            
            ct(:,:,:,:) = reshape(recvbuf,(/nx,ny,nz,3/))
            
            do i=1,nx-1
                  do j=1,ny-1
                        do k=1,nz-1
                              up2(i,j,k,1) = 0.5*(ct(i,j,k,1)+ct(i+1,j,k,1))
                              up2(i,j,k,2) = 0.5*(ct(i,j,k,2)+ct(i,j+1,k,2))
                              up2(i,j,k,3) = 0.5*(ct(i,j,k,3)+ct(i,j,k+1,3))
                        enddo
                  enddo
            enddo
            
            call boundary_vector(up2)
            
            ct(:,:,:,:) = 0.0
            
            do l=1,Ni_tot
                  i=ijkp(l,1)
                  j=ijkp(l,2)
                  k=ijkp(l,3)
                  
                  ip=i+1
                  jp=j+1
                  kp=k+1
                  
                  volb = (qx(ip)-qx(i))*(qy(jp)-qy(j))*(qz(kp)-qz(k))*beta*beta_p(l)
                  
                  if (np(i,j,k) .gt. 0.0) then
                        nvolb = np(i,j,k)*volb
                        ct(i,j,k,1) = ct(i,j,k,1) + mvp(l,1)*wght(l,1)/nvolb  
                        ct(i,j,k,2) = ct(i,j,k,2) + mvp(l,2)*wght(l,1)/nvolb
                        ct(i,j,k,3) = ct(i,j,k,3) + mvp(l,3)*wght(l,1)/nvolb
                  endif

                  if (np(ip,j,k) .gt. 0.0) then
                        nvolb = np(ip,j,k)*volb
                        ct(ip,j,k,1) = ct(ip,j,k,1) + mvp(l,1)*wght(l,2)/nvolb
                        ct(ip,j,k,2) = ct(ip,j,k,2) + mvp(l,2)*wght(l,2)/nvolb
                        ct(ip,j,k,3) = ct(ip,j,k,3) + mvp(l,3)*wght(l,2)/nvolb
                  endif

                  if (np(i,j,kp) .gt. 0.0) then
                        nvolb = np(i,j,kp)*volb
                        ct(i,j,kp,1) = ct(i,j,kp,1) + mvp(l,1)*wght(l,3)/nvolb
                        ct(i,j,kp,2) = ct(i,j,kp,2) + mvp(l,2)*wght(l,3)/nvolb
                        ct(i,j,kp,3) = ct(i,j,kp,3) + mvp(l,3)*wght(l,3)/nvolb
                  endif

                  if (np(ip,j,kp) .gt. 0.0) then
                        nvolb = np(ip,j,kp)*volb
                        ct(ip,j,kp,1) = ct(ip,j,kp,1) + mvp(l,1)*wght(l,4)/nvolb
                        ct(ip,j,kp,2) = ct(ip,j,kp,2) + mvp(l,2)*wght(l,4)/nvolb
                        ct(ip,j,kp,3) = ct(ip,j,kp,3) + mvp(l,3)*wght(l,4)/nvolb
                  endif

                  if (np(i,jp,k) .gt. 0.0) then
                        nvolb = np(i,jp,k)*volb
                        ct(i,jp,k,1) = ct(i,jp,k,1) + mvp(l,1)*wght(l,5)/nvolb
                        ct(i,jp,k,2) = ct(i,jp,k,2) + mvp(l,2)*wght(l,5)/nvolb
                        ct(i,jp,k,3) = ct(i,jp,k,3) + mvp(l,3)*wght(l,5)/nvolb
                  endif

                  if (np(ip,jp,k) .gt. 0.0) then
                        nvolb = np(ip,jp,k)*volb
                        ct(ip,jp,k,1) = ct(ip,jp,k,1) + mvp(l,1)*wght(l,6)/nvolb
                        ct(ip,jp,k,2) = ct(ip,jp,k,2) + mvp(l,2)*wght(l,6)/nvolb
                        ct(ip,jp,k,3) = ct(ip,jp,k,3) + mvp(l,3)*wght(l,6)/nvolb
                  endif

                  if (np(i,jp,kp) .gt. 0.0) then
                        nvolb = np(i,jp,kp)*volb
                        ct(i,jp,kp,1) = ct(i,jp,kp,1) + mvp(l,1)*wght(l,7)/nvolb
                        ct(i,jp,kp,2) = ct(i,jp,kp,2) + mvp(l,2)*wght(l,7)/nvolb
                        ct(i,jp,kp,3) = ct(i,jp,kp,3) + mvp(l,3)*wght(l,7)/nvolb
                  endif

                  if (np(ip,jp,kp) .gt. 0.0) then
                        nvolb = np(ip,jp,kp)*volb
                        ct(ip,jp,kp,1) = ct(ip,jp,kp,1) + mvp(l,1)*wght(l,8)/nvolb
                        ct(ip,jp,kp,2) = ct(ip,jp,kp,2) + mvp(l,2)*wght(l,8)/nvolb
                        ct(ip,jp,kp,3) = ct(ip,jp,kp,3) + mvp(l,3)*wght(l,8)/nvolb
                  endif
         
            enddo
            !Used for periodic boundary conditions
            call add_boundary_vector(ct)
            
            call boundary_vector(ct)
            
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(ct(:,:,:,:),recvbuf,count,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
            
            ct(:,:,:,:) = reshape(recvbuf,(/nx,ny,nz,3/))
            
            do i=1,nx-1
                  do j=1,ny-1
                        do k=1,nz-1
                              up_ave(i,j,k,1) = 0.5*(ct(i,j,k,1)+ct(i+1,j,k,1))
                              up_ave(i,j,k,2) = 0.5*(ct(i,j,k,2)+ct(i,j+1,k,2))
                              up_ave(i,j,k,3) = 0.5*(ct(i,j,k,3)+ct(i,j,k+1,3))
                        enddo
                  enddo
            enddo
            
            do i=1,nx
                  do j=1,ny
                        do k=1,nz
                              temp_p(i,j,k) = (1./3.)*1e6*mion*( &
                                    up2(i,j,k,1) - up_ave(i,j,k,1)**2 + &
                                    up2(i,j,k,2) - up_ave(i,j,k,2)**2 + & 
                                    up2(i,j,k,3) - up_ave(i,j,k,3)**2)
                                    
                              tp(i,j,k,1) = (1./3.)*1e6*mion*(up2(i,j,k,1) - up_ave(i,j,k,1)**2)
                              tp(i,j,k,2) = (1./3.)*1e6*mion*(up2(i,j,k,2) - up_ave(i,j,k,2)**2)
                              tp(i,j,k,3) = (1./3.)*1e6*mion*(up2(i,j,k,3) - up_ave(i,j,k,3)**2)
                        enddo
                  enddo
            enddo

            call boundary_scalar(temp_p)
            call boundary_vector(tp)
            
      end subroutine get_temperature
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     subroutine update_np_mixed()
! Weight density to eight nearest grid points
            use dimensions
            use MPI
            use grid, only: dx_grid,dy_grid,dz_grid
            use boundary
            use var_arrays, only: np_mixed,Ni_tot,ijkp,beta,beta_p,wght,mix_ind
            implicit none
            real:: volb, recvbuf(nx*ny*nz)
            integer:: i,j,k,l,ip,jp,kp,ierr,count
            
            count = nx*ny*nz
            
            do i=1,nx
                  do j=1,ny
                        do k=1,nz
                              np_mixed(i,j,k) = 0.0
                        enddo
                  enddo
            enddo

            l=1

            do while (l .le. Ni_tot)
            	do while ( (mix_ind(l) .eq. 0) .and. (l .lt. Ni_tot) )
			      l=l+1
		      end do			
                  if  (l .ge. Ni_tot) then
                        EXIT
                  endif

                  i=ijkp(l,1)
                  j=ijkp(l,2)
                  k=ijkp(l,3)
                  
                  ip=i+1
                  jp=j+1
                  kp=k+1
                  
                  volb = 1.0/(dx_grid(i)*dy_grid(j)*dz_grid(k)*beta*beta_p(l))
                  
                  np_mixed(i,j,k) = np_mixed(i,j,k) + wght(l,1)*volb
                  np_mixed(ip,j,k) = np_mixed(ip,j,k) + wght(l,2)*volb
                  np_mixed(i,j,kp) = np_mixed(i,j,kp) + wght(l,3)*volb
                  np_mixed(ip,j,kp) = np_mixed(ip,j,kp) + wght(l,4)*volb
                  np_mixed(i,jp,k) = np_mixed(i,jp,k) + wght(l,5)*volb
                  np_mixed(ip,jp,k) = np_mixed(ip,jp,k) + wght(l,6)*volb
                  np_mixed(i,jp,kp) = np_mixed(i,jp,kp) + wght(l,7)*volb
                  np_mixed(ip,jp,kp) = np_mixed(ip,jp,kp) + wght(l,8)*volb
                  
                  l=l+1

            end do
                  
            !Used for periodic boundary conditions
            call add_boundary_scalar(np_mixed)      
            
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(np_mixed(:,:,:),recvbuf,count,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
            
            np_mixed(:,:,:) = reshape(recvbuf,(/nx,ny,nz/))
            
            call boundary_scalar(np_mixed)
            
      end subroutine update_np_mixed
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
subroutine get_temperature_mixed()
            use dimensions
            use MPI
            use boundary
            use inputs, only: mion
            use grid, only: qx,qy,qz
            use var_arrays, only: vp,np_mixed,temp_p_mixed,Ni_tot,ijkp,beta,beta_p,wght,mix_ind,tp_mixed
            implicit none
            real:: recvbuf(nx*ny*nz*3),up2(nx,ny,nz,3),up_ave(nx,ny,nz,3),ct(nx,ny,nz,3),volb,nvolb,mvp(Ni_tot,3)
            integer:: i,j,k,l,m,ip,jp,kp,count,ierr
            
            count = nx*ny*nz*3
            
            up2(:,:,:,:) = 0.0
            up_ave(:,:,:,:) = 0.0
            ct(:,:,:,:) = 0.0
            
            do l=1,Ni_tot
            	do m=1,3
                        mvp(l,m) = vp(l,m)
                  enddo
            enddo
            
            l=1
            do while (l .le. Ni_tot)
            	do while ( (mix_ind(l) .eq. 0) .and. (l .lt. Ni_tot) )
			      l=l+1
		      end do
                  if  (l .ge. Ni_tot) then
                        EXIT
                  endif

                  i=ijkp(l,1)
                  j=ijkp(l,2)
                  k=ijkp(l,3)
                  
                  ip=i+1
                  jp=j+1
                  kp=k+1
                  
                  volb = (qx(ip)-qx(i))*(qy(jp)-qy(j))*(qz(kp)-qz(k))*beta*beta_p(l)
                  
                  if (np_mixed(i,j,k) .gt. 0.0) then
                        nvolb = np_mixed(i,j,k)*volb
                        ct(i,j,k,1) = ct(i,j,k,1) + mvp(l,1)**2*wght(l,1)/nvolb  
                        ct(i,j,k,2) = ct(i,j,k,2) + mvp(l,2)**2*wght(l,1)/nvolb
                        ct(i,j,k,3) = ct(i,j,k,3) + mvp(l,3)**2*wght(l,1)/nvolb
                  endif

                  if (np_mixed(ip,j,k) .gt. 0.0) then
                        nvolb = np_mixed(ip,j,k)*volb
                        ct(ip,j,k,1) = ct(ip,j,k,1) + mvp(l,1)**2*wght(l,2)/nvolb
                        ct(ip,j,k,2) = ct(ip,j,k,2) + mvp(l,2)**2*wght(l,2)/nvolb
                        ct(ip,j,k,3) = ct(ip,j,k,3) + mvp(l,3)**2*wght(l,2)/nvolb
                  endif

                  if (np_mixed(i,j,kp) .gt. 0.0) then
                        nvolb = np_mixed(i,j,kp)*volb
                        ct(i,j,kp,1) = ct(i,j,kp,1) + mvp(l,1)**2*wght(l,3)/nvolb
                        ct(i,j,kp,2) = ct(i,j,kp,2) + mvp(l,2)**2*wght(l,3)/nvolb
                        ct(i,j,kp,3) = ct(i,j,kp,3) + mvp(l,3)**2*wght(l,3)/nvolb
                  endif

                  if (np_mixed(ip,j,kp) .gt. 0.0) then
                        nvolb = np_mixed(ip,j,kp)*volb
                        ct(ip,j,kp,1) = ct(ip,j,kp,1) + mvp(l,1)**2*wght(l,4)/nvolb
                        ct(ip,j,kp,2) = ct(ip,j,kp,2) + mvp(l,2)**2*wght(l,4)/nvolb
                        ct(ip,j,kp,3) = ct(ip,j,kp,3) + mvp(l,3)**2*wght(l,4)/nvolb
                  endif

                  if (np_mixed(i,jp,k) .gt. 0.0) then
                        nvolb = np_mixed(i,jp,k)*volb
                        ct(i,jp,k,1) = ct(i,jp,k,1) + mvp(l,1)**2*wght(l,5)/nvolb
                        ct(i,jp,k,2) = ct(i,jp,k,2) + mvp(l,2)**2*wght(l,5)/nvolb
                        ct(i,jp,k,3) = ct(i,jp,k,3) + mvp(l,3)**2*wght(l,5)/nvolb
                  endif

                  if (np_mixed(ip,jp,k) .gt. 0.0) then
                        nvolb = np_mixed(ip,jp,k)*volb
                        ct(ip,jp,k,1) = ct(ip,jp,k,1) + mvp(l,1)**2*wght(l,6)/nvolb
                        ct(ip,jp,k,2) = ct(ip,jp,k,2) + mvp(l,2)**2*wght(l,6)/nvolb
                        ct(ip,jp,k,3) = ct(ip,jp,k,3) + mvp(l,3)**2*wght(l,6)/nvolb
                  endif

                  if (np_mixed(i,jp,kp) .gt. 0.0) then
                        nvolb = np_mixed(i,jp,kp)*volb
                        ct(i,jp,kp,1) = ct(i,jp,kp,1) + mvp(l,1)**2*wght(l,7)/nvolb
                        ct(i,jp,kp,2) = ct(i,jp,kp,2) + mvp(l,2)**2*wght(l,7)/nvolb
                        ct(i,jp,kp,3) = ct(i,jp,kp,3) + mvp(l,3)**2*wght(l,7)/nvolb
                  endif

                  if (np_mixed(ip,jp,kp) .gt. 0.0) then
                        nvolb = np_mixed(ip,jp,kp)*volb
                        ct(ip,jp,kp,1) = ct(ip,jp,kp,1) + mvp(l,1)**2*wght(l,8)/nvolb
                        ct(ip,jp,kp,2) = ct(ip,jp,kp,2) + mvp(l,2)**2*wght(l,8)/nvolb
                        ct(ip,jp,kp,3) = ct(ip,jp,kp,3) + mvp(l,3)**2*wght(l,8)/nvolb
                  endif
                  
                  l=l+1
         
         end do
            
            !Used for periodic boundary conditions
            call add_boundary_vector(ct)
            
            call boundary_vector(ct)
                        
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(ct(:,:,:,:),recvbuf,count,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
            
            ct(:,:,:,:) = reshape(recvbuf,(/nx,ny,nz,3/))
            
            do i=1,nx-1
                  do j=1,ny-1
                        do k=1,nz-1
                              up2(i,j,k,1) = 0.5*(ct(i,j,k,1)+ct(i+1,j,k,1))
                              up2(i,j,k,2) = 0.5*(ct(i,j,k,2)+ct(i,j+1,k,2))
                              up2(i,j,k,3) = 0.5*(ct(i,j,k,3)+ct(i,j,k+1,3))
                        enddo
                  enddo
            enddo
            
            call boundary_vector(up2)
            
            ct(:,:,:,:) = 0.0

            l=1
            do while (l .le. Ni_tot) 
            	do while ( (mix_ind(l) .eq. 0) .and. (l .lt. Ni_tot) )
			      l=l+1
		      end do
		      if  (l .ge. Ni_tot) then
			      EXIT
		      endif
                  i=ijkp(l,1)
                  j=ijkp(l,2)
                  k=ijkp(l,3)
                  
                  ip=i+1
                  jp=j+1
                  kp=k+1
                  
                  volb = (qx(ip)-qx(i))*(qy(jp)-qy(j))*(qz(kp)-qz(k))*beta*beta_p(l)
                  
                  if (np_mixed(i,j,k) .gt. 0.0) then
                        nvolb = np_mixed(i,j,k)*volb
                        ct(i,j,k,1) = ct(i,j,k,1) + mvp(l,1)*wght(l,1)/nvolb  
                        ct(i,j,k,2) = ct(i,j,k,2) + mvp(l,2)*wght(l,1)/nvolb
                        ct(i,j,k,3) = ct(i,j,k,3) + mvp(l,3)*wght(l,1)/nvolb
                  endif

                  if (np_mixed(ip,j,k) .gt. 0.0) then
                        nvolb = np_mixed(ip,j,k)*volb
                        ct(ip,j,k,1) = ct(ip,j,k,1) + mvp(l,1)*wght(l,2)/nvolb
                        ct(ip,j,k,2) = ct(ip,j,k,2) + mvp(l,2)*wght(l,2)/nvolb
                        ct(ip,j,k,3) = ct(ip,j,k,3) + mvp(l,3)*wght(l,2)/nvolb
                  endif

                  if (np_mixed(i,j,kp) .gt. 0.0) then
                        nvolb = np_mixed(i,j,kp)*volb
                        ct(i,j,kp,1) = ct(i,j,kp,1) + mvp(l,1)*wght(l,3)/nvolb
                        ct(i,j,kp,2) = ct(i,j,kp,2) + mvp(l,2)*wght(l,3)/nvolb
                        ct(i,j,kp,3) = ct(i,j,kp,3) + mvp(l,3)*wght(l,3)/nvolb
                  endif

                  if (np_mixed(ip,j,kp) .gt. 0.0) then
                        nvolb = np_mixed(ip,j,kp)*volb
                        ct(ip,j,kp,1) = ct(ip,j,kp,1) + mvp(l,1)*wght(l,4)/nvolb
                        ct(ip,j,kp,2) = ct(ip,j,kp,2) + mvp(l,2)*wght(l,4)/nvolb
                        ct(ip,j,kp,3) = ct(ip,j,kp,3) + mvp(l,3)*wght(l,4)/nvolb
                  endif

                  if (np_mixed(i,jp,k) .gt. 0.0) then
                        nvolb = np_mixed(i,jp,k)*volb
                        ct(i,jp,k,1) = ct(i,jp,k,1) + mvp(l,1)*wght(l,5)/nvolb
                        ct(i,jp,k,2) = ct(i,jp,k,2) + mvp(l,2)*wght(l,5)/nvolb
                        ct(i,jp,k,3) = ct(i,jp,k,3) + mvp(l,3)*wght(l,5)/nvolb
                  endif

                  if (np_mixed(ip,jp,k) .gt. 0.0) then
                        nvolb = np_mixed(ip,jp,k)*volb
                        ct(ip,jp,k,1) = ct(ip,jp,k,1) + mvp(l,1)*wght(l,6)/nvolb
                        ct(ip,jp,k,2) = ct(ip,jp,k,2) + mvp(l,2)*wght(l,6)/nvolb
                        ct(ip,jp,k,3) = ct(ip,jp,k,3) + mvp(l,3)*wght(l,6)/nvolb
                  endif

                  if (np_mixed(i,jp,kp) .gt. 0.0) then
                        nvolb = np_mixed(i,jp,kp)*volb
                        ct(i,jp,kp,1) = ct(i,jp,kp,1) + mvp(l,1)*wght(l,7)/nvolb
                        ct(i,jp,kp,2) = ct(i,jp,kp,2) + mvp(l,2)*wght(l,7)/nvolb
                        ct(i,jp,kp,3) = ct(i,jp,kp,3) + mvp(l,3)*wght(l,7)/nvolb
                  endif

                  if (np_mixed(ip,jp,kp) .gt. 0.0) then
                        nvolb = np_mixed(ip,jp,kp)*volb
                        ct(ip,jp,kp,1) = ct(ip,jp,kp,1) + mvp(l,1)*wght(l,8)/nvolb
                        ct(ip,jp,kp,2) = ct(ip,jp,kp,2) + mvp(l,2)*wght(l,8)/nvolb
                        ct(ip,jp,kp,3) = ct(ip,jp,kp,3) + mvp(l,3)*wght(l,8)/nvolb
                  endif

                  l=l+1
            enddo

            !Used for periodic boundary conditions
            call add_boundary_vector(ct)
            
            call boundary_vector(ct)
            
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(ct(:,:,:,:),recvbuf,count,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
            
            ct(:,:,:,:) = reshape(recvbuf,(/nx,ny,nz,3/))
            
            do i=1,nx-1
                  do j=1,ny-1
                        do k=1,nz-1
                              up_ave(i,j,k,1) = 0.5*(ct(i,j,k,1)+ct(i+1,j,k,1))
                              up_ave(i,j,k,2) = 0.5*(ct(i,j,k,2)+ct(i,j+1,k,2))
                              up_ave(i,j,k,3) = 0.5*(ct(i,j,k,3)+ct(i,j,k+1,3))
                        enddo
                  enddo
            enddo
            
            do i=1,nx
                  do j=1,ny
                        do k=1,nz
                              temp_p_mixed(i,j,k) = (1./3.)*1e6*mion*( &
                                    up2(i,j,k,1) - up_ave(i,j,k,1)**2 + &
                                    up2(i,j,k,2) - up_ave(i,j,k,2)**2 + & 
                                    up2(i,j,k,3) - up_ave(i,j,k,3)**2)
                              tp_mixed(i,j,k,1) = (1./3.)*1e6*mion*(up2(i,j,k,1) - up_ave(i,j,k,1)**2)
                              tp_mixed(i,j,k,2) = (1./3.)*1e6*mion*(up2(i,j,k,2) - up_ave(i,j,k,2)**2)
                              tp_mixed(i,j,k,3) = (1./3.)*1e6*mion*(up2(i,j,k,3) - up_ave(i,j,k,3)**2)
                        enddo
                  enddo
            enddo

            call boundary_scalar(temp_p_mixed)
            call boundary_vector(tp_mixed)
            
      end subroutine get_temperature_mixed
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     subroutine update_np_cold()
! Weight density to eight nearest grid points
            use dimensions
            use MPI
            use grid, only: dx_grid,dy_grid,dz_grid
            use boundary
            use var_arrays, only: np_cold,Ni_tot,ijkp,beta,beta_p,wght,mix_ind
            use mult_proc
            implicit none
            real:: volb, recvbuf(nx*ny*nz)
            integer:: i,j,k,l,ip,jp,kp,ierr,count
            
            count = nx*ny*nz

            do i=1,nx
                  do j=1,ny
                        do k=1,nz
                              np_cold(i,j,k) = 0.0
                        enddo
                  enddo
            enddo

	      l=1
           
            do while (l .le. Ni_tot)
            	do while ( (mix_ind(l) .ne. 0) .and. (l .lt. Ni_tot) )
			      l=l+1
		      end do
		      if  (l .ge. Ni_tot) then
			      EXIT
		      endif
		
                  i=ijkp(l,1)
                  j=ijkp(l,2)
                  k=ijkp(l,3)
                  
                  ip=i+1
                  jp=j+1
                  kp=k+1
                  
                  volb = 1.0/(dx_grid(i)*dy_grid(j)*dz_grid(k)*beta*beta_p(l))
                  
                  np_cold(i,j,k) = np_cold(i,j,k) + wght(l,1)*volb
                  np_cold(ip,j,k) = np_cold(ip,j,k) + wght(l,2)*volb
                  np_cold(i,j,kp) = np_cold(i,j,kp) + wght(l,3)*volb
                  np_cold(ip,j,kp) = np_cold(ip,j,kp) + wght(l,4)*volb
                  np_cold(i,jp,k) = np_cold(i,jp,k) + wght(l,5)*volb
                  np_cold(ip,jp,k) = np_cold(ip,jp,k) + wght(l,6)*volb
                  np_cold(i,jp,kp) = np_cold(i,jp,kp) + wght(l,7)*volb
                  np_cold(ip,jp,kp) = np_cold(ip,jp,kp) + wght(l,8)*volb
                  
                  l=l+1  

            enddo
            
            !Used for periodic boundary conditions
            call add_boundary_scalar(np_cold)
            
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(np_cold(:,:,:),recvbuf,count,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
            
            np_cold(:,:,:) = reshape(recvbuf,(/nx,ny,nz/))
            
            call boundary_scalar(np_cold)
            
      end subroutine update_np_cold
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
subroutine get_temperature_cold()
            use dimensions
            use MPI
            use boundary
            use inputs, only: mion
            use grid, only: qx,qy,qz
            use var_arrays, only: vp,np_cold,temp_p_cold,Ni_tot,ijkp,beta,beta_p,wght,mix_ind,tp_cold
            implicit none
            real:: recvbuf(nx*ny*nz*3),up2(nx,ny,nz,3),up_ave(nx,ny,nz,3),ct(nx,ny,nz,3),volb,nvolb,mvp(Ni_tot,3)
            integer:: i,j,k,l,m,ip,jp,kp,count,ierr
            
            count = nx*ny*nz*3
            
            up2(:,:,:,:) = 0.0
            up_ave(:,:,:,:) = 0.0
            ct(:,:,:,:) = 0.0
            
            
            do l=1,Ni_tot
            	do m=1,3
                        mvp(l,m) = vp(l,m)
                  enddo
            enddo

            l=1
            do while (l .le. Ni_tot)
            	do while ( (mix_ind(l) .ne. 0) .and. (l .lt. Ni_tot) )
			      l=l+1
		      end do
		      if  (l .ge. Ni_tot) then
		      	EXIT
		      endif

                  i=ijkp(l,1)
                  j=ijkp(l,2)
                  k=ijkp(l,3)
                  
                  ip=i+1
                  jp=j+1
                  kp=k+1
                  
                  volb = (qx(ip)-qx(i))*(qy(jp)-qy(j))*(qz(kp)-qz(k))*beta*beta_p(l)
                  
                  if (np_cold(i,j,k) .gt. 0.0) then
                        nvolb = np_cold(i,j,k)*volb
                        ct(i,j,k,1) = ct(i,j,k,1) + mvp(l,1)**2*wght(l,1)/nvolb  
                        ct(i,j,k,2) = ct(i,j,k,2) + mvp(l,2)**2*wght(l,1)/nvolb
                        ct(i,j,k,3) = ct(i,j,k,3) + mvp(l,3)**2*wght(l,1)/nvolb
                  endif

                  if (np_cold(ip,j,k) .gt. 0.0) then
                        nvolb = np_cold(ip,j,k)*volb
                        ct(ip,j,k,1) = ct(ip,j,k,1) + mvp(l,1)**2*wght(l,2)/nvolb
                        ct(ip,j,k,2) = ct(ip,j,k,2) + mvp(l,2)**2*wght(l,2)/nvolb
                        ct(ip,j,k,3) = ct(ip,j,k,3) + mvp(l,3)**2*wght(l,2)/nvolb
                  endif

                  if (np_cold(i,j,kp) .gt. 0.0) then
                        nvolb = np_cold(i,j,kp)*volb
                        ct(i,j,kp,1) = ct(i,j,kp,1) + mvp(l,1)**2*wght(l,3)/nvolb
                        ct(i,j,kp,2) = ct(i,j,kp,2) + mvp(l,2)**2*wght(l,3)/nvolb
                        ct(i,j,kp,3) = ct(i,j,kp,3) + mvp(l,3)**2*wght(l,3)/nvolb
                  endif

                  if (np_cold(ip,j,kp) .gt. 0.0) then
                        nvolb = np_cold(ip,j,kp)*volb
                        ct(ip,j,kp,1) = ct(ip,j,kp,1) + mvp(l,1)**2*wght(l,4)/nvolb
                        ct(ip,j,kp,2) = ct(ip,j,kp,2) + mvp(l,2)**2*wght(l,4)/nvolb
                        ct(ip,j,kp,3) = ct(ip,j,kp,3) + mvp(l,3)**2*wght(l,4)/nvolb
                  endif

                  if (np_cold(i,jp,k) .gt. 0.0) then
                        nvolb = np_cold(i,jp,k)*volb
                        ct(i,jp,k,1) = ct(i,jp,k,1) + mvp(l,1)**2*wght(l,5)/nvolb
                        ct(i,jp,k,2) = ct(i,jp,k,2) + mvp(l,2)**2*wght(l,5)/nvolb
                        ct(i,jp,k,3) = ct(i,jp,k,3) + mvp(l,3)**2*wght(l,5)/nvolb
                  endif

                  if (np_cold(ip,jp,k) .gt. 0.0) then
                        nvolb = np_cold(ip,jp,k)*volb
                        ct(ip,jp,k,1) = ct(ip,jp,k,1) + mvp(l,1)**2*wght(l,6)/nvolb
                        ct(ip,jp,k,2) = ct(ip,jp,k,2) + mvp(l,2)**2*wght(l,6)/nvolb
                        ct(ip,jp,k,3) = ct(ip,jp,k,3) + mvp(l,3)**2*wght(l,6)/nvolb
                  endif

                  if (np_cold(i,jp,kp) .gt. 0.0) then
                        nvolb = np_cold(i,jp,kp)*volb
                        ct(i,jp,kp,1) = ct(i,jp,kp,1) + mvp(l,1)**2*wght(l,7)/nvolb
                        ct(i,jp,kp,2) = ct(i,jp,kp,2) + mvp(l,2)**2*wght(l,7)/nvolb
                        ct(i,jp,kp,3) = ct(i,jp,kp,3) + mvp(l,3)**2*wght(l,7)/nvolb
                  endif

                  if (np_cold(ip,jp,kp) .gt. 0.0) then
                        nvolb = np_cold(ip,jp,kp)*volb
                        ct(ip,jp,kp,1) = ct(ip,jp,kp,1) + mvp(l,1)**2*wght(l,8)/nvolb
                        ct(ip,jp,kp,2) = ct(ip,jp,kp,2) + mvp(l,2)**2*wght(l,8)/nvolb
                        ct(ip,jp,kp,3) = ct(ip,jp,kp,3) + mvp(l,3)**2*wght(l,8)/nvolb
                  endif

                  l=l+1

         enddo
            
            !Used for periodic boundary conditions
            call add_boundary_vector(ct)

            call boundary_vector(ct)
            
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(ct(:,:,:,:),recvbuf,count,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
            
            ct(:,:,:,:) = reshape(recvbuf,(/nx,ny,nz,3/))
            
            do i=1,nx-1
                  do j=1,ny-1
                        do k=1,nz-1
                              up2(i,j,k,1) = 0.5*(ct(i,j,k,1)+ct(i+1,j,k,1))
                              up2(i,j,k,2) = 0.5*(ct(i,j,k,2)+ct(i,j+1,k,2))
                              up2(i,j,k,3) = 0.5*(ct(i,j,k,3)+ct(i,j,k+1,3))
                        enddo
                  enddo
            enddo
            
            call boundary_vector(up2)
            
            ct(:,:,:,:) = 0.0

            l=1
            do while (l .le. Ni_tot) 
            	do while ( (mix_ind(l) .ne. 0) .and. (l .lt. Ni_tot) )
			      l=l+1
		      end do
		      if  (l .ge. Ni_tot) then
			      EXIT
		      endif
		
                  i=ijkp(l,1)
                  j=ijkp(l,2)
                  k=ijkp(l,3)
                  
                  ip=i+1
                  jp=j+1
                  kp=k+1
                  
                  volb = (qx(ip)-qx(i))*(qy(jp)-qy(j))*(qz(kp)-qz(k))*beta*beta_p(l)
                  
                  if (np_cold(i,j,k) .gt. 0.0) then
                        nvolb = np_cold(i,j,k)*volb
                        ct(i,j,k,1) = ct(i,j,k,1) + mvp(l,1)*wght(l,1)/nvolb  
                        ct(i,j,k,2) = ct(i,j,k,2) + mvp(l,2)*wght(l,1)/nvolb
                        ct(i,j,k,3) = ct(i,j,k,3) + mvp(l,3)*wght(l,1)/nvolb
                  endif

                  if (np_cold(ip,j,k) .gt. 0.0) then
                        nvolb = np_cold(ip,j,k)*volb
                        ct(ip,j,k,1) = ct(ip,j,k,1) + mvp(l,1)*wght(l,2)/nvolb
                        ct(ip,j,k,2) = ct(ip,j,k,2) + mvp(l,2)*wght(l,2)/nvolb
                        ct(ip,j,k,3) = ct(ip,j,k,3) + mvp(l,3)*wght(l,2)/nvolb
                  endif

                  if (np_cold(i,j,kp) .gt. 0.0) then
                        nvolb = np_cold(i,j,kp)*volb
                        ct(i,j,kp,1) = ct(i,j,kp,1) + mvp(l,1)*wght(l,3)/nvolb
                        ct(i,j,kp,2) = ct(i,j,kp,2) + mvp(l,2)*wght(l,3)/nvolb
                        ct(i,j,kp,3) = ct(i,j,kp,3) + mvp(l,3)*wght(l,3)/nvolb
                  endif

                  if (np_cold(ip,j,kp) .gt. 0.0) then
                        nvolb = np_cold(ip,j,kp)*volb
                        ct(ip,j,kp,1) = ct(ip,j,kp,1) + mvp(l,1)*wght(l,4)/nvolb
                        ct(ip,j,kp,2) = ct(ip,j,kp,2) + mvp(l,2)*wght(l,4)/nvolb
                        ct(ip,j,kp,3) = ct(ip,j,kp,3) + mvp(l,3)*wght(l,4)/nvolb
                  endif

                  if (np_cold(i,jp,k) .gt. 0.0) then
                        nvolb = np_cold(i,jp,k)*volb
                        ct(i,jp,k,1) = ct(i,jp,k,1) + mvp(l,1)*wght(l,5)/nvolb
                        ct(i,jp,k,2) = ct(i,jp,k,2) + mvp(l,2)*wght(l,5)/nvolb
                        ct(i,jp,k,3) = ct(i,jp,k,3) + mvp(l,3)*wght(l,5)/nvolb
                  endif

                  if (np_cold(ip,jp,k) .gt. 0.0) then
                        nvolb = np_cold(ip,jp,k)*volb
                        ct(ip,jp,k,1) = ct(ip,jp,k,1) + mvp(l,1)*wght(l,6)/nvolb
                        ct(ip,jp,k,2) = ct(ip,jp,k,2) + mvp(l,2)*wght(l,6)/nvolb
                        ct(ip,jp,k,3) = ct(ip,jp,k,3) + mvp(l,3)*wght(l,6)/nvolb
                  endif

                  if (np_cold(i,jp,kp) .gt. 0.0) then
                        nvolb = np_cold(i,jp,kp)*volb
                        ct(i,jp,kp,1) = ct(i,jp,kp,1) + mvp(l,1)*wght(l,7)/nvolb
                        ct(i,jp,kp,2) = ct(i,jp,kp,2) + mvp(l,2)*wght(l,7)/nvolb
                        ct(i,jp,kp,3) = ct(i,jp,kp,3) + mvp(l,3)*wght(l,7)/nvolb
                  endif

                  if (np_cold(ip,jp,kp) .gt. 0.0) then
                        nvolb = np_cold(ip,jp,kp)*volb
                        ct(ip,jp,kp,1) = ct(ip,jp,kp,1) + mvp(l,1)*wght(l,8)/nvolb
                        ct(ip,jp,kp,2) = ct(ip,jp,kp,2) + mvp(l,2)*wght(l,8)/nvolb
                        ct(ip,jp,kp,3) = ct(ip,jp,kp,3) + mvp(l,3)*wght(l,8)/nvolb
                  endif
                  
                  l=l+1
         
         enddo
            
            !Used for periodic boundary conditions
            call add_boundary_vector(ct)
            
            call boundary_vector(ct)
            
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(ct(:,:,:,:),recvbuf,count,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
            
            ct(:,:,:,:) = reshape(recvbuf,(/nx,ny,nz,3/))
            
            do i=1,nx-1
                  do j=1,ny-1
                        do k=1,nz-1
                              up_ave(i,j,k,1) = 0.5*(ct(i,j,k,1)+ct(i+1,j,k,1))
                              up_ave(i,j,k,2) = 0.5*(ct(i,j,k,2)+ct(i,j+1,k,2))
                              up_ave(i,j,k,3) = 0.5*(ct(i,j,k,3)+ct(i,j,k+1,3))
                        enddo
                  enddo
            enddo
            
            do i=1,nx
                  do j=1,ny
                        do k=1,nz
                              temp_p_cold(i,j,k) = (1./3.)*1e6*mion*( &
                                    up2(i,j,k,1) - up_ave(i,j,k,1)**2 + &
                                    up2(i,j,k,2) - up_ave(i,j,k,2)**2 + & 
                                    up2(i,j,k,3) - up_ave(i,j,k,3)**2)
                                    
                              tp_cold(i,j,k,1) = (1./3.)*1e6*mion*(up2(i,j,k,1) - up_ave(i,j,k,1)**2)
                              tp_cold(i,j,k,2) = (1./3.)*1e6*mion*(up2(i,j,k,2) - up_ave(i,j,k,2)**2)
                              tp_cold(i,j,k,3) = (1./3.)*1e6*mion*(up2(i,j,k,3) - up_ave(i,j,k,3)**2)
                        enddo
                  enddo
            enddo

	    
            call boundary_scalar(temp_p_cold)
            call boundary_vector(tp_cold)
            
      end subroutine get_temperature_cold
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_pindex(l)
            use dimensions
            use inputs, only: dx,dy
            use grid, only: qx,qy,qz
            use var_arrays, only: ijkp, xp
            implicit none
            integer, intent(in):: l
            integer:: hi,mid,i,j,k

            i=1
            hi = nx
            do
                  mid = (i+hi)/2
                  if (xp(l,1) .lt. qx(mid)) then
                        hi=mid
                  else
                        i=mid
                  endif
                  if (i+1 .ge. hi) exit
            enddo

            ijkp(l,1)=i

            j = floor(xp(l,2)/dy)
            ijkp(l,2) = j

            k=1
            hi = nz
            do
                  mid = (k+hi)/2
                  if (xp(l,3) .lt. qz(mid)) then
                        hi=mid
                  else
                        k=mid
                  endif
                  if (k+1 .ge. hi) exit
            enddo
            ijkp(l,3) = k


            
      end subroutine get_pindex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine particle_boundary()
            use dimensions
            use boundary
            use inputs, only: boundx, PI, vsw, dx, dy, km_to_m, &
            beta_particle, kboltz, mion, nf_init, &
            b0_init,mu0,boundx, q, mO, va_f, vth, va!, removed
            use var_arrays, only: xp, vp, Ni_tot, b0, mix_ind
            use grid, only: qx,qy,qz
            implicit none
            integer:: l
            REAL:: vx, vy, vz, va_x
            
            if (boundx .eq. 1) then   !Fully periodic
                  where (xp(:,1) .gt. qx(nx-1))
                        xp(:,1) = qx(1) + (xp(:,1) - qx(nx-1))
                  endwhere
                  where (xp(:,1) .le. qx(1))
                        xp(:,1) = qx(nx-1) - (qx(1) - xp(:,1))
                  endwhere
                  where (xp(:,2) .gt. qy(ny-1))
                        xp(:,2) = qy(1) + (xp(:,2) - qy(ny-1))
                  endwhere
                  where (xp(:,2) .le. qy(1))
                        xp(:,2) = qy(ny-1) - (qy(1) - xp(:,2))
                  endwhere
                  where (xp(:,3) .gt. qz(nz-1))
                        xp(:,3) = qz(1) + (xp(:,3) - qz(nz-1))
                  endwhere
                  where (xp(:,3) .le. qz(1))
                        xp(:,3) = qx(nx-1) - (qx(1) - xp(:,1))
                  endwhere
                  
            endif
            
            if (boundx .eq. 2) then      ! Periodic in x and y
                  where (xp(:,1) .gt. qx(nx-1))
                        xp(:,1) = qx(1) + (xp(:,1) - qx(nx-1))
                  endwhere

                  where (xp(:,1) .le. qx(1))
                        xp(:,1) = qx(nx-1) - (qx(1) - xp(:,1))
                  endwhere


                  where (xp(:,2) .gt. qy(ny-1))
                        xp(:,2) = qy(1) + (xp(:,2) - qy(ny-1))
                  endwhere

                  where (xp(:,2) .le. qy(1))
                        xp(:,2) = qy(ny-1) - (qy(1) - xp(:,2))
                  endwhere
                  
                  ! Particles are put back in at opposite velocity in a random location in x and y
                  do l=1,Ni_tot
                        if (xp(l,3) .le. qz(1)) then
                              vp(l,3) = -vp(l,3)
                              xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
                              xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
                              xp(l,3) = qz(1)+(qz(1) - xp(l,3))
                        else if (xp(l,3) .ge. qz(nz-1)) then
                              vp(l,3) = -vp(l,3)
                              xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
                              xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
                              xp(l,3) = qz(nz-1)-(xp(l,3) - qz(nz-1))
                  endif
                  enddo

            endif


            ! Local Simulations, with injection of FS on right side
            if (boundx .eq. 4) then   !Periodic in y only, remove particles that are outside of x and z plane

            !Y boundary
                  where (xp(:,2) .gt. qy(ny-1))
                        xp(:,2) = qy(1) + (xp(:,2) - qy(ny-1))
                  endwhere
                  where (xp(:,2) .le. qy(1))
                        xp(:,2) = qy(ny-1) - (qy(1) - xp(:,2))
                  endwhere
                  
      
                  !The particle that leaves z boundaries is reflected but placed at a different spot along the same boundary,
                  !particles that leave the x boundaries are removed
                  do l=Ni_tot,1,-1
                  
                        !x boundaries
                        if ( (xp(l,1) .lt. 0) )then!  .and. (vp(l,1) .lt. 0) )  then
                              call remove_ion(l)
                        else if ( ( xp(l,1) .gt. qx(nx) )) then! .and. ( vp(l,1) .gt. 0 ) .and. (mix_ind(l) .eq. 1) ) then
                              call remove_ion(l)     
                        
                        !z boundaries      
                        else if (xp(l,3) .ge. (qz(nz-1)) ) then
                              vp(l,3) = -vp(l,3)

                              xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
                              xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
                              xp(l,3) = qz(nz-1)-(xp(l,3) - qz(nz-1))

                        else if (xp(l,3) .le. qz(1)) then
                              vp(l,3) = -vp(l,3)

                              xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
                              xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
                              xp(l,3) = qz(1)+(qz(1) - xp(l,3))

                        endif
            
                  enddo

            endif

            ! Planar Shock simuilations, with shock generation on right side.
            if (boundx .eq. 5) then  !Periodic in y and z, for shock generation boundary - particles that leave left boundary are removed, right boundary are reflected
                  
                  !Y Direction
                  where (xp(:,2) .gt. qy(ny-1))
                        xp(:,2) = qy(1) + (xp(:,2) - qy(ny-1))
                  endwhere
                  where (xp(:,2) .le. qy(1))
                        xp(:,2) = qy(ny-1) - (qy(1) - xp(:,2))
                  endwhere
                  
                  !Z Direction
                  where (xp(:,3) .gt. qz(nz-1))
                        xp(:,3) = qz(1) + (xp(:,3) - qz(nz-1))
                  endwhere
                  where (xp(:,3) .le. qz(1))
                        xp(:,3) = qz(nz-1) - (qz(1) - xp(:,3))
                  endwhere

                  do l=Ni_tot,1,-1
                        !x boundaries
                        if ( (xp(l,1) .lt. (qx(2)-qx(1)) )  .and. ( vp(l,1) .lt. 0  ) )then!  .and. (vp(l,1) .lt. 0) )  then
                              call remove_ion(l)                 
                        else if ( ( xp(l,1) .ge. qx(nx-1) )) then! .and. ( vp(l,1) .gt. 0 ) .and. (mix_ind(l) .eq. 1) ) then      
                              vp(l,1) = -vp(l,1)
                              xp(l,1) = qx(nx-1)-(xp(l,1) - qx(nx-1))
                              mix_ind(l) = 1 !reflected solar wind ions are not marked as "mixed" or as foreshock ions
                        endif
            
                  enddo

            endif             

      end subroutine particle_boundary
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module gutsp
