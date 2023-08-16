module gutsf
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine f_update_tlev(b1,b12,b1p2,bt,b0) !loops run 1 o n since values are only being copied
            use dimensions
            use boundary
            use inputs, only: q, mO,pi,ddthickness
            use grid
            implicit none
            real, intent(in):: b0(nx,ny,nz,3)
            real, intent(inout):: b1(nx,ny,nz,3),b1p2(nx,ny,nz,3)
            real, intent(out):: bt(nx,ny,nz,3), b12(nx,ny,nz,3)
            
            integer:: i,j,k,m
            real:: eoverm,dtheta!,mPerCell
            
            dtheta = 2.0*pi / 4.0 /(2.0*ddthickness)
            
            eoverm = q/mO
            
            do i=1,nx
                  do j=1,ny
                        do k= 1,nz
                              bt(i,j,k,1) = b1p2(i,j,k,1) + b0(i,j,k,1)
                              bt(i,j,k,2) = b1p2(i,j,k,2) + b0(i,j,k,2)
                              bt(i,j,k,3) = b1p2(i,j,k,3) + b0(i,j,k,3) 
                              do m=1,3
                                    b12(i,j,k,m)= b1(i,j,k,m)
                                    b1(i,j,k,m) = b1p2(i,j,k,m)
                              enddo
                        enddo
                  enddo
            enddo
            
            call obstacle_boundary_B(bt)
            call obstacle_boundary_B(b12)
            call obstacle_boundary_B(b1)

      end subroutine f_update_tlev
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine crossf2(aa,btc,cc)
!The cross product is formed at the main cell center.  aa and btc must be given already 
!extrapolated to the main cell center.
            use dimensions
            use boundary
            implicit none
            real, intent(inout):: aa(nx,ny,nz,3), btc(nx,ny,nz,3)
            real, intent(out):: cc(nx,ny,nz,3)
            real:: ax,ay,az,bx,by,bz, ct(nx,ny,nz,3)
            integer:: i,j,k,ip,jp,kp

            !real, allocatable :: ct(:,:,:,:)
            !allocate(ct(nx,ny,nz,3))


            call boundary_vector(aa)
            call boundary_vector(btc)

            do i=2,nx-1
                  do j=2,ny-1
                        do k=2,nz-1
                              ax = aa(i,j,k,1)
                              bx = btc(i,j,k,1)
                              
                              ay = aa(i,j,k,2)
                              by = btc(i,j,k,2)
                              
                              az = aa(i,j,k,3)
                              bz = btc(i,j,k,3)
                              
                              ct(i,j,k,1) = ay*bz - az*by
                              ct(i,j,k,2) = az*bx - ax*bz
                              ct(i,j,k,3) = ax*by - ay*bx
                              
                        enddo
                  enddo
            enddo
            
            call boundary_vector(ct)
            
! Extrapolate back to the main cell contravariant positions.
! Just average across cells since cell edges are centered about the grid points

            do i=2,nx-1
                  do j=2,ny-1
                        do k = 2, nz-1
                              ip = i+1
                              jp = j+1
                              kp = k+1
                              
                              cc(i,j,k,1) = 0.5*(ct(i,j,k,1)+ct(ip,j,k,1))
                              cc(i,j,k,2) = 0.5*(ct(i,j,k,2)+ct(i,jp,k,2))
                              cc(i,j,k,3) = 0.5*(ct(i,j,k,3)+ct(i,j,kp,3))
                        enddo
                  enddo
            enddo
            
            call boundary_vector(cc)
            
      end subroutine crossf2
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine curlB(b1,np,aj)
! Calculates curl B / n*alpha.  The resulting "current" is called aj
! which is used in several other places in the code.  This curl is 
! performed on the main cell where B is covarient.  The resulting
! current is main cell contravarient.  Note that dz_cell is used for
! the cell dimensions since dz_grid is not equal to dz_cell on non-
! uniform grid.
            use dimensions
            use boundary
            use grid, only: dx_cell,dy_cell,dz_cell
            use inputs, only: alpha
            use var_arrays, only: curlBcurrent
            implicit none
            real, intent(in):: np(nx,ny,nz)
            real, intent(inout):: b1(nx,ny,nz,3)
            real, intent(out) :: aj(nx,ny,nz,3)
            real:: curl_B(3), ntot(3)
            integer:: i,j,k,m,ip,jp,kp
            
            call boundary_vector(b1)

            do i=2,nx-1
                  do j=2,ny-1
                        do k=2,nz-1
                                    
                              ip = i+1
                              jp = j+1
                              kp = k+1
                              
                              ntot(1) = 0.5*(np(i,j,k)+np(ip,j,k))
                              ntot(2) = 0.5*(np(i,j,k)+np(i,jp,k))
                              ntot(3) = 0.5*(np(i,j,k)+np(i,j,kp))
                              
                              
                              curl_B(1) = (b1(i,j,k,3) - &
                                    b1(i,j-1,k,3))/dy_cell(j) + &
                                    (b1(i,j,k-1,2) - b1(i,j,k,2))/dz_cell(k) 
                              curl_B(2) = (b1(i,j,k,1) - &
                                    b1(i,j,k-1,1))/dz_cell(k) + &
                                    (b1(i-1,j,k,3) - b1(i,j,k,3))/dx_cell(i)
                              curl_B(3) = (b1(i,j,k,2) - &
                                    b1(i-1,j,k,2))/dx_cell(i) + &
                                    (b1(i,j-1,k,1) - b1(i,j,k,1))/dy_cell(j)
                                   
                                    
                              do m = 1,3
                                    aj(i,j,k,m) = curl_B(m)/(ntot(m)*alpha)
                                    curlBcurrent(i,j,k,m) = curl_B(m)
                              enddo
                        enddo
                  enddo
            enddo

      end subroutine curlB
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine curlE(E,curl_E)
! E is dual cell covarient, and curl_E will be returned as main
! cell covarient...as all magnetic fields are.  All i,j,k exclude
! boundaries.  Boundaries are taken care of in main fluid code.
            use dimensions
            use grid, only: dx_grid,dy_grid,dz_grid
            implicit none
            real, intent(in):: E(nx,ny,nz,3)
            real, intent(out):: curl_E(nx,ny,nz,3)
            integer:: i,j,k
      
            do i=2,nx-1
                  do j=2,ny-1
                        do k=2,nz-1
                              curl_E(i,j,k,1) = (E(i,j+1,k,3) - E(i,j,k,3))/dy_grid(j) &
                                    + (E(i,j,k,2) - E(i,j,k+1,2))/dz_grid(k)
                              curl_E(i,j,k,2) = (E(i,j,k,3) - E(i+1,j,k,3))/dx_grid(i) &
                                    + (E(i,j,k+1,1) - E(i,j,k,1))/dz_grid(k)
                              curl_E(i,j,k,3) = (E(i,j,k,1) - E(i,j+1,k,1))/dy_grid(j) &
                                    + (E(i+1,j,k,2) - E(i,j,k,2))/dx_grid(i)
                              
                        enddo
                  enddo
            enddo
                        
      end subroutine curlE
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_E(E,bt,aj,up,nu)
! E must be at time level m. We have uf at levels m-1/2 and m+1/2, so
! the average value is used for uf in the calculation of ui.
            use dimensions
            use boundary
            use grid_interp
            use var_arrays, only: gradP, ue
            implicit none
            real, intent(in):: up(nx,ny,nz,3), nu(nx,ny,nz)
            real, intent(inout):: aj(nx,ny,nz,3), bt(nx,ny,nz,3)
            real, intent(out):: E(nx,ny,nz,3)
            real:: gradPmf(3)
            integer:: i,j,k,m
            real :: a(nx,ny,nz,3), c(nx,ny,nz,3), aa(nx,ny,nz,3), btc(nx,ny,nz,3)
            !real, allocatable, dimension(:,:,:,:) :: a, c, aa, btc
            !allocate(a(nx,ny,nz,3))
            !allocate(c(nx,ny,nz,3))
            !allocate(aa(nx,ny,nz,3))
            !allocate(btc(nx,ny,nz,3))

            call face_to_center(aj,aa)
            
            do i=2,nx-1
                  do j=2,ny-1
                        do k=2,nz-1
                              do m = 1,3
                                    a(i,j,k,m) = aa(i,j,k,m) - up(i,j,k,m)
                                    ue(i,j,k,m) = a(i,j,k,m)
                              enddo
                        enddo
                  enddo
            enddo
            
            call face_to_center(a,aa)
            call edge_to_center(bt,btc)
            call crossf2(a,btc,c)
            
            do i=2,nx-1
                  do j=2,ny-1
                        do k=2,nz-1
                              gradPmf(1) = 0.5*(gradP(i,j,k,1) + gradP(i+1,j,k,1))
                              gradPmf(2) = 0.5*(gradP(i,j,k,2) + gradP(i,j+1,k,2))
                              gradPmf(3) = 0.5*(gradP(i,j,k,3) + gradP(i,j,k+1,3))
                              do m =1,3
                                    E(i,j,k,m) = c(i,j,k,m) + nu(i,j,k)*aj(i,j,k,m) - gradPmf(m)        ! Add in electron pressure
                              enddo
                        enddo
                  enddo
            enddo
            
            call boundary_vector(E)
            call obstacle_boundary_E(E)
      end subroutine get_E
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine predict_B(b12,b1p2,bt,E,aj,up,nu,dtsub)
! Predictor step in magnetic field update
            use dimensions
            use boundary
            implicit none
            real, intent(in):: b12(nx,ny,nz,3),up(nx,ny,nz,3),nu(nx,ny,nz),dtsub
            real, intent(inout):: aj(nx,ny,nz,3), bt(nx,ny,nz,3)
            real, intent(out):: b1p2(nx,ny,nz,3), E(nx,ny,nz,3)
            integer:: i,j,k,m
            real :: curl_E(nx,ny,nz,3)
            !real, allocatable :: curl_E(:,:,:,:)
            !allocate(curl_E(nx,ny,nz,3))

            call get_E(E,bt,aj,up,nu)
            call curlE(E,curl_E)
            
            do i=2,nx-1
                  do j=2,ny-1
                        do k=2,nz-1
                              do m=1,3
                                    b1p2(i,j,k,m) = b12(i,j,k,m) - 2.0*dtsub*curl_E(i,j,k,m)
                              enddo
                        enddo
                  enddo
            enddo
            
            call boundary_vector(b1p2)
            call obstacle_boundary_B(b1p2)
      end subroutine predict_B
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_Ep1(E,b1,b1p2,aj,up,np,nu)
! The main feature here is that E must be calculated at time level
! m + 1/2.  That means that we need B at m + 1/2.  So b1p1 is
! calculated as 0.5*(b1 + b1p2).  uf and np are already at time level
! m + 1/2, so they are used as is.
            use dimensions
            use grid_interp
            use boundary
            use var_arrays, only: gradP
            implicit none
            real, intent(in):: b1(nx,ny,nz,3), b1p2(nx,ny,nz,3), up(nx,ny,nz,3), nu(nx,ny,nz), &
                               np(nx,ny,nz)
            real, intent(out):: E(nx,ny,nz,3), aj(nx,ny,nz,3)
            integer:: i,j,k,m
            real :: gradPmf(3)
            real :: b1p1(nx,ny,nz,3), btc(nx,ny,nz,3), a(nx,ny,nz,3), aa(nx,ny,nz,3), c(nx,ny,nz,3)
            !real, allocatable, dimension(:,:,:,:) :: b1p1, btc, a, aa, c
            !allocate(b1p1(nx,ny,nz,3))
            !allocate(btc(nx,ny,nz,3))
            !allocate(a(nx,ny,nz,3))
            !allocate(aa(nx,ny,nz,3))
            !allocate(c(nx,ny,nz,3))

            do i=1,nx
                  do j=1,ny
                        do k=1,nz
                              b1p1(i,j,k,1) = 0.5*(b1p2(i,j,k,1) + b1(i,j,k,1))
                              b1p1(i,j,k,2) = 0.5*(b1p2(i,j,k,2) + b1(i,j,k,2))
                              b1p1(i,j,k,3) = 0.5*(b1p2(i,j,k,3) + b1(i,j,k,3))
                        enddo
                  enddo
            enddo
            
            call curlB(b1p1,np,aj)
            call obstacle_boundary_B(b1p1)
            call face_to_center(aj,aa)

             do m=1,3
                  do k=2,nz-1
                        do j=2,ny-1
                              do i=2,nx-1
                                    a(i,j,k,m) = aa(i,j,k,m) - up(i,j,k,m)
                              enddo
                        enddo
                  enddo
            enddo
            
            call face_to_center(a,aa)
            call edge_to_center(b1p1,btc)
            
            call crossf2(a,btc,c)
            
            do i=2,nx-1
                  do j=2,ny-1
                        do k=2,nz-1
                              gradPmf(1) = 0.5*(gradP(i,j,k,1) + gradP(i+1,j,k,1))
                              gradPmf(2) = 0.5*(gradP(i,j,k,2) + gradP(i,j+1,k,2))
                              gradPmf(3) = 0.5*(gradP(i,j,k,3) + gradP(i,j,k+1,3))
                              do m=1,3
                                    E(i,j,k,m) = c(i,j,k,m) + nu(i,j,k)*aj(i,j,k,m) - gradPmf(m)
                              enddo
                        enddo
                  enddo
            enddo
            
            call boundary_vector(E)
            call obstacle_boundary_E(E)
            
      end subroutine get_Ep1
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine correct_B(b1,b1p2,E,aj,up,np,nu,dtsub)
! Corrector step in magnetic field update
            use dimensions
            use boundary
            use inputs, only: lww1,lww2
            implicit none
            real, intent(in):: up(nx,ny,nz,3),np(nx,ny,nz),nu(nx,ny,nz),dtsub
            real, intent(inout)::  b1(nx,ny,nz,3),b1p2(nx,ny,nz,3)
            real, intent(out):: E(nx,ny,nz,3),aj(nx,ny,nz,3)
            integer:: i,j,k,m
            real :: curl_E(nx,ny,nz,3)
            !real, allocatable :: curl_E(:,:,:,:)
            !allocate(curl_E(nx,ny,nz,3))

            call get_Ep1(E,b1,b1p2,aj,up,np,nu) !E at time level m
            call curlE(E,curl_E)
            
            do i=2,nx-1
                  do j=2,ny-1
                        do k=2,nz-1
                              do m=1,3
                                    !grid averaging
                                    b1p2(i,j,k,m) = lww1*(b1(i+1,j,k,m)+b1(i-1,j,k,m)+  &
                                    b1(i,j+1,k,m)+b1(i,j-1,k,m)+ &
                                    b1(i,j,k+1,m)+b1(i,j,k-1,m))+ &
                                    lww2*b1(i,j,k,m) - &
                                    dtsub*curl_E(i,j,k,m)
 
                                    !no grid averaging
                                    !b1p2(i,j,k,m) = b1(i,j,k,m) - dtsub*curl_E(i,j,k,m)
                              enddo
                        enddo
                  enddo
            enddo
            
            call boundary_vector(b1p2)
            call obstacle_boundary_B(b1p2)
      end subroutine correct_B
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine check_time_step(bt,np,dtsub,ntf)
            use dimensions
            use grid, only: dz_grid
            use inputs, only: alpha
            implicit none
            real, intent(in):: bt(nx,ny,nz,3), np(nx,ny,nz)
            integer, intent(inout):: ntf
            real, intent(inout):: dtsub
            real:: ak,a1,a2,womega,phi,deltat,btot
            integer:: i,j,k
            
            do i=1,nx
                  do j=1,ny
                        do k=1,nz
                              ak = 2.0/dz_grid(k)
                              btot = sqrt(bt(i,j,k,1)**2 + bt(i,j,k,2)**2 + &
                                    bt(i,j,k,3)**2)
                              a1 = ak**2*Btot/(alpha*(np(i,j,k)))
                              a2 = (ak*Btot)**2/(alpha*(np(i,j,k)))
                              
                              womega = 0.5*(a1 + sqrt(a1**2 + 4*a2))
                              phi = womega/ak
                              deltat = dz_grid(k)/phi

                              if(deltat .le. 2.0*dtsub) then 
                                    write(*,*) 'time stepping error...',i,j,k
                                    dtsub = dtsub/2.0
                                    ntf = ntf*2
                              endif
                        enddo
                  enddo
            enddo
            
      end subroutine check_time_step
      
end module gutsf
