module gutsf
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fieldMove_Boundary(b0,m) !loops run 1 o n since values are only being copied
            use dimensions
            use boundary
            use inputs, only: q, mO,b0_init,pi,ddthickness, ByConeAngle, dt, va, delz
            use grid
            implicit none
            !real, intent(in):: 
            real, intent(inout):: b0(nx,ny,nz,3)!,b1(nx,ny,nz,3),b1p2(nx,ny,nz,3)
            !real, intent(out):: bt(nx,ny,nz,3), b12(nx,ny,nz,3)
            
            integer:: i,j,k,m
            real:: eoverm,dtheta
            integer:: Bsetup, zMove,mPerCell
            
            
            
            Bsetup= 1
            dtheta = 2.0*pi / 4.0 /(2.0*ddthickness)
            
            eoverm = q/mO
            
            	!do i=1,nx
                !  do j=1,ny
                !        do k=1,nz            
		!		if (Bsetup .eq. 1) then !Finite Thickness TD
		!		!for rotating B direction, constant Btot
		!			if (k .le. nz/2.0-ddthickness) then
   		!				 b1(i,j,k,1) = b1(i,j,k,1) + b0_init*eoverm
   		!				 b1(i,j,k,2) = b1(i,j,k,2) + 0.0 
   		!				 b1(i,j,k,3) = b1(i,j,k,3) + 0.0
		!			endif
		!			if ((k .gt. nz/2.0-ddthickness) .and. (k .lt. nz/2.0+ddthickness)) then
   		!	 			 b1(i,j,k,1) = b1(i,j,k,1) + b0_init*eoverm*cos(dtheta*(k-nz/2+ddthickness))
   		!				 b1(i,j,k,2) = b1(i,j,k,2) + -b0_init*eoverm*sin(dtheta*(k-nz/2+ddthickness))
   		!				 b1(i,j,k,3) = b1(i,j,k,3) + 0.0
		!			endif
		!			if (k .ge. nz/2.0+ddthickness) then
    		!				b1(i,j,k,1) = b1(i,j,k,1) + 0.0
    		!				b1(i,j,k,2) = b1(i,j,k,2) + -b0_init*eoverm
    		!				b1(i,j,k,3) = b1(i,j,k,3) + 0.0!b0_init*eoverm
		!			endif
        	!		else
!
 !                             		b1(i,j,k,1) = b0_init*eoverm
  !                           		 b1(i,j,k,2) = 0.0
   !                           		b1(i,j,k,3) = 0.0
    !                          	endif
     !                   enddo
      !            enddo
       !     enddo
            
            
            
            
          mPerCell = floor((qx(2)-qx(1))/(1.0*va*dt))
          zMove = (m/mPerCell)
            
            
            write(*,*) 'm,mPerCell, zMove', m, mPerCell,zMove
            
            
          do i=1,2
                 do j=1,ny
                        do k= 1,nz
                              
                                    
                                    
                                  if (k .le. nz/2.0) then
                              		if (k .le. nz/3.0) then !BL bottom 
   						b0(i,j,k,1) =  (cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm +  (cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( ( qz(nz/3.0)-qz(k+zMove))/(ddthickness*delz))
   						b0(i,j,k,2) =  (sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm +  (sin(ByConeAngle/180.0*pi)-sin((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( ( qz(nz/3.0)-qz(k+zMove))/(ddthickness*delz))
   						b0(i,j,k,3) = 0.0
					endif
					if (k .gt. nz/3.0) then !BL top
    						b0(i,j,k,1) = -(cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm +  (cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( (qz(k+zMove)-qz(nz/3.0) )/(ddthickness*delz))
    						b0(i,j,k,2) =  (sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm + (-sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( (qz(k+zMove)-qz(nz/3.0) )/(ddthickness*delz))
    						b0(i,j,k,3) = 0.0
					endif
				endif
				 
				 if (k .gt. nz/2.0) then			
				 		
				 	if (k .gt. 2.0*nz/3.0) then !BL top
    						b0(i,j,k,1) =  (sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm + (-sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( (qz(k+zMove)-qz(2.0*nz/3.0) )/(ddthickness*delz))
    						b0(i,j,k,2) = (cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm -  (cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( (qz(k+zMove)-qz(2.0*nz/3.0) )/(ddthickness*delz))
    						
    						b0(i,j,k,3) = 0.0
					endif
					
                              		if (k .le. 2.0*nz/3.0) then !BL bottom 
   						b0(i,j,k,1) =  -(sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm -  (sin(ByConeAngle/180.0*pi)-sin((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( ( qz(2.0*nz/3.0)-qz(k+zMove))/(ddthickness*delz))
   						b0(i,j,k,2) =  (cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm +  (cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( ( qz(2.0*nz/3.0)-qz(k+zMove))/(ddthickness*delz))
   						
   						b0(i,j,k,3) = 0.0
					endif
				 endif   
                                    
                                    
                                    
                                    
                              !enddo
                       enddo
                  enddo
            enddo
            
           
            
            
            
   !         do i=1,nx
    !              do j=1,ny
     !                   do k= 1,nz
      !                        bt(i,j,k,1) = b1p2(i,j,k,1)!! + b0(i,j,k,1)
       !                       bt(i,j,k,2) = b1p2(i,j,k,2)!! + b0(i,j,k,2)
        !                      bt(i,j,k,3) = b1p2(i,j,k,3)!! + b0(i,j,k,3) 
         !                     do m=1,3
          !                          b12(i,j,k,m)= b1(i,j,k,m)
           !                         b1(i,j,k,m) = b1p2(i,j,k,m)
            !                  enddo
             !           enddo
              !    enddo
            !enddo
            
            !call fix_normal_b(bt)
            !call fix_normal_b(b1)
            !call fix_normal_b(b12)
            !call fix_normal_b(b1p2)

      end subroutine fieldMove_Boundary
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine f_update_tlev(b1,b12,b1p2,bt,b0) !loops run 1 o n since values are only being copied
            use dimensions
            use boundary
            use inputs, only: q, mO,b0_init,pi,ddthickness!, ByConeAngle, dt, va, delz
            use grid
            implicit none
            real, intent(in):: b0(nx,ny,nz,3)
            real, intent(inout):: b1(nx,ny,nz,3),b1p2(nx,ny,nz,3)
            real, intent(out):: bt(nx,ny,nz,3), b12(nx,ny,nz,3)
            
            integer:: i,j,k,m
            real:: eoverm,dtheta!,mPerCell
            integer:: Bsetup!, zMove
            
            
            
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
            
            !call fix_normal_b(bt)
            !call fix_normal_b(b1)
            !call fix_normal_b(b12)
            !call fix_normal_b(b1p2)

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
            real:: ct(nx,ny,nz,3)
            real:: ax,ay,az,bx,by,bz
            integer:: i,j,k,ip,jp,kp
            
            call boundary_vector(aa)
            call boundary_vector(btc)
!            call periodic(aa)
!            call periodic(btc)
            
            do i=2,nx-1
                  do j=2,ny-1
                        do k=2,nz-1
!                              im = i-1
!                              jm = j-1
!                              km = k-1
                              
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
!            call periodic(ct)
            
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
!            call periodic(cc)
            
      end subroutine crossf2
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine curlB2(b0,b1,np,aj)
! Calculates curl B / n*alpha. The resulting 'current' is called aj which is used in several other placed in the code.
! This curl is performed on the main cell where B is covariant.  The resulting current is main cell contravariant.
! Nnote that dz_cell is used for the cell dimension since dz_grid is not equal to dz_cell on non-uniform grid.

            use dimensions
            use boundary
            use grid, only: dx_cell,dy_cell,dz_cell
            use inputs, only: dx,dy,alpha
            use var_arrays, only: curlBcurrent
            implicit none
            real, intent(in):: np(nx,ny,nz)
            real, intent(inout):: b0(nx,ny,nz,3), b1(nx,ny,nz,3)
            real, intent(out):: aj(nx,ny,nz,3)
            real:: curl_B(3), ntot(3)
            integer:: i,j,k,m,ip,jp,kp
            
            call boundary_vector(b1)
             do i=1,1
                 do j=1,ny
                        do k= 1,nz
                              do m=1,3

                                    !b1(i,j,k,m) = b0(i,j,k,m)
 
                              enddo
                        enddo
                  enddo
             enddo
!            call periodic(b1)
!            call fix_normal_b(b1)
            do i=2,nx-1
                  do j=2,ny-1
                        do k=2,nz-1
                              ip = i+1
                              jp = j+1
                              kp = k+1
                              
                              ntot(1) = 0.5*(np(i,j,k)+np(ip,j,k))
                              ntot(2) = 0.5*(np(i,j,k)+np(i,jp,k))
                              ntot(3) = 0.5*(np(i,j,k)+np(i,j,kp))
                              
                              curl_B(1) = (b1(i,j,k,3)/dy) - (b1(i,j-1,k,3)/dy) &
                                    + (b1(i,j,k-1,2)/dz_cell(k)) &
                                    - (b1(i,j,k,2)/dz_cell(k))
                              curl_B(2) = (b1(i,j,k,1)/dz_cell(k)) &
                                    - (b1(i,j,k-1,1)/dz_cell(k)) &
                                    - (b1(i,j,k,3)/dx) + (b1(i-1,j,k,3)/dx)
                              curl_B(3) = (b1(i,j,k,2)/dx) - (b1(i-1,j,k,2)/dx) &
                                    + (b1(i,j-1,k,1)/dy) - (b1(i,j,k,1)/dy)
                                    
                              do m=1,3
                                    aj(i,j,k,m) = curl_B(m)/(ntot(m)*alpha)
                                    curlBcurrent(i,j,k,m) = curl_B(m)
                              enddo
                        enddo
                  enddo
            enddo
            
      end subroutine curlB2
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine curlB(b0,b1,np,aj)
! Calculates curl B / n*alpha.  The resulting "current" is called aj
! which is used in several other places in the code.  This curl is 
! performed on the main cell where B is covarient.  The resulting
! current is main cell contravarient.  Note that dz_cell is used for
! the cell dimensions since dz_grid is not equal to dz_cell on non-
! uniform grid.
            use dimensions
            use boundary
            use grid, only: dx_cell,dy_cell,dz_cell
            use inputs, only: alpha,ddthickness
            use var_arrays, only: curlBcurrent
            implicit none
            real, intent(in):: np(nx,ny,nz)
            real, intent(inout):: b0(nx,ny,nz,3),b1(nx,ny,nz,3)
            real, intent(out) :: aj(nx,ny,nz,3)
            real:: curl_B(3), ntot(3)
            integer:: i,j,k,m,ip,jp,kp
            
            call boundary_vector(b1)
            do i=1,1
                 do j=1,ny
                        do k= 1,nz
                              do m=1,3

                                    !b1(i,j,k,m) = b0(i,j,k,m)
 
                              enddo
                        enddo
                  enddo
             enddo
!            call periodic(b1)
            !write(*,*) 'b1(x,:,:),b0', b1(1,2,51,2),b0(1,2,51,2)
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
                                    !if ((i.eq.100) .and. (k .eq. 151)) then
                                       !write(*,*) 'bxj bxk', (b1(i,j,k,3) - &
                                    !b1(i,j-1,k,3))/dy_cell(j),(b1(i,j,k-1,2) - b1(i,j,k,2))/dz_cell(k) 
					!write(*,*) 'aj,curl_Bx', aj(100,2,151,1),curl_B(1)
					!write(*,*) 'bj,bj-1,bk-1,bk', b1(i,j,k,3), b1(i,j-1,k,3),b1(i,j,k-1,2),b1(i,j,k,2)
				    !endif
                              enddo
                        enddo
                  enddo
            enddo
                        
      do i=1,nx      
          do j=2,ny-1
                  do k=1,nz
                        do m=1,3
                        
                        	!if (i .ge. nx-1) then
                        	!	if (k .lt. nz/2.0-ddthickness) then
						!aj(i,j,k,m) = 0.0
				!	endif
				!	if (k .gt. nz/2.0-ddthickness) then
						!aj(i,j,k,m) = 0.0
				!	endif
                              	!	
                              	!endif
                              	
                              	!if (i .le. 1) then
                              	!	aj(i,j,k,m) = 0.0
                              	!endif
                              	!if (k .lt. nz/2.0-ddthickness) then
				!	aj(2,j,k,m) = 0.0
				!endif
				!if (k .gt. nz/2.0-ddthickness) then
				!	aj(2,j,k,m) = 0.0
				!endif
                              	
                              	!if (k .le. 2) then
                              	!	aj(i,j,k,m) = 0.0
                              	!endif
                              	
                              	!if (k .ge. nz-1) then
                              	!	aj(i,j,k,m) = 0.0
                              	!endif
                              	
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
            use var_arrays, only: grav, gradP, ue
            use inputs, only: mion
            implicit none
            real, intent(in):: up(nx,ny,nz,3), nu(nx,ny,nz)
            real, intent(inout):: bt(nx,ny,nz,3)
            real, intent(out):: E(nx,ny,nz,3)
            real:: aj(nx,ny,nz,3), a(nx,ny,nz,3), c(nx,ny,nz,3), aa(nx,ny,nz,3), btc(nx,ny,nz,3), gravc(nx,ny,nz), gradPmf(3)
            integer:: i,j,k,m
            
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
            call grav_to_center(grav,gravc)
            
            
            
            do i=2,nx-1
                  do j=2,ny-1
                        do k=2,nz-1
                              gradPmf(1) = 0.5*(gradP(i,j,k,1) + gradP(i+1,j,k,1))
                              gradPmf(2) = 0.5*(gradP(i,j,k,2) + gradP(i,j+1,k,2))
                              gradPmf(3) = 0.5*(gradP(i,j,k,3) + gradP(i,j,k+1,3))
                              do m =1,3
                                    E(i,j,k,m) = c(i,j,k,m) + nu(i,j,k)*aj(i,j,k,m) - gradPmf(m)        ! Add in electron pressure
                              enddo
                                    !E(i,j,k,3) = c(i,j,k,3) + nu(i,j,k)*aj(i,j,k,3) + gravc(i,j,k) - gradPmf(3) ! Add in gravity term and electron pressure
!                                    write(*,*) 'Electric field.................', c(i,j,k,3) + nu(i,j,k)*aj(i,j,k,3)
!                                    write(*,*) 'Gravity field..................', gravc(i,j,k)
                        enddo
                  enddo
            enddo
            
            call boundary_vector(E)
!            call periodic(E)
	call fix_tangential_E(E)
      
      end subroutine get_E
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine predict_B(b0,b12,b1p2,bt,E,aj,up,nu,dtsub)
! Predictor step in magnetic field update
            use dimensions
            use boundary
            implicit none
            real, intent(in):: b0(nx,ny,nz,3), b12(nx,ny,nz,3), aj(nx,ny,nz,3),up(nx,ny,nz,3),nu(nx,ny,nz),dtsub
            real, intent(inout):: bt(nx,ny,nz,3)
            real, intent(out):: b1p2(nx,ny,nz,3), E(nx,ny,nz,3)
            real:: curl_E(nx,ny,nz,3)
            integer:: i,j,k,m
            
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
            do i=1,1
                 do j=1,ny
                        do k= 1,nz
                              do m=1,3

              !                      b1p2(i,j,k,m) = b0(i,j,k,m)
 
                              enddo
                        enddo
                  enddo
             enddo
!            call periodic(b1p2)
            
      end subroutine predict_B
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_Ep1(E,b0,b1,b1p2,aj,up,np,nu)
! The main feature here is that E must be calculated at time level
! m + 1/2.  That means that we need B at m + 1/2.  So b1p1 is
! calculated as 0.5*(b1 + b1p2).  uf and np are already at time level
! m + 1/2, so they are used as is.
            use dimensions
            use grid_interp
            use boundary
            use var_arrays, only: grav, gradP
            use inputs, only: mion
            implicit none
            real, intent(in):: b1(nx,ny,nz,3), b1p2(nx,ny,nz,3), up(nx,ny,nz,3), nu(nx,ny,nz), &
                               np(nx,ny,nz)
            real, intent(out):: E(nx,ny,nz,3), aj(nx,ny,nz,3)
            real,intent(inout):: b0(nx,ny,nz,3)
            real:: b1p1(nx,ny,nz,3), &          !b1 at time level m+1/2
!                   btp1(nx,ny,nz,3), &          !bt at time level m+1/2
!                   btp1mf(nx,ny,nz,3), &        !btp1 at contravariant position
                   btc(nx,ny,nz,3), a(nx,ny,nz,3), aa(nx,ny,nz,3), c(nx,ny,nz,3), gravc(nx,ny,nz), gradPmf(3)
            integer:: i,j,k,m
            
            do i=1,nx
                  do j=1,ny
                        do k=1,nz
!                              btp1(i,j,k,1) = 0.5*(b1p2(i,j,k,1) + b1(i,j,k,1))
                              b1p1(i,j,k,1) = 0.5*(b1p2(i,j,k,1) + b1(i,j,k,1))
!                              btp1(i,j,k,2) = 0.5*(b1p2(i,j,k,2) + b1(i,j,k,2))
                              b1p1(i,j,k,2) = 0.5*(b1p2(i,j,k,2) + b1(i,j,k,2))
!                              btp1(i,j,k,3) = 0.5*(b1p2(i,j,k,3) + b1(i,j,k,3))
                              b1p1(i,j,k,3) = 0.5*(b1p2(i,j,k,3) + b1(i,j,k,3))
                        enddo
                  enddo
            enddo
            
            call curlB(b0,b1p1,np,aj)
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
!            call edge_to_face(btp1,btp1mf)
!            call face_to_center(btp1mf,btc)
            call edge_to_center(b1p1,btc)
            call grav_to_center(grav,gravc)
            
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
                                  
                                  !  E(i,j,k,3) = c(i,j,k,3) + nu(i,j,k)*aj(i,j,k,3) + gravc(i,j,k) - gradPmf(3)! Add in gravity term and electron pressure
!                                    write(*,*) 'Electric field.................', c(i,j,k,3) + nu(i,j,k)*aj(i,j,k,3)
!                                    write(*,*) 'Gravity field..................', gravc(i,j,k)
                        enddo
                  enddo
            enddo
            
            call boundary_vector(E)
            

	    call fix_tangential_E(E)
!            call periodic(E)
            
      end subroutine get_Ep1
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine correct_B(b0,b1,b1p2,E,aj,up,np,nu,dtsub)
! Corrector step in magnetic field update
            use dimensions
            use boundary
            use inputs, only: lww1,lww2
            implicit none
            real, intent(in):: up(nx,ny,nz,3),np(nx,ny,nz),nu(nx,ny,nz),dtsub
            real, intent(inout)::  b0(nx,ny,nz,3), b1(nx,ny,nz,3),b1p2(nx,ny,nz,3)
            real, intent(out):: E(nx,ny,nz,3),aj(nx,ny,nz,3)
            real:: curl_E(nx,ny,nz,3)
            integer:: i,j,k,m
            
            call get_Ep1(E,b0,b1,b1p2,aj,up,np,nu) !E at time level m
            call curlE(E,curl_E)
            do i=1,1
                 do j=1,ny
                        do k= 1,nz
                              do m=1,3

                    !                b1(i,j,k,m) = b0(i,j,k,m)
 
                              enddo
                        enddo
                  enddo
             enddo
            
            do i=2,nx-1
                  do j=2,ny-1
                        do k=2,nz-1
                              do m=1,3
                              !if ((i .ge. nx-3) .or. (i .le. 4) ) then
                              		b1p2(i,j,k,m)=lww1*(b1(i+1,j,k,m)+b1(i-1,j,k,m)+  &
                                          b1(i,j+1,k,m)+b1(i,j-1,k,m)+ &
                                          b1(i,j,k+1,m)+b1(i,j,k-1,m))+ &
                                          lww2*b1(i,j,k,m) - &
                                          dtsub*curl_E(i,j,k,m)
                              !else
                                    !grid averaging
                            		!b1p2(i,j,k,m)=lww1*(b1(i+1,j,k,m)+b1(i-1,j,k,m)+  &
                                        !  b1(i,j+1,k,m)+b1(i,j-1,k,m)+ &
                                        !  b1(i,j,k+1,m)+b1(i,j,k-1,m))+ &
                                        !  1.0*b1(i,j,k,m) - &
                                        !  dtsub*curl_E(i,j,k,m)
                              !            b1p2(i,j,k,m) = b1(i,j,k,m) - dtsub*curl_E(i,j,k,m)
                              !endif
                                          
                                          
                                    !no grid averaging
                                !    b1p2(i,j,k,m) = b1(i,j,k,m) - dtsub*curl_E(i,j,k,m)
                              enddo
                        enddo
                  enddo
            enddo
            
            call boundary_vector(b1p2)
           do i=1,1
                 do j=1,ny
                        do k= 1,nz
                              do m=1,3

                      !              b1p2(i,j,k,m) = b0(i,j,k,m)
 
                              enddo
                        enddo
                  enddo
             enddo
!            call periodic(b1p2)
            
      end subroutine correct_B
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine check_time_step(bt,np,dtsub,ntf)
            use dimensions
            use grid, only: dz_grid
            use inputs, only: dx,alpha
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
                                    ntf = ntf*2.0
                              endif
                        enddo
                  enddo
            enddo
            
      end subroutine check_time_step

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Momentum_diag(up,uf,np,nf,E,b1,pup,puf,peb,input_p)
            use dimensions
            use inputs, only: epsilon,mO,q,mBa
            use grid, only: dx_cell, dy_cell, dz_cell
            implicit none
            real, intent(out):: pup(3),puf(3),peb(3)
            real, intent(in):: up(nx,ny,nz,3), uf(nx,ny,nz,3), np(nx,ny,nz), &
                               nf(nx,ny,nz), input_p(3)
            real, intent(inout):: E(nx,ny,nz,3), b1(nx,ny,nz,3)
            real:: exb(nx,ny,nz,3)                   
            real:: npave(3), nfave(3), vol
            integer:: i,j,k,m,ip,jp,kp
            
            call crossf2(E,b1,exb)
            
            do m=1,3
                  pup(m) = 0
                  puf(m) = 0
                  peb(m) = 0
            enddo
            
            do i = 2,nx-1
                  do j = 2, ny - 1
                        do k = 2, nz-1
                              ip = i+1
                              jp = j+1
                              kp = k+1
                              if (ip .eq. nx) ip = nx-1
                              if (jp .eq. ny) jp = ny-1
                              if (kp .eq. nz) kp = nz-1
                              vol = dx_cell(i)*dy_cell(j)*dz_cell(k)
!                              npave(1) = 0.5*(np(i,j,k) + np(ip,j,k))
!                              npave(2) = 0.5*(np(i,j,k) + np(i,jp,k))
!                              npave(3) = 0.5*(np(i,j,k) + np(i,j,kp))
!                              nfave(1) = 0.5*(np(i,j,k) + nf(ip,j,k))
!                              nfave(2) = 0.5*(np(i,j,k) + nf(i,jp,k))
!                              nfave(3) = 0.5*(np(i,j,k) + nf(i,j,kp))
                              do m = 1,3
                                    pup(m) = pup(m) + np(i,j,k)*vol*mBa*up(i,j,k,m)
                                    puf(m) = puf(m) + nf(i,j,k)*vol*mO*uf(i,j,k,m)
                                    peb(m) = peb(m) + epsilon*1e3*exb(i,j,k,m)*vol*(mO/q)
                              enddo
                              
                        enddo
                  enddo
            enddo
            
      end subroutine Momentum_diag
      
end module gutsf
