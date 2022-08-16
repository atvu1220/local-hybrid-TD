module initial
      use dimensions
      implicit none
      save
      contains
      
      
      subroutine grd6_setup(b0,bt,b12,b1,b1p2,nu,input_Eb)
            use inputs, only: q, mO, PI, b0_top, b0_bottom, b0_init, nu_init, km_to_m, mu0,dx,dy,delz, ddthickness, pi, dx_frac, magneticShear,ByConeAngle, boundx
            use grid, only: dx_cell, dy_cell, dz_cell,qx,qy,qz
            implicit none
            real, intent(out):: b0(nx,ny,nz,3), &
                                bt(nx,ny,nz,3), &
                                b12(nx,ny,nz,3), &
                                b1(nx,ny,nz,3), &
                                b1p2(nx,ny,nz,3), &
                                nu(nx,ny,nz)
            real, intent(inout):: input_Eb                    
                                
            real:: eoverm, mO_q, vol, b0eoverm
            real:: b0_1x, b0_2x, b0_1y, b0_2y, phi, dtheta
            integer:: i,j,k,m,Bsetup
            
            if (ByConeAngle .gt. 1.0) then
            	Bsetup= 15
            	if (boundx .eq. 5) then
            		Bsetup = 17 !2 TD, reflecting, periodic
            	endif
            else	
            	Bsetup= 8
            	if (boundx .eq. 5) then
            		Bsetup = 16 !2 TD, reflecting, periodic
            	endif
            endif
            
            !Bsetup = 18
            
            dtheta = 2*pi / 4 /(2*ddthickness)
            
            eoverm = q/mO
            mO_q = mO/q
            b0eoverm=b0_init*eoverm
            phi = 2.0*PI/180.0
            	do i=1,nx
                  do j=1,ny
                        do k=1,nz            
				if (Bsetup .eq. 1) then !Finite Thickness TD
				!for rotating B direction, constant Btot
					if (k .le. nz/2.0-ddthickness) then
   						 b0(i,j,k,1) = b0_init*eoverm
   						 b0(i,j,k,2) = 0.0
   						 b0(i,j,k,3) = 0.0
					endif
					if ((k .gt. nz/2.0-ddthickness) .and. (k .lt. nz/2.0+ddthickness)) then
   			 			b0(i,j,k,1) = b0_init*eoverm*cos(dtheta*(k-nz/2+ddthickness))
   						 b0(i,j,k,2) = b0_init*eoverm*sin(dtheta*(k-nz/2+ddthickness))
   						 b0(i,j,k,3) = 0.0
					endif
					if (k .ge. nz/2.0+ddthickness) then
    						b0(i,j,k,1) = 0.0
    						b0(i,j,k,2) = b0_init*eoverm
    						b0(i,j,k,3) = 0.0!b0_init*eoverm
					endif
        			endif
        			if (Bsetup .eq. 2) then
					if (k .le. nz/2.0) then
   						 b0(i,j,k,1) = b0_init*eoverm!0.0
   						 b0(i,j,k,2) = 0.0!-b0_init*eoverm
   						 b0(i,j,k,3) = 0.0
					endif
					if (k .ge. nz/2.0) then
    						b0(i,j,k,1) = 0.0
    						b0(i,j,k,2) = b0_init*eoverm
    						b0(i,j,k,3) = 0.0!b0_init*eoverm
					endif
                              		!b0(i,j,k,1) = b0_init*eoverm
                             		 !b0(i,j,k,2) = 0.0
                              		!b0(i,j,k,3) = 0.0
                              	endif
                              	if (Bsetup .eq. 3) then
                              		if (k .le. nz/2.0) then
   						b0(i,j,k,1) = b0_init*eoverm*tanh( ( qz(nz/2.0)-qz(k))/(ddthickness*delz))
   						b0(i,j,k,2) = 0.0
   						b0(i,j,k,3) = 0.0
					endif
					if (k .ge. nz/2.0) then
    						b0(i,j,k,1) = 0.0
    						b0(i,j,k,2) = b0_init*eoverm*tanh( (qz(k)-qz(nz/2.0) )/(ddthickness*delz))
    						b0(i,j,k,3) = 0.0
					endif
                              	endif
                              	
                               if (Bsetup .eq. 4) then !RD
                              		if (i .le. 2*nx/10.0) then
   						b0(i,j,k,1) = b0_init*eoverm/(sqrt(2.0))
   						b0(i,j,k,2) = b0_init*eoverm/(sqrt(2.0))
   						b0(i,j,k,3) = 0.0
					endif
					if (i .ge. 2*nx/10.0) then
    						b0(i,j,k,1) = b0_init*eoverm
    						b0(i,j,k,2) = 0.0
    						b0(i,j,k,3) = 0.0!b0_init*eoverm
					endif
                              	endif
                              if (Bsetup .eq. 5) then !BLMN Coordinates, +y
                              !Constant B everywhere, BM
                              b0(i,j,k,1) = b0_init*eoverm*0.25
                              b0(i,j,k,2) = b0_init*eoverm*0.25
                              b0(i,j,k,3) = 0.0
                              		if (k .le. nz/2.0) then !BL bottom
   						b0(i,j,k,1) = b0(i,j,k,1) + b0_init*eoverm*0.25*tanh( ( qz(nz/2.0)-qz(k))/(ddthickness*delz))
   						b0(i,j,k,2) = b0(i,j,k,2) - b0_init*eoverm*0.25*tanh( ( qz(nz/2.0)-qz(k))/(ddthickness*delz))
   						b0(i,j,k,3) = 0.0
					endif
					if (k .gt. nz/2.0) then !BL top
    						b0(i,j,k,1) = b0(i,j,k,1) - b0_init*eoverm*0.25*tanh( (qz(k)-qz(nz/2.0) )/(ddthickness*delz))
    						b0(i,j,k,2) = b0(i,j,k,2) + b0_init*eoverm*0.25*tanh( (qz(k)-qz(nz/2.0) )/(ddthickness*delz))
    						b0(i,j,k,3) = 0.0
					endif
					!write(*,*) 'b0_init, eoverm', b0_init, eoverm, b0_init*eoverm
					!write(*,*) 'i,j,k',i,j,k,b0(i,j,k,1),b0(i,j,k,2)
                              endif
                              if (Bsetup .eq. 6) then !BLMN Coordinates, but -y 
                              !Constant B everywhere, BM
                              b0(i,j,k,1) = b0_init*eoverm*0.25
                              b0(i,j,k,2) = -b0_init*eoverm*0.25
                              b0(i,j,k,3) = 0.0
                              		if (k .le. nz/2.0) then !BL bottom
   						b0(i,j,k,1) = b0(i,j,k,1) + b0_init*eoverm*0.25*tanh( ( qz(nz/2.0)-qz(k))/(ddthickness*delz))
   						b0(i,j,k,2) = b0(i,j,k,2) + b0_init*eoverm*0.25*tanh( ( qz(nz/2.0)-qz(k))/(ddthickness*delz))
   						b0(i,j,k,3) = 0.0
					endif
					if (k .gt. nz/2.0) then !BL top
    						b0(i,j,k,1) = b0(i,j,k,1) - b0_init*eoverm*0.25*tanh( (qz(k)-qz(nz/2.0) )/(ddthickness*delz))
    						b0(i,j,k,2) = b0(i,j,k,2) - b0_init*eoverm*0.25*tanh( (qz(k)-qz(nz/2.0) )/(ddthickness*delz))
    						b0(i,j,k,3) = 0.0
					endif
                              endif
                              if (Bsetup .eq. 7) then !BLMN Coordinates
                              !Constant B everywhere, BM
                              b0(i,j,k,1) = b0_init*eoverm*0.25*(1.0+0.5*sqrt(2.0))!0.25
                              b0(i,j,k,2) = b0_init*eoverm*0.5*(1.0/4.0)*sqrt(2.0)!0.25
                              b0(i,j,k,3) = 0.0
                              		if (k .le. nz/2.0) then !BL bottom
   						b0(i,j,k,1) = b0(i,j,k,1) + 0.25*(1.0-0.5*sqrt(2.0))*b0_init*eoverm*tanh( ( qz(nz/2.0)-qz(k))/(ddthickness*delz))
   						b0(i,j,k,2) = b0(i,j,k,2) - 0.5*(1.0/4.0)*sqrt(2.0)*b0_init*eoverm*tanh( ( qz(nz/2.0)-qz(k))/(ddthickness*delz))
   						b0(i,j,k,3) = 0.0
					endif
					if (k .gt. nz/2.0) then !BL top
    						b0(i,j,k,1) = b0(i,j,k,1) - 0.25*(1.0-0.5*sqrt(2.0))*b0_init*eoverm*tanh( (qz(k)-qz(nz/2.0) )/(ddthickness*delz))
    						b0(i,j,k,2) = b0(i,j,k,2) + 0.5*(1.0/4.0)*sqrt(2.0)*b0_init*eoverm*tanh( (qz(k)-qz(nz/2.0) )/(ddthickness*delz))
    						b0(i,j,k,3) = 0.0
					endif
                              endif
                               if (Bsetup .eq. 8) then !BLMN Coordinates, with variable magneticShear from input.dat
                              !Constant B everywhere, BM
                              b0(i,j,k,1) = b0_init*eoverm*0.5*cos(magneticShear/180*pi/2)
                              b0(i,j,k,2) = -(b0_init*eoverm*0.5*sin(magneticShear/180*pi/2))
                              b0(i,j,k,3) = 0.0
                              		if (k .le. nz/2.0) then !BL bottom
   						b0(i,j,k,1) = b0(i,j,k,1) + (cos(0.0) - cos(magneticShear/180.0*pi/2.0) ) * b0_init*eoverm*0.5*tanh( ( qz(nz/2.0)-qz(k))/(ddthickness*delz))
   						b0(i,j,k,2) = (b0(i,j,k,2) + (         sin(magneticShear/180.0*pi/2.0) ) * b0_init*eoverm*0.5*tanh( ( qz(nz/2.0)-qz(k))/(ddthickness*delz)))
   						b0(i,j,k,3) = 0.0
					endif
					if (k .gt. nz/2.0) then !BL top
    						b0(i,j,k,1) = b0(i,j,k,1) - (cos(magneticShear/180.0*pi/2) - cos(magneticShear/180.0*pi)   ) * b0_init*eoverm*0.5*tanh( (qz(k)-qz(nz/2.0) )/(ddthickness*delz))
    						b0(i,j,k,2) = (b0(i,j,k,2) -(sin(magneticShear/180.0*pi)   - sin(magneticShear/180.0*pi/2.0) ) * b0_init*eoverm*0.5*tanh( (qz(k)-qz(nz/2.0) )/(ddthickness*delz)))
    						b0(i,j,k,3) = 0.0
					endif
					!write(*,*) 'b0_init, eoverm', b0_init, eoverm, b0_init*eoverm
					!write(*,*) 'i,j,k',i,j,k,b0(i,j,k,1),b0(i,j,k,2)
                              endif
                              if (Bsetup .eq. 9) then !BLMN Coordinates, with variable magneticShear frin input.dat and made force-free by Bz, 7/16/2021 "Kinetic model of force-free current sheet Kolotkov et al. 2015"
                              !Constant B everywhere, BM
                              b0(i,j,k,1) = b0_init*eoverm*0.5*cos(magneticShear/180*pi/2)
                              b0(i,j,k,2) = b0_init*eoverm*0.5*sin(magneticShear/180*pi/2)
                              b0(i,j,k,3) = 0.0
                              		if (k .le. nz/2.0) then !BL bottom
   						b0(i,j,k,1) = b0(i,j,k,1) + (cos(0.0) - cos(magneticShear/180.0*pi/2.0) ) * b0_init*eoverm*0.5*tanh( ( qz(nz/2.0)-qz(k))/(ddthickness*delz))
   						b0(i,j,k,2) = b0(i,j,k,2) - (         sin(magneticShear/180.0*pi/2.0) ) * b0_init*eoverm*0.5*tanh( ( qz(nz/2.0)-qz(k))/(ddthickness*delz))
   						b0(i,j,k,3) = sqrt(2.0)/2.0*b0_init*eoverm*1.0/(cosh( ( qz(nz/2.0)-qz(k) ) / (ddthickness*delz) ) )
					endif
					if (k .gt. nz/2.0) then !BL top
    						b0(i,j,k,1) = b0(i,j,k,1) - (cos(magneticShear/180.0*pi/2) - cos(magneticShear/180.0*pi)   ) * b0_init*eoverm*0.5*tanh( (qz(k)-qz(nz/2.0) )/(ddthickness*delz))
    						b0(i,j,k,2) = b0(i,j,k,2) + (sin(magneticShear/180.0*pi)   - sin(magneticShear/180.0*pi/2.0) ) * b0_init*eoverm*0.5*tanh( (qz(k)-qz(nz/2.0) )/(ddthickness*delz))
    						b0(i,j,k,3) = sqrt(2.0)/2.0*b0_init*eoverm*1.0/(cosh( (qz(k)-qz(nz/2.0) )/(ddthickness*delz)))
					endif
					!write(*,*) 'b0_init, eoverm', b0_init, eoverm, b0_init*eoverm
					!write(*,*) 'i,j,k',i,j,k,b0(i,j,k,1),b0(i,j,k,2)
                              endif
                              if (Bsetup .eq. 10) then !BLMN Coordinates, +y, Force-Free with z
                              !Constant B everywhere, BM
                              b0(i,j,k,1) = b0_init*eoverm*0.25
                              b0(i,j,k,2) = b0_init*eoverm*0.25
                              b0(i,j,k,3) = 0.0
                              		if (k .le. nz/2.0) then !BL bottom
   						b0(i,j,k,1) = b0(i,j,k,1) + b0_init*eoverm*0.25*tanh( ( qz(nz/2.0)-qz(k))/(ddthickness*delz))
   						b0(i,j,k,2) = b0(i,j,k,2) - b0_init*eoverm*0.25*tanh( ( qz(nz/2.0)-qz(k))/(ddthickness*delz))
   						b0(i,j,k,3) = sqrt(2.0)/2.0*b0_init*eoverm*0.5*1.0/(cosh( ( qz(nz/2.0)-qz(k))/(ddthickness*delz)))
					endif
					if (k .gt. nz/2.0) then !BL top
    						b0(i,j,k,1) = b0(i,j,k,1) - b0_init*eoverm*0.25*tanh( (qz(k)-qz(nz/2.0) )/(ddthickness*delz))
    						b0(i,j,k,2) = b0(i,j,k,2) + b0_init*eoverm*0.25*tanh( (qz(k)-qz(nz/2.0) )/(ddthickness*delz))
    						b0(i,j,k,3) = sqrt(2.0)/2.0*b0_init*eoverm*0.5*1.0/(cosh( (qz(k)-qz(nz/2.0) )/(ddthickness*delz)))
					endif
					!write(*,*) 'b0_init, eoverm', b0_init, eoverm, b0_init*eoverm
					!write(*,*) 'i,j,k',i,j,k,b0(i,j,k,1),b0(i,j,k,2)
                              endif
                              if (Bsetup .eq. 12) then !BLMN Coordinates, +y 30 degrees, +90 rotation
                              !Constant B everywhere, BM
                              !b0(i,j,k,1) = 0.0
                              !b0(i,j,k,2) = 0.0
                              !b0(i,j,k,3) = 0.0
                              		if (k .le. nz/2.0) then !BL bottom
   						b0(i,j,k,1) =  (cos(30.0/180.0*pi)-cos(120.0/180.0*pi))*0.25*b0_init*eoverm +  (cos(30.0/180.0*pi)+cos(120.0/180.0*pi))*b0_init*eoverm*0.25*tanh( ( qz(nz/2.0)-qz(k))/(ddthickness*delz))
   						b0(i,j,k,2) =  (sin(30.0/180.0*pi)+sin(120.0/180.0*pi))*0.25*b0_init*eoverm +  (sin(30.0/180.0*pi)-sin(120.0/180.0*pi))*b0_init*eoverm*0.25*tanh( ( qz(nz/2.0)-qz(k))/(ddthickness*delz))
   						b0(i,j,k,3) = 0.0
					endif
					if (k .gt. nz/2.0) then !BL top
    						b0(i,j,k,1) = -(cos(30.0/180.0*pi)-cos(120.0/180.0*pi))*0.25*b0_init*eoverm +  (cos(30.0/180.0*pi)+cos(120.0/180.0*pi))*b0_init*eoverm*0.25*tanh( (qz(k)-qz(nz/2.0) )/(ddthickness*delz))
    						b0(i,j,k,2) =  (sin(30.0/180.0*pi)+sin(120.0/180.0*pi))*0.25*b0_init*eoverm + (-sin(30.0/180.0*pi)+sin(120.0/180.0*pi))*b0_init*eoverm*0.25*tanh( (qz(k)-qz(nz/2.0) )/(ddthickness*delz))
    						b0(i,j,k,3) = 0.0
					endif
					!write(*,*) 'b0_init, eoverm', b0_init, eoverm, b0_init*eoverm
					!write(*,*) 'i,j,k',i,j,k,b0(i,j,k,1),b0(i,j,k,2)
                              endif
                              if (Bsetup .eq. 13) then !BLMN Coordinates, -y -30 degrees, +90 rotation
                              !Constant B everywhere, BM
                              !b0(i,j,k,1) = 0.0
                              !b0(i,j,k,2) = 0.0
                              !b0(i,j,k,3) = 0.0
                              		if (k .le. nz/2.0) then !BL bottom
   						b0(i,j,k,1) =  (cos(-30.0/180.0*pi)-cos(60.0/180.0*pi))*0.25*b0_init*eoverm +  (cos(-30.0/180.0*pi)+cos(60.0/180.0*pi))*b0_init*eoverm*0.5*tanh( ( qz(nz/2.0)-qz(k))/(ddthickness*delz))
   						b0(i,j,k,2) =  (sin(-30.0/180.0*pi)+sin(60.0/180.0*pi))*0.25*b0_init*eoverm +  (sin(-30.0/180.0*pi)-sin(60.0/180.0*pi))*b0_init*eoverm*0.5*tanh( ( qz(nz/2.0)-qz(k))/(ddthickness*delz))
   						b0(i,j,k,3) = 0.0
					endif
					if (k .gt. nz/2.0) then !BL top
    						b0(i,j,k,1) = -(cos(-30.0/180.0*pi)-cos(60.0/180.0*pi))*0.25*b0_init*eoverm +  (cos(-30.0/180.0*pi)+cos(60.0/180.0*pi))*b0_init*eoverm*0.25*tanh( (qz(k)-qz(nz/2.0) )/(ddthickness*delz))
    						b0(i,j,k,2) =  (sin(-30.0/180.0*pi)+sin(60.0/180.0*pi))*0.25*b0_init*eoverm + (-sin(-30.0/180.0*pi)+sin(60.0/180.0*pi))*b0_init*eoverm*0.25*tanh( (qz(k)-qz(nz/2.0) )/(ddthickness*delz))
    						b0(i,j,k,3) = 0.0
					endif
					!write(*,*) 'b0_init, eoverm', b0_init, eoverm, b0_init*eoverm
					!write(*,*) 'i,j,k',i,j,k,b0(i,j,k,1),b0(i,j,k,2)
                              endif
                              if (Bsetup .eq. 14) then !BLMN Coordinates, -y -30 degrees, -90 rotation
                              !Constant B everywhere, BM
                              !b0(i,j,k,1) = 0.0
                              !b0(i,j,k,2) = 0.0
                              !b0(i,j,k,3) = 0.0
                              		if (k .le. nz/2.0) then !BL bottom
   						b0(i,j,k,1) =  (cos(-30.0/180.0*pi)-cos(-120.0/180.0*pi))*0.25*b0_init*eoverm +  (cos(-30.0/180.0*pi)+cos(-120.0/180.0*pi))*b0_init*eoverm*0.5*tanh( ( qz(nz/2.0)-qz(k))/(ddthickness*delz))
   						b0(i,j,k,2) =  (sin(-30.0/180.0*pi)+sin(-120.0/180.0*pi))*0.25*b0_init*eoverm +  (sin(-30.0/180.0*pi)-sin(-120.0/180.0*pi))*b0_init*eoverm*0.5*tanh( ( qz(nz/2.0)-qz(k))/(ddthickness*delz))
   						b0(i,j,k,3) = 0.0
					endif
					if (k .gt. nz/2.0) then !BL top
    						b0(i,j,k,1) = -(cos(-30.0/180.0*pi)-cos(-120.0/180.0*pi))*0.25*b0_init*eoverm +  (cos(-30.0/180.0*pi)+cos(-120.0/180.0*pi))*b0_init*eoverm*0.25*tanh( (qz(k)-qz(nz/2.0) )/(ddthickness*delz))
    						b0(i,j,k,2) =  (sin(-30.0/180.0*pi)+sin(-120.0/180.0*pi))*0.25*b0_init*eoverm + (-sin(-30.0/180.0*pi)+sin(-120.0/180.0*pi))*b0_init*eoverm*0.25*tanh( (qz(k)-qz(nz/2.0) )/(ddthickness*delz))
    						b0(i,j,k,3) = 0.0
					endif
					!write(*,*) 'b0_init, eoverm', b0_init, eoverm, b0_init*eoverm
					!write(*,*) 'i,j,k',i,j,k,b0(i,j,k,1),b0(i,j,k,2)
                              endif
                              if (Bsetup .eq. 15) then !BLMN Coordinates,+y 60 degrees 9/28, corrected coefficients.
                              !Constant B everywhere, BM
                              !b0(i,j,k,1) = 0.0
                              !b0(i,j,k,2) = 0.0
                              !b0(i,j,k,3) = 0.0
                              		if (k .le. nz/2.0) then !BL bottom 
   						b0(i,j,k,1) =  (cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+90)/180.0*pi))*0.25*b0_init*eoverm +  (cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.25*tanh( ( qz(nz/2.0)-qz(k))/(ddthickness*delz))
   						b0(i,j,k,2) =  (sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*0.25*b0_init*eoverm +  (sin(ByConeAngle/180.0*pi)-sin((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.25*tanh( ( qz(nz/2.0)-qz(k))/(ddthickness*delz))
   						b0(i,j,k,3) = 0.0
					endif
					if (k .gt. nz/2.0) then !BL top
    						b0(i,j,k,1) = -(cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+90)/180.0*pi))*0.25*b0_init*eoverm +  (cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.25*tanh( (qz(k)-qz(nz/2.0) )/(ddthickness*delz))
    						b0(i,j,k,2) =  (sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*0.25*b0_init*eoverm + (-sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.25*tanh( (qz(k)-qz(nz/2.0) )/(ddthickness*delz))
    						b0(i,j,k,3) = 0.0
					endif
					!write(*,*) 'b0_init, eoverm', b0_init, eoverm, b0_init*eoverm
					!write(*,*) 'i,j,k',i,j,k,b0(i,j,k,1),b0(i,j,k,2)
                              endif
                              
                              
                              
                              if (Bsetup .eq. 16) then !BLMN Coordinates, with variable magneticShear from input.dat// Periodic boundary along z, reflecting x boundary. 2 TDs
                              !Constant B everywhere, BM
                              b0(i,j,k,1) = b0_init*eoverm*0.5*cos(magneticShear/180*pi/2)
                              b0(i,j,k,2) = b0_init*eoverm*0.5*sin(magneticShear/180*pi/2)
                              b0(i,j,k,3) = 0.0
                              
                              
                              if (k .le. nz/2.0) then
                              
                              		if (k .le. nz/3.0) then !BL bottom
   						b0(i,j,k,1) = b0(i,j,k,1) + (cos(0.0) - cos(magneticShear/180.0*pi/2.0) ) * b0_init*eoverm*0.5*tanh( ( qz(nz/3.0)-qz(k))/(ddthickness*delz))
   						b0(i,j,k,2) = b0(i,j,k,2) - (         sin(magneticShear/180.0*pi/2.0) ) * b0_init*eoverm*0.5*tanh( ( qz(nz/3.0)-qz(k))/(ddthickness*delz))
   						b0(i,j,k,3) = 0.0
					endif
					if (k .gt. nz/3.0) then !BL top
    						b0(i,j,k,1) = b0(i,j,k,1) - (cos(magneticShear/180.0*pi/2) - cos(magneticShear/180.0*pi)   ) * b0_init*eoverm*0.5*tanh( (qz(k)-qz(nz/3.0) )/(ddthickness*delz))
    						b0(i,j,k,2) = b0(i,j,k,2) + (sin(magneticShear/180.0*pi)   - sin(magneticShear/180.0*pi/2.0) ) * b0_init*eoverm*0.5*tanh( (qz(k)-qz(nz/3.0) )/(ddthickness*delz))
    						b0(i,j,k,3) = 0.0
					endif
			    endif
			    
			    if (k .gt. nz/2.0) then
			    
                              		if (k .ge. 2.0*nz/3.0) then !BL top
    						b0(i,j,k,1) = b0(i,j,k,1) + (sin(magneticShear/180.0*pi)   - sin(magneticShear/180.0*pi/2.0) ) * b0_init*eoverm*0.5*tanh( (qz(k)-qz(2.0*nz/3.0) )/(ddthickness*delz))
    						b0(i,j,k,2) = b0(i,j,k,2) - (cos(magneticShear/180.0*pi/2) - cos(magneticShear/180.0*pi)   ) * b0_init*eoverm*0.5*tanh( (qz(k)-qz(2.0*nz/3.0) )/(ddthickness*delz))
    						b0(i,j,k,3) = 0.0
					endif
					
					
                              		if (k .lt. 2.0*nz/3.0) then !BL bottom
   						b0(i,j,k,1) = b0(i,j,k,1) - (         sin(magneticShear/180.0*pi/2.0) ) * b0_init*eoverm*0.5*tanh( ( qz(2.0*nz/3.0)-qz(k))/(ddthickness*delz))
   						b0(i,j,k,2) = b0(i,j,k,2) + (cos(0.0) - cos(magneticShear/180.0*pi/2.0) ) * b0_init*eoverm*0.5*tanh( ( qz(2.0*nz/3.0)-qz(k))/(ddthickness*delz))
   						b0(i,j,k,3) = 0.0
					endif

			    endif
			    	
					
					
					
					!write(*,*) 'b0_init, eoverm', b0_init, eoverm, b0_init*eoverm
					!write(*,*) 'i,j,k',i,j,k,b0(i,j,k,1),b0(i,j,k,2)
                            endif
                            
                            if (Bsetup .eq. 17) then !BLMN Coordinates,+y 60 degrees 9/28, corrected coefficients. peridoic z, 2 TD, reflecting boundary, (changed coeff from 0.25 to 0.5 , 1/29/2022)

				if (k .le. nz/2.0) then
                              		if (k .le. nz/3.0) then !BL bottom 
   						b0(i,j,k,1) =  (cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm +  (cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( ( qz(nz/3.0)-qz(k))/(ddthickness*delz))
   						b0(i,j,k,2) =  (sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm +  (sin(ByConeAngle/180.0*pi)-sin((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( ( qz(nz/3.0)-qz(k))/(ddthickness*delz))
   						b0(i,j,k,3) = 0.0
					endif
					if (k .gt. nz/3.0) then !BL top
    						b0(i,j,k,1) = +(cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm -  (cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( (qz(k)-qz(nz/3.0) )/(ddthickness*delz))
    						b0(i,j,k,2) =  (sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm + (-sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( (qz(k)-qz(nz/3.0) )/(ddthickness*delz))
    						b0(i,j,k,3) = 0.0
					endif
				endif
				 
				 if (k .gt. nz/2.0) then			
				 		
				 	if (k .gt. 2.0*nz/3.0) then !BL top
    						b0(i,j,k,1) =  (-sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm + (sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( (qz(k)-qz(2.0*nz/3.0) )/(ddthickness*delz))
    						b0(i,j,k,2) = (cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm -  (cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( (qz(k)-qz(2.0*nz/3.0) )/(ddthickness*delz))
    						
    						b0(i,j,k,3) = 0.0
					endif
					
                              		if (k .le. 2.0*nz/3.0) then !BL bottom 
   						b0(i,j,k,1) =  (-sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm -  (sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( ( qz(2.0*nz/3.0)-qz(k))/(ddthickness*delz))
   						b0(i,j,k,2) =  (cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm +  (cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( ( qz(2.0*nz/3.0)-qz(k))/(ddthickness*delz))
   						
   						b0(i,j,k,3) = 0.0
					endif
				 endif 
                            endif
                            
                            
                            if (Bsetup .eq. 18) then !BLMN Coordinates,+y 60 degrees, peridoic z, 2 TD, reflecting boundary, initially quasi-PERP, time vary quasi-parallel with +z SW later. 
				if (i .le. 50) then !(1*exp(-(qx(i))**2/(50.0/2.0*dx)**2))
				if (k .le. nz/2.0) then
				
                              		if (k .gt. nz/3.0) then !BL bottom 
   						b0(i,j,k,1) =  (cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm +  (-cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( (qz(k)-qz(nz/3.0) )/(ddthickness*delz))
   						b0(i,j,k,2) =  (sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm +  (sin(ByConeAngle/180.0*pi)-sin((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( (qz(k)-qz(nz/3.0) )/(ddthickness*delz))
   						b0(i,j,k,3) = 0.0
					endif
					if (k .le. nz/3.0) then !BL top
    						b0(i,j,k,1) = -(cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm +  (cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( ( qz(nz/3.0)-qz(k))/(ddthickness*delz))
    						b0(i,j,k,2) =  (sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm + (-sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( ( qz(nz/3.0)-qz(k))/(ddthickness*delz))
    						b0(i,j,k,3) = 0.0
					endif
				endif
				 
				 if (k .gt. nz/2.0) then			
				 		
				 	if (k .le. 2.0*nz/3.0) then !BL top
    						b0(i,j,k,1) =  (sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm + (-sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( ( qz(2.0*nz/3.0)-qz(k))/(ddthickness*delz))
    						b0(i,j,k,2) = (cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm -  (cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( ( qz(2.0*nz/3.0)-qz(k))/(ddthickness*delz))
    						
    						b0(i,j,k,3) = 0.0
					endif
					
                              		if (k .gt. 2.0*nz/3.0) then !BL bottom 
   						b0(i,j,k,1) =  -(sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm -  (sin(ByConeAngle/180.0*pi)-sin((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( (qz(k)-qz(2.0*nz/3.0) )/(ddthickness*delz))
   						b0(i,j,k,2) =   (cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm +  (cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( (qz(k)-qz(2.0*nz/3.0) )/(ddthickness*delz))
   						
   						b0(i,j,k,3) = 0.0
					endif
				 endif 
				
				else
				if (k .le. nz/2.0) then
                              		if (k .le. nz/3.0) then !BL bottom 
   						b0(i,j,k,1) =  -(cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm +  (cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( ( qz(nz/3.0)-qz(k))/(ddthickness*delz))
   						
   						!(cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+90)/180.0*pi))*0.25*b0_init*eoverm +  (cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.25*
   						
   						b0(i,j,k,2) =   (sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm + (-sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( ( qz(nz/3.0)-qz(k))/(ddthickness*delz))
   						!(sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*0.25*b0_init*eoverm +  (sin(ByConeAngle/180.0*pi)-sin((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.25*
   						
   						b0(i,j,k,3) = 0.0
					endif
					if (k .gt. nz/3.0) then !BL top
    						b0(i,j,k,1) = -(cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm +  (cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( ( qz(nz/2.0)-qz(k))/(ddthickness*delz))
    						b0(i,j,k,2) =  (sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm + (-sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( ( qz(nz/2.0)-qz(k))/(ddthickness*delz))
    						b0(i,j,k,3) = 0.0
					endif
				endif
				 
				 if (k .gt. nz/2.0) then			
				 		
				 	if (k .gt. 2.0*nz/3.0) then !BL top
    						b0(i,j,k,1) =  -(sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm -  (sin(ByConeAngle/180.0*pi)-sin((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( (qz(k)-qz(2.0*nz/3.0) )/(ddthickness*delz))
    						!(sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*0.25*b0_init*eoverm + (-sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.25*
    						b0(i,j,k,2) =  +(cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm +  (cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( (qz(k)-qz(2.0*nz/3.0) )/(ddthickness*delz))
    						!(cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+90)/180.0*pi))*0.25*b0_init*eoverm -  (cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.25*
    						
    						b0(i,j,k,3) = 0.0
					endif
					
                              		if (k .le. 2.0*nz/3.0) then !BL bottom 
   						b0(i,j,k,1) =  +(-sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm -  (sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( (qz(k)-qz(nz/2.0) )/(ddthickness*delz))
   						b0(i,j,k,2) =  +( cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm +  (cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( (qz(k)-qz(nz/2.0) )/(ddthickness*delz))
   						
   						b0(i,j,k,3) = 0.0
					endif
				 endif
				 endif
                            endif

    				!Setup for testing TD in plane shock vertical Bz
    				 if (Bsetup .eq. 19) then !BLMN Coordinates,+y 60 degrees 9/28, corrected coefficients. peridoic z, 2 TD, reflecting boundary, (changed coeff from 0.25 to 0.5 , 1/29/2022)

                              		if (i .le. nx/2.0) then !BL bottom 
   						b0(i,j,k,1) =  0.0
   						b0(i,j,k,2) =  (sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm +  (sin(ByConeAngle/180.0*pi)-sin((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( ( qx(nx/2.0)-qx(i))/(ddthickness*delz))
   						b0(i,j,k,3) = (cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm +  (cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( ( qx(nx/2.0)-qx(i))/(ddthickness*delz))
					endif
					if (i .gt. nx/2.0) then !BL top
    						b0(i,j,k,1) =  0.0
    						b0(i,j,k,2) =  (sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm + (-sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( (qx(i)-qx(nx/2.0) )/(ddthickness*delz))
    						b0(i,j,k,3) = -(cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm +  (cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( (qx(i)-qx(nx/2.0) )/(ddthickness*delz))
					endif
				endif
    				
    				!Setup for testing TD in plane shock
    				 if (Bsetup .eq. 20) then !BLMN Coordinates,+y 60 degrees 9/28, corrected coefficients. peridoic z, 2 TD, reflecting boundary, (changed coeff from 0.25 to 0.5 , 1/29/2022)

                              		if (k .le. nz/2.0) then !BL bottom 
   						b0(i,j,k,1) =  (cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm +  (cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( ( qz(nz/2.0)-qz(k))/(ddthickness*delz))
   						b0(i,j,k,2) =  (sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm +  (sin(ByConeAngle/180.0*pi)-sin((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( ( qz(nz/2.0)-qz(k))/(ddthickness*delz))
   						b0(i,j,k,3) = 0.0
					endif
					if (k .gt. nz/2.0) then !BL top
    						b0(i,j,k,1) =  -(cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm +  (cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( (qz(k)-qz(nz/2.0) )/(ddthickness*delz))
    						b0(i,j,k,2) =  (sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*0.5*b0_init*eoverm + (-sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*b0_init*eoverm*0.5*tanh( (qz(k)-qz(nz/2.0) )/(ddthickness*delz))
    						b0(i,j,k,3) = 0.0
					endif
				endif
    						
    						
    						
                              
                        enddo
                  enddo
            enddo
        !write(*,*) 'during initialization b1(x,:,:),b0', b1(1,2,51,2),b0(1,2,51,2)
!            input_Eb = 0.0
            do i=1,nx
                  do j=1,ny
                        do k= 1,nz
                        	!nu(i,j,k) = 2*nu_init*b0_init*eoverm*&
                        	!( ( exp(-(qx(nx)-qx(i)))**2/(2*dx)**2 ) + (exp(+(qx(1)-qx(i)))**2/(2*dx)**2) ) + nu_init
                        
                                !nu(i,j,k) = 2*nu_init*b0eoverm*&
                                !(1*exp(-(qx(nx)-qx(i))**2/(2.0*dx)**2) + &
                                !exp(-(qx(1)-qx(i))**2/(2.0*dx)**2)) + nu_init
                                
                                
                                !if (i .le. 5) then
                        	!	if ((k .ge. nz/2.0-ddthickness) .and. (k .le. nz/2.0+ddthickness)) then
				!	nu(i,j,k) = 2*nu_init*b0eoverm*&
                                !(20*exp(-(qx(nx-2)-qx(nx/2))**2/(5.0*dx)**2) + &
                                !exp(-(qx(2)-qx(nx/2))**2/(5.0*dx)**2)) + nu_init
				!	!write(*,*) 'x...nu...n0...', i,nu(i,j,k),nu(nx/2,2,nz/2)
				!	endif
                              	!endif
                                
                                !if (i .ge. nx-5) then
                        	!	if ((k .ge. nz/2.0-ddthickness) .and. (k .le. nz/2.0+ddthickness)) then
				!	nu(i,j,k) = 2*nu_init*b0eoverm*&
                                !(20*exp(-(qx(nx-2)-qx(nx/2))**2/(5.0*dx)**2) + &
                                !exp(-(qx(2)-qx(nx/2))**2/(5.0*dx)**2)) + nu_init
				!	!write(*,*) 'x...nu...n0...', i,nu(i,j,k),nu(nx/2,2,nz/2)
				!	endif
                              	!endif
                                
                               ! if ((j .eq. 2) .and. (k .eq. nz/2) ) then
                               ! write(*,*) 'x...nu...', i,nu(i,j,k)
                               ! write(*,*) 'b,eoverm,beoverm', b0_init,eoverm,b0_init*eoverm
                               ! endif
                              nu(i,j,k) = nu_init
                          
                              do m = 1,3
                                    bt(i,j,k,m) = b0(i,j,k,m)
                                    b12(i,j,k,m) = b0(i,j,k,m)
                                    b1(i,j,k,m) = b0(i,j,k,m)
                                    b1p2(i,j,k,m) = b0(i,j,k,m)
                                    vol = dx_cell(i)*dy_cell(j)*dz_cell(k)*km_to_m**3
                                    input_Eb = input_Eb + (vol/(2.0*mu0))*(mO_q*b0(i,j,k,m))**2
                              enddo
                        enddo
                  enddo
            enddo
            
            
!            open(40,file='b0.dat',status='unknown',form='unformatted')
!            write(40) nz
!            write(40) b0
!            close(40)
            
      end subroutine grd6_setup
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 
      subroutine grd7()
            use grid
            use mult_proc, only: my_rank
            use inputs, only: dx,dy,delz,out_dir,zsf!,nrgrd
            implicit none
            integer, parameter:: nrgrd = 750
            integer:: i,j,k,ind
            real:: xsf,zplus,zminus,xplus,xminus,yplus,yminus
            
            rk = nz/2
            rj= ny/2
            ri = nx/2
!!!!!!!!!Unstretched grids!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
            do j=1,ny
                  qy(j) = j*dy
                  dy_grid(j) = dy
            enddo
            
            do i = 1,nx!-1
                  qx(i) = i*dx
                  dx_grid(i) = dx
            enddo
            !qx(nx) = dx*(nx+50) !10 is 10 additional distance for the ghost cell...
            !dx_grid(nx) = 10*dx

            do k = 1,nz
                qz(k) = k*delz
                dz_grid(k) = delz
            enddo

            
!!!!!!!Stretch X direction!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!            xsf = 0.0   !x stretch factor
!            !up from center
!            do i = ri, ri+nrgrd
!                  dx_grid(i) = dx
!            enddo
!            do i = ri + nrgrd + 1, nx
!                  dx_grid(i) = dx + xsf*dx*(i-(ri+nrgrd+1))/(nx-(ri+nrgrd+1))
!            enddo
!            !down from center
!            do i=ri-nrgrd, ri-1
!                  dx_grid(i) = dx
!            enddo
!            do i = 1, ri-nrgrd-1
!                  ind = ri-nrgrd-i
!                  dx_grid(ind) = dx + xsf*dx*(ri-nrgrd-1-ind)/(ri-nrgrd-1)
!            enddo
            
!            qx(1) = dx
            
!            do i = 2,nx
!                  qx(i) = qx(i-1)+dx_grid(i)
!            enddo
            
!            do i=1, nx-1
!                  dx_grid(i) = qx(i+1)-qx(i)
!            enddo
            
!            dx_grid(nx) = dx_grid(nx-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
            
            
!!!!!!!Stretch Z direction!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
!            zsf = 0.0   !z stretch factor
            !up from center
            
!            do k = rk, rk + nrgrd
!                  dz_grid(k) = delz
!            enddo
                      
            
!            do k = rk + nrgrd + 1, nz
!                  dz_grid(k) = delz + zsf*delz*(k-(rk+nrgrd+1))/(nz-(rk+nrgrd+1))
!            enddo
            
!            !down from center
!            do k = rk - nrgrd, rk - 1
!                  dz_grid(k) = delz
!            enddo
!            do k = 1, rk - nrgrd - 1
!                  ind = rk-nrgrd-k
!                  dz_grid(ind) =  delz + zsf*delz*(rk-nrgrd-1-ind)/(rk-nrgrd-1)
!            enddo
            
!            qz(1) = delz
            
!            do k=2,nz
!                  qz(k) = qz(k-1)+dz_grid(k)
 !           enddo
            
 !           do k=1, nz-1
 !                 dz_grid(k) = qz(k+1)-qz(k)
 !           enddo
            
 !           dz_grid(nz) = dz_grid(nz-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            dz_cell(1) = dz_grid(1)
            dz_cell(nz) = dz_grid(nz)
            zrat(1) = 0.5
            zrat(nz) = 0.5
            do k=2, nz-1
                  dz_cell(k) = ((qz(k+1) + qz(k))/2.0) - ((qz(k) + qz(k-1))/2.0)
                  zplus = (qz(k+1) + qz(k))/2.0
                  zminus = (qz(k) + qz(k-1))/2.0
                  zrat(k) = (qz(k) - zminus)/(zplus - zminus)
            enddo
            
            dx_cell(1) = dx_grid(1)
            dx_cell(nx) = dx_grid(nx)
            xrat(1) = 0.5
            xrat(nx) = 0.5
            do i=2, nx-1
                  dx_cell(i) = ((qx(i+1) + qx(i))/2.0) - ((qx(i) + qx(i-1))/2.0)
                  xplus = (qx(i+1) + qx(i))/2.0
                  xminus = (qx(i) + qx(i-1))/2.0
                  xrat(i) = (qx(i) - xminus)/(xplus - xminus)
            enddo
            
            dy_cell(1) = dy_grid(1)
            dy_cell(ny) = dy_grid(ny)
            yrat(1) = 0.5
            yrat(ny) = 0.5
            do j=2, ny-1
                  dy_cell(j) = ((qy(j+1) + qy(j))/2.0) - ((qy(j) + qy(j-1))/2.0)
                  yplus = (qy(j+1) + qy(j))/2.0
                  yminus = (qy(j) + qy(j-1))/2.0
                  yrat(j) = (qy(j) - yminus)/(yplus - yminus)
            enddo
            
            if (my_rank .eq. 0) then
                  open(40,file=trim(out_dir)//'c.coord.dat',status='unknown',form='unformatted')
                  
                  write(40) nx
                  write(40) ny
                  write(40) nz
                  write(40) qx
                  write(40) qy
                  write(40) qz
                  write(40) dz_grid
                  write(40) dz_cell
                  close(40)
                  
            endif
            
      end subroutine grd7
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      subroutine grid_gaussian()
! Generates a grid that scales the z axis so that the spacing between grid points is a multiple of lambda_i
! using an analytical gaussian.
            use grid
            use mult_proc, only: my_rank
            use inputs, only: dx,dy,delz,out_dir,zsf,nrgrd,grad,amp,q,mion,nf_init
            implicit none
            integer:: i,j,k,ind
            real:: xsf,zplus,zminus,xplus,xminus,yplus,yminus
            
            rk = nz/2
            rj= ny/2
            ri = nx/2
!!!!!!!!!Unstretched grids!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
            do j=1,ny
                  qy(j) = j*dy
                  dy_grid(j) = dy
            enddo
            
            do i = 1,nx
                  qx(i) = i*dx
                  dx_grid(i) = dx
            enddo
            
!!!!!!!!!Stretch z!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            qz(nz/2) = 0
            do k=rk+1, nz
                  qz(k) = qz(k-1) + zsf*3e8/1e3*sqrt(8.85e-12*mion/(q*q* &
                  (nf_init/1e9+nf_init/1e9*(amp-1.0)*exp(-(qz(k-1)/(grad*delz))**2))))
            enddo
            
!            do k = 1, rk-1
!                  ind = rk-k
!                  qz(ind) = qz(ind+1) - zsf*3e8/1e3*sqrt(8.85e-12*mion/(q*q* &
!                  (nf_init/1e9+nf_init/1e9*(amp-1.0)*exp(-(qz(ind+1)/(grad*delz))**2))))
!            enddo

            do k = 1 , rk-1
                  ind = rk - k
                  qz(ind) = -qz(rk+k)
            enddo
            zplus = qz(nz-1)
            do k=1,nz
                  qz(k) = qz(k) + zplus
            enddo
            
            do k=1, nz-1
                  dz_grid(k) = qz(k+1)-qz(k)
            enddo
            
            dz_grid(nz) = dz_grid(nz-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            dz_cell(1) = dz_grid(1)
            dz_cell(nz) = dz_grid(nz)
            zrat(1) = 0.5
            zrat(nz) = 0.5
            do k=2, nz-1
                  dz_cell(k) = ((qz(k+1) + qz(k))/2.0) - ((qz(k) + qz(k-1))/2.0)
                  zplus = (qz(k+1) + qz(k))/2.0
                  zminus = (qz(k) + qz(k-1))/2.0
                  zrat(k) = (qz(k) - zminus)/(zplus - zminus)
            enddo
            
            dx_cell(1) = dx_grid(1)
            dx_cell(nx) = dx_grid(nx)
            xrat(1) = 0.5
            xrat(nx) = 0.5
            do i=2, nx-1
                  dx_cell(i) = ((qx(i+1) + qx(i))/2.0) - ((qx(i) + qx(i-1))/2.0)
                  xplus = (qx(i+1) + qx(i))/2.0
                  xminus = (qx(i) + qx(i-1))/2.0
                  xrat(i) = (qx(i) - xminus)/(xplus - xminus)
            enddo
            
            dy_cell(1) = dy_grid(1)
            dy_cell(ny) = dy_grid(ny)
            yrat(1) = 0.5
            yrat(ny) = 0.5
            do j=2, ny-1
                  dy_cell(j) = ((qy(j+1) + qy(j))/2.0) - ((qy(j) + qy(j-1))/2.0)
                  yplus = (qy(j+1) + qy(j))/2.0
                  yminus = (qy(j) + qy(j-1))/2.0
                  yrat(j) = (qy(j) - yminus)/(yplus - yminus)
            enddo
            
            if (my_rank .eq. 0) then
                  open(40,file=trim(out_dir)//'c.coord.dat',status='unknown',form='unformatted')
                  
                  write(40) nx
                  write(40) ny
                  write(40) nz
                  write(40) qx
                  write(40) qy
                  write(40) qz
                  write(40) dz_grid
                  write(40) dz_cell
                  close(40)
                  
            endif

      end subroutine grid_gaussian

      
end module initial
            
