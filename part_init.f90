module part_init
      implicit none
      save
      contains
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine load_foreshock_Maxwellian(vth,Ni_tot_1,Ni_tot_2,mass,mratio,population)
            use dimensions
            use boundary
            use inputs, only: PI, vsw, dx, dy, km_to_m, beta_particle, kboltz, mion, nf_init,b0_init,mu0,boundx, q, mO, va_f, delz,ddthickness, plasma_beta,FSBeamWidth, FSThermalRatio,ddthickness,FSDriftSpeed,ForeshockBeta, TDpressureBalance
            use grid, only: qx,qy,qz,dz_grid
            use gutsp
            use var_arrays, only: np,vp,vp1,xp,up,Ni_tot,ijkp,m_arr,mrat,beta,beta_p,wght,temp_p,mix_ind,b0, b1,E,ExB,bt,up_cold
            implicit none
            integer(4), intent(in):: Ni_tot_1, Ni_tot_2
            real, intent(in):: mratio, mass, vth
            real:: Lo_y,eoverm
            integer:: population
            integer:: disp,montecarlo
            real:: vx, vy, vz, va, va_x, Temp, Tempcalc, pl_beta(nx,ny,nz), fitdist,randpop, EvBx, EvBy,EvBz,theta2,phi2, vxring, vyring, vzring,vring,vr,vtheta
            integer:: l,m,i,j,k

            disp = 0
            do i=1,nx
                  do j=1,ny
                        do k=1,nz
                              pl_beta(i,j,k) = plasma_beta
                        enddo
                  enddo
            enddo


            do l = Ni_tot_1,Ni_tot_2

                  ! X Component
                  if ((population .eq. 1) .or. (population .eq. 5) )then !Solar Wind
                        mix_ind(l) = 0
                        
                        if (boundx .eq. 5) then
                              xp(l,1) = 0.0*qx(1)+(1.0-pad_ranf())*(qx(2)-qz(1)) +(1.0)*(qx(2)-qz(1)) !Shock real, to propagate SW/field
                        endif
                        if (boundx .eq. 4) then
                              xp(l,1) = 1.0*qx(1)+(1.0-pad_ranf())*(qx(2)-qz(1)) !Normal runs, open boundary
                        endif
                                                
                  else if (population .eq. 2 .or. population .eq. 6) then !Foreshock Right
                  
                        mix_ind(l) = 1
                        xp(l,1) = qx(nx-1)+(1.0-pad_ranf())*(qx(2)-qx(1))
                        ! xp(l,1) = (1.0-pad_ranf())*(qx(1))!+qx(1)
                        !xp(l,1) = (1.0-pad_ranf())*(qx(1))!+qx(1))
                        
                  else if ((population .eq. 0) .or. (population .eq. 4)) then !Solar Wind, Everywhere
                  
                        mix_ind(l) = 0
                        xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)+qx(1)) !nx-1 for periodic, nx for nonperiodic in x
                        !xp(l,1) = (1.0-pad_ranf())*(qx(nx)-qx(1))
                        
                  else if (population .eq. 3) then !Foreshock, initial Beam
                  
                        mix_ind(l) = 1
                        !xp(l,1) = qx(nx/2 - nx/5)+(1.0-pad_ranf())*(qx(nx)-qx(nx/2-nx/5))
                        !xp(l,1) = qx(nx/2)+(1.0-pad_ranf())*(qx(nx)-qx(nx/2))
                        !xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(2)-qz(1))
                        xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)+qx(1)) !nx-1 for periodic, nx for nonperiodic in x
                        
                  else if ((population .eq. 7) .or. (population .eq. 7) )then !Foreshock Left
                  
                        mix_ind(l) = 1
                        xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(2)-qz(1))
                        
                  endif




                  ! Y Component
                  xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1) )




                  ! Z Component
                  if ((population .eq. 2) .or. (population .eq. 7) .or. (population .eq. 3) ) then !Lower Beam
                        !Foreshock ions within TD have slower Vx.
                        if ((population .eq. 2) .or. (population .eq. 7)) then
                        montecarlo = 0
                        
                              do 22 while (montecarlo .eq. 0)
                              
                                    xp(l,3) = qz((nz)/2) - qz(1)*FSBeamWidth +(1.0-pad_ranf())*(qz(1)*FSBeamWidth)

                                    call get_pindex(i,j,k,l)

                                    fitdist = (   FSDriftSpeed*va_f*( b0(i,j,k+int(TDpressureBalance),1) / (  sqrt( b0(i,j,k+int(TDpressureBalance),1)**2+b0(i,j,k+int(TDpressureBalance),2)**2+b0(i,j,k+int(TDpressureBalance),3)**2 ) )  ) - ( (1.0-0*TDpressureBalance)*(va_f)*b0(i,j,k+int(TDpressureBalance),2)*b0(i,j,k+int(TDpressureBalance),2) / ( ( b0(i,j,k+int(TDpressureBalance),1)**2+b0(i,j,k+int(TDpressureBalance),2)**2+b0(i,j,k+int(TDpressureBalance),3)**2 ) ) ) /(FSDriftSpeed*va_f))
                                    fitdist = (   ( b0(i,j,k+8*int(TDpressureBalance),1) / (  sqrt( b0(i,j,k+8*int(TDpressureBalance),1)**2+b0(i,j,k+8*int(TDpressureBalance),2)**2+b0(i,j,k+8*int(TDpressureBalance),3)**2 ) )  ) - ( (1.0-0*TDpressureBalance)*b0(i,j,k+8*int(TDpressureBalance),2)*b0(i,j,k+8*int(TDpressureBalance),2) / ( ( b0(i,j,k+8*int(TDpressureBalance),1)**2+b0(i,j,k+8*int(TDpressureBalance),2)**2+b0(i,j,k+8*int(TDpressureBalance),3)**2 ) ) ))
                                    fitdist= b0(i,j,k+0*int(TDpressureBalance),1)/sqrt( b0(i,j,k+0*int(TDpressureBalance),1)**2+b0(i,j,k+0*int(TDpressureBalance),2)**2+b0(i,j,k+0*int(TDpressureBalance),3)**2 )
                                    
                                    if (pad_ranf() .le. fitdist) montecarlo = 1

                              22 continue

                        endif

                  else if ((population .eq. 4) .or. (population .eq. 5))then !TD
                        !Fill in a bit more in solar wind within TD
                        montecarlo = 0
                        do 20 while (montecarlo .eq. 0) 

                              xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz-1)-qz(1))
                              
                              if ( xp(l,3) .le. qz(nz/2) ) then
                                    fitdist = 1 - 0.5 - 0.5*tanh( (qz(nz/2) - xp(l,3)) / (ddthickness*delz) )
                              else if ( xp(l,3) .gt. qz(nz/2) ) then
                                    fitdist = 1 - 0.5 - 0.5*tanh( (xp(l,3)- qz(nz/2)) / (ddthickness*delz) )
                              endif

                              if (pad_ranf() .le. fitdist) montecarlo = 1

                        20 continue
                        
                  else if (population .eq. 6) then !Upper Beam
                        montecarlo = 0
                        do 23 while (montecarlo .eq. 0)
                              xp(l,3) = qz((nz)/2) + qz(1)*FSBeamWidth -(1.0-pad_ranf())*(qz(1)*FSBeamWidth)
                              call get_pindex(i,j,k,l)
                              
                              fitdist= b0(i,j,k+0*int(TDpressureBalance),2)/sqrt( b0(i,j,k+0*int(TDpressureBalance),1)**2+b0(i,j,k+0*int(TDpressureBalance),2)**2+b0(i,j,k+0*int(TDpressureBalance),3)**2 )!
                             
                              if (pad_ranf() .le. fitdist) montecarlo = 1     

                        23 continue
      
                  else !All other Z
                        xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz-1)-qz(1))
                  endif

                  m_arr(l) = mass
                  mrat(l) = mratio

                  va = b0_init/sqrt(mu0*mion*nf_init/1e9)/1e3
                  eoverm = q/mO

                  vx = vth*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
                  vy = vth*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
                  vz = vth*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())

                  
                  
                  if ((population .eq. 1) .or. (population .eq. 5) .or. (population .eq. 0) .or. (population .eq. 4)) then !Solar Wind, Left Edge
                        vp(l,1) = va_f*va+vx
                        vp(l,2) = vy
                        vp(l,3) = vz !- 1.0*va
                        
                        beta_p(l) = beta_particle

                  else if ((population .eq. 2) .or. (population .eq. 3) .or. (population .eq. 7) .or. (population .eq. 6)) then !Foreshock Ions, Right Edge
                        call get_pindex(i,j,k,l)

                        i=i-1
                        EvBx = -( ( up(i,j,k,2)*b1(i,j,k,3) - up(i,j,k,3)*b1(i,j,k,2) )  )
                        EvBy =  ( ( up(i,j,k,1)*b1(i,j,k,3) - up(i,j,k,3)*b1(i,j,k,1) )  )
                        EvBz = -( ( up(i,j,k,1)*b1(i,j,k,2) - up(i,j,k,2)*b1(i,j,k,1) )  )
                                  
                              
                        !Constant solar wind bulk flow for ExB
                        !*vp(l,1) = -FSDriftSpeed*va_f*va*( b0(i,j,k,1) / (  sqrt( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) ) ) + FSThermalRatio*vx + ( (va_f*va)*b0(i,j,k,2) )*b0(i,j,k,2) / ( ( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) )
                        !*vp(l,2) = -FSDriftSpeed*va_f*va*( b0(i,j,k,2) / (  sqrt( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) ) ) + FSThermalRatio*vy - ( (va_f*va)*b0(i,j,k,2) )*b0(i,j,k,1) / ( ( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) )
                        !*vp(l,3) = FSThermalRatio*vz    
                        
                              
                        !Current B (local) for ExB using current (local) ion flow. (2022)
                        vp(l,1) = -FSDriftSpeed*va_f*va*( b1(i,j,k,1) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + FSThermalRatio*vx + ( EvBy*b1(i,j,k,3) - EvBz*b1(i,j,k,2) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
                        vp(l,2) = -FSDriftSpeed*va_f*va*( b1(i,j,k,2) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + FSThermalRatio*vy - ( EvBx*b1(i,j,k,3) - EvBz*b1(i,j,k,1) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
                        vp(l,3) = -FSDriftSpeed*va_f*va*( b1(i,j,k,3) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + FSThermalRatio*vz + ( EvBx*b1(i,j,k,2) - EvBy*b1(i,j,k,1) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
                        
                        !Temperature Anisotrophy
                        ! vp(l,1) = vp(l,1) - 4*vx*  ( b1(i,j,k,1) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) 
                        ! vp(l,2) = vp(l,2) - 4*vx*  ( b1(i,j,k,2) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) 
                        ! vp(l,3) = vp(l,3) - 4*vx*  ( b1(i,j,k,3) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) 
                  
                  
                        !Gyrophase Bunched with local B !Updated 12/22/2021 sqrt(8.0*8.0-4.0*4.0), using an additional Vperp
                        !vp(l,1) = -FSDriftSpeed*va_f*va*( b1(i,j,k,1) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + FSThermalRatio*vx                                 + ( EvBy*b1(i,j,k,3) - EvBz*b1(i,j,k,2) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
                        !vp(l,2) = -FSDriftSpeed*va_f*va*( b1(i,j,k,2) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + FSThermalRatio*vy + 0*sqrt(8.0*8.0-4.0*4.0)*vth   - ( EvBx*b1(i,j,k,3) - EvBz*b1(i,j,k,1) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
                        !vp(l,3) = -FSDriftSpeed*va_f*va*( b1(i,j,k,3) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + FSThermalRatio*vz + 6.0*vth	                  + ( EvBx*b1(i,j,k,2) - EvBy*b1(i,j,k,1) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
                        
                        !Gyrophase Bunched with local B
                        !vp(l,1) = -FSDriftSpeed*va_f*va*( b1(i,j,k,1) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + FSThermalRatio*vx 				           + ( EvBy*b1(i,j,k,3) - EvBz*b1(i,j,k,2) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
                        !vp(l,2) = -FSDriftSpeed*va_f*va*( b1(i,j,k,2) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + FSThermalRatio*vy - sqrt(2.0)/2.0*6*vth           - ( EvBx*b1(i,j,k,3) - EvBz*b1(i,j,k,1) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
                        !vp(l,3) = -FSDriftSpeed*va_f*va*( b1(i,j,k,3) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + FSThermalRatio*vz + sqrt(2.0)/2.0*6*vth           + ( EvBx*b1(i,j,k,2) - EvBy*b1(i,j,k,1) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
                        

                        !Gyrophase Bunched with constant B
                        !EvBx = -( ( 0.0*b0(i,j,k,3) - 0.0*b0(i,j,k,2) )  )
                        !EvBy =  ( ( va_f*va*b0(i,j,k,3) - 0.0*b0(i,j,k,1) )  )
                        !EvBz = -( ( va_f*va*b0(i,j,k,2) - 0.0*b0(i,j,k,1) )  )
                        
                        !vp(l,1) = -FSDriftSpeed*va_f*va*( b0(i,j,k,1) / (  sqrt( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) ) ) + FSThermalRatio*vx                                    + ( EvBy*b0(i,j,k,3) - EvBz*b0(i,j,k,2) ) / ( ( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) )
                        !vp(l,2) = -FSDriftSpeed*va_f*va*( b0(i,j,k,2) / (  sqrt( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) ) ) + FSThermalRatio*vy + sqrt(2.0)/2.0*FSThermalRatio*vth - ( EvBx*b0(i,j,k,3) - EvBz*b0(i,j,k,1) ) / ( ( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) )
                        !vp(l,3) = -FSDriftSpeed*va_f*va*( b0(i,j,k,3) / (  sqrt( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) ) ) + FSThermalRatio*vz + sqrt(2.0)/2.0*FSThermalRatio*vth + ( EvBx*b0(i,j,k,2) - EvBy*b0(i,j,k,1) ) / ( ( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) )
                        
                        
                        !Shell distribution beam +exB
                        !theta2 = pad_ranf()*2*PI
                        !phi2   = pad_ranf()*PI   
                        !vp(l,1) = -FSDriftSpeed*va_f*va*( b1(i,j,k,1) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + FSThermalRatio*vth*cos(theta2)*sin(phi2) + ( EvBy*b1(i,j,k,3) - EvBz*b1(i,j,k,2) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) 
                        !vp(l,2) = -FSDriftSpeed*va_f*va*( b1(i,j,k,2) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + FSThermalRatio*vth*sin(theta2)*sin(phi2) - ( EvBx*b1(i,j,k,3) - EvBz*b1(i,j,k,1) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
                        !vp(l,3) = -FSDriftSpeed*va_f*va*( b1(i,j,k,3) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + FSThermalRatio*vth*cos(phi2) + ( EvBx*b1(i,j,k,2) - EvBy*b1(i,j,k,1) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
                        
                        
                        
                        !!!!!!!!!  
                        
                        !Ring beam +exB Updated 3/19/2022 for vperp parameter scan, constant B
                  !     i=i-1
                  !   	EvBx = -( ( up(i,j,k,2)*b0(i,j,k,3) - up(i,j,k,3)*b0(i,j,k,2) )  )
                  !   	EvBy =  ( ( up(i,j,k,1)*b0(i,j,k,3) - up(i,j,k,3)*b0(i,j,k,1) )  )
                  !   	EvBz = -( ( up(i,j,k,1)*b0(i,j,k,2) - up(i,j,k,2)*b0(i,j,k,1) )  )
                        
                  !	vring = 10.0
                  !	vr = sqrt(vy*vy + vz*vz)
                  !	vtheta = atan(vz/vy)
                        
                  !	vyring = (FSThermalRatio*vr+vring*vth)*cos(vtheta)
                  !	vzring = (FSThermalRatio*vr+vring*vth)*sin(vtheta)
                  !	vy = sign(vyring,vy)
                  !	vz = vzring
                              
                  !	vp(l,1) = -FSDriftSpeed*va_f*va*( b1(i,j,k,1) / (  sqrt( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) ) ) + 0*FSThermalRatio*vx + ( EvBy*b0(i,j,k,3) - EvBz*b0(i,j,k,2) ) / ( ( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) )
                  !     vp(l,2) = -FSDriftSpeed*va_f*va*( b1(i,j,k,2) / (  sqrt( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) ) ) + vy                - ( EvBx*b0(i,j,k,3) - EvBz*b0(i,j,k,1) ) / ( ( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) )
                  !     vp(l,3) = -FSDriftSpeed*va_f*va*( b1(i,j,k,3) / (  sqrt( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) ) ) + vz                + ( EvBx*b0(i,j,k,2) - EvBy*b0(i,j,k,1) ) / ( ( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) )
                        
                  !!!!!!!!!1
                        
                        !Ring beam +exB Updated 12/22/2021
                  !	vring = sqrt(8.0*8.0-4.0*4.0)
                  !	vr = sqrt(vy*vy + vz*vz)
                  !	vtheta = atan(vz/vy)
                        
                  !	vyring = (FSThermalRatio*vr+vring*vth)*cos(vtheta)
                  !	vzring = (FSThermalRatio*vr+vring*vth)*sin(vtheta)
                  !	vy = sign(vyring,vy)
                  !	vz = vzring
                              
                  !	vp(l,1) = -FSDriftSpeed*va_f*va*( b1(i,j,k,1) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + FSThermalRatio*vx + ( EvBy*b1(i,j,k,3) - EvBz*b1(i,j,k,2) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
                  !     vp(l,2) = -FSDriftSpeed*va_f*va*( b1(i,j,k,2) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + vy - ( EvBx*b1(i,j,k,3) - EvBz*b1(i,j,k,1) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
                  !     vp(l,3) = -FSDriftSpeed*va_f*va*( b1(i,j,k,3) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + vz + ( EvBx*b1(i,j,k,2) - EvBy*b1(i,j,k,1) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
                        
                        
                        !Ring Beam, Vth=0, Vperp = 8
                  !     vring = 8.0
                  !	vr = sqrt(vy*vy + vz*vz)
                  !	vtheta = atan(vz/vy)
                        
                  !	vyring = (vring*vth)*cos(vtheta)
                  !	vzring = (vring*vth)*sin(vtheta)
                  !	vy = sign(vyring,vy)
                  !	vz = vzring
                        
                  !     vp(l,1) = -FSDriftSpeed*va_f*va*( b1(i,j,k,1) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + 0  + ( EvBy*b1(i,j,k,3) - EvBz*b1(i,j,k,2) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
                  !     vp(l,2) = -FSDriftSpeed*va_f*va*( b1(i,j,k,2) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + vy - ( EvBx*b1(i,j,k,3) - EvBz*b1(i,j,k,1) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
                  !     vp(l,3) = -FSDriftSpeed*va_f*va*( b1(i,j,k,3) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + vz + ( EvBx*b1(i,j,k,2) - EvBy*b1(i,j,k,1) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
                        
                        
                        
                        
                        !Ring beam
                  !        vxring = 0*0.5*FSThermalRatio*vth *( b1(i,j,k,1) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) )
                  !        vyring = 0.5*1.0*FSThermalRatio*vth !1.5*1.0
                  !        vzring = 0.5*1.0*FSThermalRatio*vth !1.5*1.0	


                  !        montecarlo = 0
                  !        do 24 while (montecarlo .eq. 0)
                  !        	vx = vth*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
                  !        	if (FSThermalRatio*vx .ge. vxring) montecarlo = 1 !cut-off is half of thermal.
                  !        24 continue
                        
                  !        montecarlo = 0
                  !        do 25 while (montecarlo .eq. 0)
                  !        	vy = vth*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
                  !        	if (abs(FSThermalRatio*vy) .ge. abs(vyring)) montecarlo = 1 !cut-off is half of thermal.
                  !        25 continue
                        
                  !        montecarlo = 0
                  !        do 26 while (montecarlo .eq. 0)
                  !        	vz = vth*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
                  !        	if (abs(FSThermalRatio*vz) .ge. abs(vzring)) montecarlo = 1 !cut-off is half of thermal.
                  !        26 continue
                        
                  !        vp(l,1) = -FSDriftSpeed*va_f*va*( b1(i,j,k,1) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + FSThermalRatio*vx + ( EvBy*b1(i,j,k,3) - EvBz*b1(i,j,k,2) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
                  !        vp(l,2) = -FSDriftSpeed*va_f*va*( b1(i,j,k,2) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + FSThermalRatio*vy - ( EvBx*b1(i,j,k,3) - EvBz*b1(i,j,k,1) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
                  !        vp(l,3) = -FSDriftSpeed*va_f*va*( b1(i,j,k,3) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + FSThermalRatio*vz + ( EvBx*b1(i,j,k,2) - EvBy*b1(i,j,k,1) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )

                        beta_p(l) = (ForeshockBeta)*beta_particle !2 for generating half of z side and 10 for 10% of FS ions.
                  endif

                  do m=1,3
                        vp1(l,m) = vp(l,m)
                  enddo

            enddo

            call get_interp_weights()

            call update_np()

            call update_up(vp)

      end subroutine load_foreshock_Maxwellian
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module part_init
