module part_init
      implicit none
      save
      contains

      subroutine Energy_diag(Evp,Euf,EB1,EB1x,EB1y,EB1z,EE,EeP)
            use dimensions
            use mpi
            use mult_proc, only: my_rank
            use grid, only: dx_cell,dy_cell,dz_cell
            use inputs, only: mion,q,mu0,mO,km_to_m,epsilon
            use var_arrays, only: xp,vp,b0,b1,E,nu,up,np,Ni_tot,beta,beta_p,input_E,prev_Etot,bndry_Eflux,m_arr,ionPos
            implicit none
            real, intent(out):: Euf,EB1,EB1x,EB1y,EB1z,EE,EeP,Evp
            real:: denf,m_q,recvbuf,total_E,aveEvp,norm_E,vol
            real:: S_Evp,S_input_E
            integer:: count, ierr
            integer:: i,j,k,m,l
            
            count = 1
            m_q = mion/q
            
            Euf = 0.0
            EB1 = 0.0
            EB1x = 0.0
            EB1y = 0.0
            EB1z = 0.0
            EE = 0.0
            EeP = 0.0
            
            do i=1,nx-1
!                  do j = 1,ny-1
                   j=2
                        do k=1,nz-1
                              vol= dx_cell(i)*dy_cell(j)*dz_cell(k)*km_to_m**3
                              EB1x=EB1x + (vol/(2.0*mu0))*(m_q*b1(i,j,k,1))**2
                              EB1y=EB1y + (vol/(2.0*mu0))*(m_q*b1(i,j,k,2))**2
                              EB1z=EB1z + (vol/(2.0*mu0))*(m_q*b1(i,j,k,3))**2
                                    do m=1,3
                                          denf = np(i,j,k)/(km_to_m**3)
                                          Euf = Euf + 0.5*mO*denf*vol*(up(i,j,k,m)*km_to_m)**2
                                          EB1 = EB1 + (vol/(2.0*mu0))*(m_q*(b1(i,j,k,m)-b0(i,j,k,m)))**2
                                          EE = EE + (epsilon*vol/2.0)*(m_q*E(i,j,k,m)*km_to_m)**2
                                    enddo
                        enddo
!                  enddo
            enddo
            
            Evp = 0.0
            do l=1, Ni_tot
                  do m=1,3
                        Evp = Evp + 0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 / (beta*beta_p(l))
                  enddo
            enddo
            
            
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            
            call MPI_ALLREDUCE(Evp,recvbuf,count,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
            S_Evp = recvbuf
            
            call MPI_ALLREDUCE(input_E,recvbuf,count,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
            S_input_E = recvbuf
            
            total_E = S_Evp + EE + EB1
            aveEvp = S_Evp/S_input_E
            
            
            if (my_rank .eq. 0) then
                  !write(*,*) 'Normalized energy.................',total_E/S_input_E
                  !write(*,*) 'Normalized energy (bndry).........',total_E/(S_input_E+bndry_Eflux)
                 ! write(*,*) '          particles...............',S_Evp/S_input_E       
                  !write(*,*) '          field...................',(EE+EB1)/S_input_E
                  write(501) m
                  write(501) S_Evp/S_input_E
                  write(501) (EE+EB1)/S_input_E
            endif
            
            norm_E = total_E/S_input_E
            prev_Etot = norm_E

      end subroutine Energy_diag
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine load_Maxwellian(pbeta,Ni_tot_1,mass,mratio)
            use dimensions
            use boundary
            use inputs, only: PI, vsw, dx, dy, km_to_m, beta_particle, kboltz, mion, amp, grad, nf_init,b0_init,mu0
            use grid, only: qx,qy,qz,dz_grid
            use gutsp
            use var_arrays, only: np,vp,vp1,xp,input_p,up,Ni_tot,input_E,ijkp,m_arr,mrat,beta,beta_p,wght,grav,temp_p
            implicit none
            integer(4), intent(in):: Ni_tot_1 
            real, intent(in):: mratio, mass, pbeta
                                  
            integer:: disp
            real:: vx, vy, vz, Temp, Tempcalc, vth
            integer:: l,m,i,j,k
            
            disp = 0 !Displacement of gradient
!            amp = 100.0
!            grad = 100.0 ! density gradient (larger = more gradual
            

!            va = b0_init/sqrt(mu0*mion*nf_init/1e9)/1e3
      !      vth = sqrt(pbeta*va**2)
            
            do l = Ni_tot_1,Ni_tot
                  xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
                  xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
                  xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz-1)-qz(1))
                  m_arr(l) = mass
                  mrat(l) = mratio

!                  beta_p(l) = 1.0/(beta_particle+beta_particle*amp*exp(-((xp(l,3)-qz(nz/2-disp))/ &
!                        (grad*dz_grid(nz/2-disp)))**2))
                  beta_p(l) = beta_particle

                  call get_pindex(i,j,k,l)
!                  vth2=sqrt(vth*vth*beta_p(l)) !thermal speed dependent on np to set up pressure balance for density gradient

!                  vth2=va*sqrt(pl_beta(ijkp(l,1),ijkp(l,2),ijkp(l,3)))
                 

                  
                  vx = 40.0*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf()) !remember to add in vsw to get the flow velocity
                  vy = 40.0*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
                  vz = 40.0*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
                  
!                  ii = ijkp(l,1)
!                  kk = ijkp(l,3)
                  
!                  vp(l,1) = -0.0*(exp(-(xp(l,3)-qz(nz/2))**2/(10.*delz)**2)
!               x        *exp(-(xp(l,1)-qx(nx/2))**2/(10.*dx)**2))+vx
                  vp(l,1) = vx+57.0*exp(-(xp(l,3)-qz(nz/2))**2/(5*dz_grid(nz/2))**2) !Gaussian velocity perturbation (20)
                  vp(l,2) = vy 
                  vp(l,3) = vz 
                  
                  do m=1,3
                        vp1(l,m) = vp(l,m)
                        input_E = input_E + 0.5*m_arr(l)*(vp(l,m)*km_to_m)**2/(beta * beta_p(l))
                        input_p(m) = input_p(m) + m_arr(l) * vp(l,m) / (beta * beta_p(l))
                  enddo
                  
            enddo
            
            call get_interp_weights()
            call update_np()
            call update_up(vp)

            
            ! Add a centrifugal gravity term to keep the plasma confined to the torus.  Use T * dn/dz = nmg.  
            ! Depends on the density gradient.  Currently set as a gaussian.
            
!            Temp = vth**2/(3*kboltz)*mion*1.48*10-23!8.61738e-5
!            write(*,*) 'vth.................', vth
!            write(*,*) 'boltzman............', kboltz
!            write(*,*) 'temperature(analytic)..', Temp
!            call get_temperature()
!            Tempcalc = sum(temp_p(2,2,1:(nz-1)))/1e6/(nz-1) !in kg km^2/s^2
!            write(*,*) 'temperature (2,2,100)..', temp_p(2,2,2:10)/1.6e-19
!            stop
            
            do i=1,nx
            do j=1,ny
            do k=1,nz
                  ! Gravity is based on the analytical expression for the density profile (look at beta_p)
                  ! np = const/(beta*beta_p), and grav = const * (dn/dx) / n
                  

                        
                       
                        
!                  grav(i,j,k) = -2.0*Tempcalc/(mion*(grad*dz_grid(nz/2-disp))**2 &
!                        *(1.0+amp*exp(-((qz(k)-qz(nz/2-disp))/(grad*dz_grid(nz/2-disp)))**2))) &
!                        *amp*(qz(k)-qz(nz/2-disp))*exp(-((qz(k)-qz(nz/2-disp))/(grad*dz_grid(nz/2-disp)))**2)   
                 
       
                 grav(i,j,k) = 0.0
                  
                  
            enddo
            enddo
            enddo
           ! write(*,*) 'gravity...', grav(2,2,nz/2+50), grav(2,2,nz/2-50)
           ! stop

           call count_ppc() 
      end subroutine load_Maxwellian
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine load_ring_beam(vring,dNi,mass,mratio)
            use dimensions
            use boundary
            use inputs, only: PI, vsw, dx, dy, km_to_m, beta_pu, ion_amu,m_pu,beta_particle, amp, grad
            use grid, only: qx,qy,qz,dz_grid
            use gutsp
            use var_arrays, only: np,vp,vp1,xp,input_p,up,Ni_tot,input_E,ijkp,m_arr,mrat,beta,beta_p,wght
            implicit none
            integer(4), intent(in):: dNi 
            real, intent(in):: vring, mass,mratio
                                  
            integer:: disp, flg, l1
            real:: vth2, vx, vy, vz, rand1, theta2
            integer:: i,j,k,l,m
            
            disp = 0 !Displacement of gradient
!            amp = 100.0
!            grad = 400.0 ! density gradient (larger = more gradual
            
!            v1=1.0
            l1=Ni_tot+1
            
            do l = l1, l1+dni-1
                 ! xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
                  xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
                  flg=0
                        do while (flg .eq. 0)
                            !  xp(l,3) = qz(nz/2-20) + (1.0-pad_ranf())*(qz(nz/2+20)-qz(nz/2-20))
                              xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz-1)-qz(1))
                              rand1=pad_ranf()
                              if (exp(-(xp(l,3)-qz(nz/2))**2/(10*dz_grid(nz/2)**2)) .gt. rand1) then
                                    flg = 1
                              endif
                        enddo
                        
                  flg=0
                        do while (flg .eq. 0)
                              xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
                              rand1=pad_ranf()
                              if (exp(-(xp(l,1)-qx(nx/2))**2/(10*dx**2)) .gt. rand1) then
                                    flg = 1
                              endif
                        enddo
                        
                        beta_p(l) = beta_particle/10.0
                        m_arr(l) = mass
                        mrat(l) = mratio
                  
                  
                  call get_pindex(i,j,k,l)
                  
                  
!                  ii = ijkp(l,1)
!                  kk = ijkp(l,3)
                  
!                  vp(l,1) = -0.0*(exp(-(xp(l,3)-qz(nz/2))**2/(10.*delz)**2)
!               x        *exp(-(xp(l,1)-qx(nx/2))**2/(10.*dx)**2))+vx

            !   Ring beam velocity initializtion
!                  theta2 = pad_ranf()*2*PI
!                  vp(l,1) = vring*cos(theta2)
!                  vp(l,2) = vring*sin(theta2)
!                  vp(l,3) = 0.0
                  
            !   Maxwellian thermal distribution 
            
                  vth2=100.0;
                  
                  vx = vth2*sqrt(-log(pad_ranf()))*cos(2*PI*pad_ranf()) !remember to add in vsw to get the flow velocity
                  vy = vth2*sqrt(-log(pad_ranf()))*cos(2*PI*pad_ranf())
                  vz = vth2*sqrt(-log(pad_ranf()))*cos(2*PI*pad_ranf())
                  
                  
                  do m=1,3
                        vp1(l,m) = vp(l,m)
                        input_E = input_E + 0.5*m_arr(l)*(vp(l,m)*km_to_m)**2/(beta * beta_p(l))
                        input_p(m) = input_p(m) + m_arr(l) * vp(l,m) / (beta * beta_p(l))
                  enddo
                  
            enddo
            Ni_tot=Ni_tot+dNi
            
            call get_interp_weights()
            call update_np()
            call update_up(vp)
            
      end subroutine load_ring_beam
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine load_const_ppc(pbeta,begin,mass,mratio)
            use mult_proc, only: my_rank, procnum
            use dimensions
            use boundary
            use inputs, only: PI, vsw, dx, dy, km_to_m, beta_particle, kboltz, mion, amp, grad, nf_init,b0_init,mu0,ppc
            use grid, only: qx,qy,qz,dz_grid,dy_grid,dx_grid
            use gutsp
            use var_arrays, only: np,vp,vp1,xp,input_p,up,Ni_tot,Ni_tot_sys,Ni_init,input_E,ijkp,m_arr,mrat,beta,beta_p,wght,grav,temp_p
            use misc
            implicit none
            integer(4), intent(in):: begin 
            real, intent(in):: mratio, mass, pbeta
                                  
            integer:: disp,ppcpp,extra_ppc
            integer(4):: Ni_tot_1
            real:: vx, vy, vz, vth, Temp, Tempcalc,vol,new
            integer:: l,m,i,j,k
            
            disp = 0 !Displacement of gradient
            
            ppcpp=int(ppc/procnum)
        !    vth=sqrt(pbeta*va**2)
            
            Ni_tot_1 = begin
            do i=1,nx-2
            do j=1,ny-2
            do k=1,nz-2
            vol = dx_grid(i)*dy_grid(j)*dz_grid(k)  !km^3
            new = 10.0*exp(-((qz(k)-qz(nz/2-disp))/(grad*dz_grid(nz/2-disp)))**2)
            if (new .lt. 1.0) then
                  if (new .gt. pad_ranf()) then
                        new = 1.0
                  else 
                        new = 0.0
                  endif
            endif
            extra_ppc = nint(new)
            
            do l = Ni_tot_1, Ni_tot_1+ppcpp+extra_ppc-1
                  xp(l,1) = qx(i)+(1.0-pad_ranf())*(qx(i+1)-qx(i))
                  xp(l,2) = qy(j)+(1.0-pad_ranf())*(qy(j+1)-qy(j))
                  xp(l,3) = qz(k)+(1.0-pad_ranf())*(qz(k+1)-qz(k))
                  m_arr(l) = mass
                  mrat(l) = mratio

!                  beta_p(l) = 1.0/(beta_particle+beta_particle*amp*exp(-((xp(l,3)-qz(nz/2-disp))/ &
!                        (grad*dz_grid(nz/2-disp)))**2))

!                  beta_p(l) = beta_particle
                  beta_p(l) = (ppcpp+extra_ppc)*procnum/(nf_init+nf_init*(amp-1.0)*exp(-((qz(k)-qz(nz/2-disp))/(grad*dz_grid(nz/2-disp)))**2))/vol
!                   beta_p(l) = ppcpp*procnum/nf_init/vol

!                  call get_pindex(i,j,k,l)
!                  vth2=sqrt(vth*vth*beta_p(l)) !thermal speed dependent on np to set up pressure balance for density gradient

!                  vth2=va*sqrt(pl_beta(ijkp(l,1),ijkp(l,2),ijkp(l,3)))


                  
                  vx = 40.0*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf()) !remember to add in vsw to get the flow velocity
                  vy = 40.0*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
                  vz = 40.0*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
                  
!                  ii = ijkp(l,1)
!                  kk = ijkp(l,3)
                  
!                  vp(l,1) = -0.0*(exp(-(xp(l,3)-qz(nz/2))**2/(10.*delz)**2)
!               x        *exp(-(xp(l,1)-qx(nx/2))**2/(10.*dx)**2))+vx
                  vp(l,1) = vx+57.0!*exp(-(xp(l,3)-qz(nz/2))**2/(30*dz_grid(nz/2))**2) !Gaussian velocity perturbation (20)
                  vp(l,2) = vy 
                  vp(l,3) = vz 
                  
                  
                  
            enddo
            
            Ni_tot_1 = Ni_tot_1+ppcpp+extra_ppc
            Ni_tot=Ni_tot_1-1
            
            enddo
            enddo
            enddo
            
            Ni_tot_sys = Ni_tot*procnum
            beta=1.0;
            
            do l = begin,Ni_tot
                  call get_pindex(i,j,k,l)
                  do m=1,3
                        vp1(l,m) = vp(l,m)
                        input_E = input_E + 0.5*m_arr(l)*(vp(l,m)*km_to_m)**2/(beta * beta_p(l))
                        input_p(m) = input_p(m) + m_arr(l) * vp(l,m) / (beta * beta_p(l))
                  enddo
            enddo
          
            if (my_rank .eq. 0) then
                  write(*,*) 'Particles per processor.....', Ni_tot
                  write(*,*) 'Particles per cell..........', Ni_tot_sys/((nx-2)*(ny-2)*(nz-2)), ppcpp*procnum
                  write(*,*) 'Total particles.............', Ni_tot_sys
                  write(*,*) '***************************************************'
            endif
            
            call get_interp_weights()
            call update_np()
            call update_up(vp)

            
            ! Add a centrifugal gravity term to keep the plasma confined to the torus.  Use T * dn/dz = nmg.  
            ! Depends on the density gradient.  Currently set as a gaussian.
            
!            Temp = vth**2/(3*kboltz)*mion*1.48*10-23!8.61738e-5
!            write(*,*) 'vth.................', vth
!            write(*,*) 'boltzman............', kboltz
!            write(*,*) 'temperature(analytic)..', Temp
            call get_temperature()
            Tempcalc = sum(temp_p(2,2,1:(nz-1)))/1e6/(nz-1) !in kg km^2/s^2
!            write(*,*) 'temperature (2,2,100)..', temp_p(2,2,2:10)/1.6e-19
!            stop
            
            do i=1,nx
            do j=1,ny
            do k=1,nz
                  ! Gravity is based on the analytical expression for the density profile (look at beta_p)
                  ! np = const/(beta*beta_p), and grav = const * (dn/dx) / n.  Units are in km/s^2.
                  
            !      grav = (T*dn/dz)/(n*m)  where n = No*(1+(amp-1)*exp(-(z-z0)^2/dz^2))
                     
            !     grav(i,j,k) = 0.0
                 grav(i,j,k) = -2.0*Tempcalc*(amp-1.0)*(qz(k)-qz(nz/2-disp))/(grad*dz_grid(nz/2-disp))**2*exp(-((qz(k)-qz(nz/2-disp))/ &
                        (grad*dz_grid(nz/2-disp)))**2)/(mion*(1.0+(amp-1.0)*exp(-((qz(k)-qz(nz/2-disp))/(grad*dz_grid(nz/2-disp)))**2)))
                        
                  
            enddo
            enddo
            enddo
            Ni_init = Ni_tot
            call count_ppc()
      
      end subroutine load_const_ppc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine load_den_grad(begin,mass,mratio)
            use mult_proc, only: my_rank, procnum
            use dimensions
            use boundary
            use inputs, only: PI, vsw, dx, dy, km_to_m, beta_particle, kboltz, mion, amp, grad, nf_init,b0_init,mu0,ppc
            use grid, only: qx,qy,qz,dz_grid,dy_grid,dx_grid
            use gutsp
            use var_arrays, only: np,vp,vp1,xp,input_p,up,Ni_tot,Ni_tot_sys,input_E,ijkp,m_arr,mrat,beta,beta_p,wght,grav,temp_p
            use misc
            implicit none
            integer(4), intent(in):: begin 
            real, intent(in):: mratio, mass
                                  
            integer:: disp,ppcpp,extra_ppc
            integer(4):: Ni_tot_1
            real:: vx, vy, vz, vth, Temp, Tempcalc,vol,new
            integer:: l,m,i,j,k
            
            
            disp = 0 !Displacement of gradient
            beta = real(ppc)/(nf_init*dx_grid(1)*dy_grid(1)*dz_grid(1))
        !    vth=sqrt(pbeta*va**2)
            
            write(*,*) 'step 4', beta
            
            Ni_tot_1 = begin
            do i=1,nx-2
            do j=1,ny-2
            do k=1,nz-2
            vol = dx_grid(i)*dy_grid(j)*dz_grid(k)  !km^3
            ppcpp = int(ppc/procnum*beta*(nf_init+nf_init*(amp-1.0)*exp(-((qz(k)-qz(nz/2-disp))/(grad*dz_grid(nz/2-disp)))**2))*vol)
            write(*,*) 'step 5',ppcpp 
            new = 0.0
            if (new .lt. 1.0) then
                  if (new .gt. pad_ranf()) then
                        new = 1.0
                  else 
                        new = 0.0
                  endif
            endif
            extra_ppc = nint(new)
            
            do l = Ni_tot_1, Ni_tot_1+ppcpp+extra_ppc-1
                  xp(l,1) = qx(i)+(1.0-pad_ranf())*(qx(i+1)-qx(i))
                  xp(l,2) = qy(j)+(1.0-pad_ranf())*(qy(j+1)-qy(j))
                  xp(l,3) = qz(k)+(1.0-pad_ranf())*(qz(k+1)-qz(k))
                  m_arr(l) = mass
                  mrat(l) = mratio
                  beta_p(l) = 1.0
                  
                  vx = 40.0*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf()) !remember to add in vsw to get the flow velocity
                  vy = 40.0*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
                  vz = 40.0*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
                  
!                  ii = ijkp(l,1)
!                  kk = ijkp(l,3)
                  
!                  vp(l,1) = -0.0*(exp(-(xp(l,3)-qz(nz/2))**2/(10.*delz)**2)
!               x        *exp(-(xp(l,1)-qx(nx/2))**2/(10.*dx)**2))+vx
                  vp(l,1) = vx+57.0!*exp(-(xp(l,3)-qz(nz/2))**2/(30*dz_grid(nz/2))**2) !Gaussian velocity perturbation (20)
                  vp(l,2) = vy 
                  vp(l,3) = vz 
                  
                  
                  
            enddo
            
            Ni_tot_1 = Ni_tot_1+ppcpp+extra_ppc
            Ni_tot=Ni_tot_1-1
            
            enddo
            enddo
            enddo
            
            
            Ni_tot_sys = Ni_tot*procnum
            
            do l = begin,Ni_tot
                  call get_pindex(i,j,k,l)
                  do m=1,3
                        vp1(l,m) = vp(l,m)
                        input_E = input_E + 0.5*m_arr(l)*(vp(l,m)*km_to_m)**2/(beta * beta_p(l))
                        input_p(m) = input_p(m) + m_arr(l) * vp(l,m) / (beta * beta_p(l))
                  enddo
            enddo
            
            
            call get_interp_weights()
            call update_np()
            call update_up(vp)
            
            call get_temperature()
            Tempcalc = sum(temp_p(2,2,1:(nz-1)))/1e6/(nz-1) !in kg km^2/s^2
            do i=1,nx
                  do j=1,ny
                        do k=1,nz
                        grav(i,j,k) = -2.0*Tempcalc*(amp-1.0)*(qz(k)-qz(nz/2-disp))/(grad*dz_grid(nz/2-disp))**2*exp(-((qz(k)-qz(nz/2-disp))/ &
                              (grad*dz_grid(nz/2-disp)))**2)/(mion*(1.0+(amp-1.0)*exp(-((qz(k)-qz(nz/2-disp))/(grad*dz_grid(nz/2-disp)))**2)))
                        enddo
                  enddo
            enddo
     !       call count_ppc()
            
      end subroutine load_den_grad
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      subroutine load_aniso_Maxwellian(vth,Ni_tot_1)
            use dimensions
            use boundary
            use grid, only: qx,qy,qz,dz_grid
            use inputs, only: mion, dx, dy,delz, vsw, km_to_m, PI, beta_particle
            use gutsp
            use var_arrays, only: np,vp,vp1,xp,input_p,up,Ni_tot,input_E,ijkp,m_arr,mrat,beta,beta_p,wght
            implicit none
            real, intent(in):: vth
            integer(4), intent(in):: Ni_tot_1
            real:: vx,vy,vz
            real:: aniso_frac, vthx, vthy, vthz
            integer:: i,j,k,l,m
            
            aniso_frac = 0.06
            
            
            do l = 1, Ni_tot_1
                  xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
                  xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
                  xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz-1)-qz(1))
                  m_arr(l) = mion
                  mrat(l) = 1.0
                  beta_p(l) = beta_particle
                  
                  call get_pindex(i,j,k,l)
                  
                  vx = vsw+vth*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
                  vy = vth*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
                  vz = vth*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())

!                  ii = ijkp(l,1)
!                  kk = ijkp(l,3)
                  vp(l,1) = -0.0*(exp(-(xp(l,3)-qz(nz/2))**2/(10.*delz)**2) &
                        *exp(-(xp(l,1)-qx(nx/2))**2/(10.*dx)**2))+vx
                  vp(l,2) = vy 
                  vp(l,3) = vz 
                  
                  do m = 1,3
                        vp1(l,m) = vp(l,m)
                        input_E = input_E + 0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 / (beta * beta_p(l))
                        input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / (beta * beta_p(l))
                  enddo
                  
            enddo
            
            vthx = 1200.0
            vthy = 1200.0
            vthz = 500.0
            
            do l = Ni_tot_1+1, Ni_tot_1 + aniso_frac*Ni_tot_1
                  xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
                  xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
                  xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz-1)-qz(1))
                  m_arr(l) = mion
                  mrat(l) = 1.0
                  beta_p(l) = beta_particle
                  
                  call get_pindex(i,j,k,l)
                  
                  vx = vsw+vthx*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
                  vy = vthy*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
                  vz = vthz*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())

!                  ii = ijkp(l,1)
!                  kk = ijkp(l,3)
                  vp(l,1) = -0.0*(exp(-(xp(l,3)-qz(nz/2))**2/(10.*delz)**2) &
                        *exp(-(xp(l,1)-qx(nx/2))**2/(10.*dx)**2))+vx
                  vp(l,2) = vy 
                  vp(l,3) = vz 
                  
                  do m = 1,3
                        vp1(l,m) = vp(l,m)
                        input_E = input_E + 0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 / (beta * beta_p(l))
                        input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / (beta * beta_p(l))
                  enddo
                  
            enddo
            
            Ni_tot = Ni_tot_1 + aniso_frac*Ni_tot_1
            
            do l = Ni_tot + 1, Ni_tot + aniso_frac*Ni_tot
                  xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
                  xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
                  xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz-1)-qz(1))
                  m_arr(l) = 2.0*mion
                  mrat(l) = 0.5
                  beta_p(l) = beta_particle
                  
                  call get_pindex(i,j,k,l)
                  
                  vx = vsw+vth*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
                  vy = vth*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
                  vz = vth*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())

!                  ii = ijkp(l,1)
!                  kk = ijkp(l,3)
                  vp(l,1) = -0.0*(exp(-(xp(l,3)-qz(nz/2))**2/(10.*delz)**2) &
                        *exp(-(xp(l,1)-qx(nx/2))**2/(10.*dx)**2))+vx
                  vp(l,2) = vy 
                  vp(l,3) = vz 
                  
                  do m = 1,3
                        vp1(l,m) = vp(l,m)
                        input_E = input_E + 0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 / (beta * beta_p(l))
                        input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / (beta * beta_p(l))
                  enddo
                  
            enddo
            
            Ni_tot = Ni_tot + aniso_frac*Ni_tot
            
            call get_interp_weights()
            call update_np()
            call update_up(vp)
      
      end subroutine load_aniso_Maxwellian
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine sw_part_setup_maxwl_mix()
            use dimensions
            use inputs, only: mion, vth_bottom
            use gutsp
            use var_arrays, only: np,vp,vp1,xp,input_p,up,np_t_flg,np_b_flg,Ni_tot,input_E,ijkp,m_arr,mrat,beta,beta_p,wght
            implicit none
            integer(4):: Ni_tot_1
            
            
            np_t_flg = 0
            np_b_flg = 0
            
!       Add cold populations first

            Ni_tot_1 = 1
            call load_Maxwellian(vth_bottom,Ni_tot_1,mion,1.0) !mass ratio
! add He++

!      Ni_tot_1 = Ni_tot + 1
!      Ni_tot = 2.0*Ni_tot_0
      
!      call load_Maxwellian(np,vp,vp1,xp,input_p,up,
!     x     vth_bottom, Ni_tot_1, 2.0*mion, 1.0/2.0, 10.0) !inverse ration, this is a mass 2 particle (q/m)
         

! add pickup distribution

!         Ni_tot_1 = Ni_tot + 1
!         Ni_tot = 3*Ni_tot_0

!         do 69 l = Ni_tot_1,Ni_tot
!
!            xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
!            xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
!
!           xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz-1)-qz(1))
!
!            m_arr(l) = mion
!            mrat(l) = 1.0
!            beta_p(l) = beta_pu

!            i=0
! 71         continue
!            i = i + 1
!            if (xp(l,1) .gt. qx(i)) go to 71 !find i on non-uniform 
!            i = i-1
!            ijkp(l,1)= i


!            ijkp(l,2) = floor(xp(l,2)/dy) 
            
!            k=0
! 70         continue
!            k = k + 1
!            if (xp(l,3) .gt. qz(k)) go to 70 !find k on non-uniform 
!            k = k-1
!            ijkp(l,3)= k

!            theta = pad_ranf()*2*PI
            
!            vp(l,1) = vsw+vsw*cos(theta) !+dvx
!            vp(l,2) = vsw*sin(theta) !+dvz 
!            vp(l,3) = 0.0

!            if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
!            if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
            
!            do  m=1,3
!               vp1(l,m) = vp(l,m)
!               input_E = input_E + 
!     x              0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /(beta*beta_p(l))
!               input_p(m) = input_p(m) + m_arr(l)*vp(l,m)/(beta*beta_p(l))
!            enddo
            

! 69      enddo
            call get_interp_weights()
            call update_np()
            call update_up(vp)
            
      end subroutine sw_part_setup_maxwl_mix
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine load_foreshock_Maxwellian(vth,Ni_tot_1,Ni_tot_2,mass,mratio,population0,swfsRatio)
use dimensions
use boundary
use inputs, only: PI, vsw, dx, dy, km_to_m, beta_particle, kboltz, mion, amp, grad, nf_init,b0_init,mu0,boundx, Lo, q, mO, va_f, delz,ddthickness, plasma_beta,FSBeamWidth, FSThermalRatio,ddthickness,FSDriftSpeed,ForeshockBeta, TDpressureBalance
use grid, only: qx,qy,qz,dz_grid
use gutsp
use var_arrays, only: np,vp,vp1,xp,input_p,up,Ni_tot,input_E,ijkp,m_arr,mrat,beta,beta_p,wght,grav,temp_p,mix_ind,b0, b1,E,ExB,bt,up_cold
implicit none
integer(4), intent(in):: Ni_tot_1, Ni_tot_2,population0
real, intent(in):: mratio, mass, vth,swfsRatio
real:: Lo_y,eoverm
integer:: population
integer:: disp,montecarlo
real:: vth2, vx, vy, vz, va, va_x, Temp, Tempcalc, pl_beta(nx,ny,nz), fitdist,randpop, EvBx, EvBy,EvBz,theta2,phi2, vxring, vyring, vzring,vring,vr,vtheta
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



if (population0 .eq. 1) then
	randpop = (1.0-pad_ranf())
	if (randpop .le. swfsRatio) then
		population = 1
	else
		population = 2
	endif
else
	population = population0
endif

    ! X Component
    if ((population .eq. 1) .or. (population .eq. 5) )then !Solar Wind
    
        mix_ind(l) = 0
        
        if (boundx .eq. 5) then
        xp(l,1) = 0.0*qx(1)+(1.0-pad_ranf())*(qx(2)-qz(1)) +(1.0)*(qx(2)-qz(1)) !Shock real, to propagate SW/field
        endif
        if (boundx .eq. 4) then
        xp(l,1) = 1.0*qx(1)+(1.0-pad_ranf())*(qx(2)-qz(1)) !Normal runs, open boundary
        endif
        
        !xp(l,1) = (1.0-pad_ranf())*(qx(1))!+qx(1)) !too significant beam
        
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
    if ((population .eq. 2) .or. (population .eq. 7) .or. (population .eq. 3) ) then !Middle Beam
    	!!xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz)-qz(1))
        !xp(l,3) = (qz(nz/2)+qz(1)-ddthickness*(qz(2)-qz(1))) - (1.0-pad_ranf())*(float(FSBeamWidth)*(qz(2)-qz(1)))
        !!!xp(l,3) = (qz(nz/2)-qz(1)-ddthickness*(qz(2)-qz(1))) - (1.0-pad_ranf())*(float(FSBeamWidth)*(qz(2)-qz(1))) !Old below TD
        !xp(l,3) = (qz(nz/2-1)-2*ddthickness*(qz(2)-qz(1))) - (1.0-pad_ranf())*(float(FSBeamWidth)*(qz(2)-qz(1)))
        
        !below TD.
!*        if (population .eq. 3) then
!*         xp(l,3) = (qz(nz/2-1)) - (1.0-pad_ranf())*(float(FSBeamWidth)*(qz(2)-qz(1)))
!*         endif
        
        
        
        
        
        
        
        
        !foreshock ions within TD have slower Vx.
	if ((population .eq. 2) .or. (population .eq. 7)) then
         montecarlo = 0
        
         do 22 while (montecarlo .eq. 0)
        
         	xp(l,3) = qz((nz)/2) - qz(1)*FSBeamWidth +(1.0-pad_ranf())*(qz(1)*FSBeamWidth)
         	call get_pindex(i,j,k,l)
!        	k=k+int(ddthickness/2)!ddthickness/2
!         	fitdist = ( b0(i,j,k,1) / (  sqrt( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) ) )
!         	!fitdist = tanh( ( qz((nz-2)/2.0)-xp(l,3))/((1.0/3.0)*ddthickness*delz))
		!!*fitdist = (FSDriftSpeed*va_f*( b0(i,j,k-int(TDpressureBalance),1) / (  sqrt( b0(i,j,k-int(TDpressureBalance),1)**2+b0(i,j,k-int(TDpressureBalance),2)**2+b0(i,j,k-int(TDpressureBalance),3)**2 ) )  ) - ( (1.0-TDpressureBalance)*(va_f)*b0(i,j,k-int(TDpressureBalance),2) )*b0(i,j,k-int(TDpressureBalance),2) / ( ( b0(i,j,k-int(TDpressureBalance),1)**2+b0(i,j,k-int(TDpressureBalance),2)**2+b0(i,j,k-int(TDpressureBalance),3)**2 ) ) ) /(FSDriftSpeed*va_f)
		fitdist = (   FSDriftSpeed*va_f*( b0(i,j,k+int(TDpressureBalance),1) / (  sqrt( b0(i,j,k+int(TDpressureBalance),1)**2+b0(i,j,k+int(TDpressureBalance),2)**2+b0(i,j,k+int(TDpressureBalance),3)**2 ) )  ) - ( (1.0-0*TDpressureBalance)*(va_f)*b0(i,j,k+int(TDpressureBalance),2)*b0(i,j,k+int(TDpressureBalance),2) / ( ( b0(i,j,k+int(TDpressureBalance),1)**2+b0(i,j,k+int(TDpressureBalance),2)**2+b0(i,j,k+int(TDpressureBalance),3)**2 ) ) ) /(FSDriftSpeed*va_f))
		fitdist = (   ( b0(i,j,k+8*int(TDpressureBalance),1) / (  sqrt( b0(i,j,k+8*int(TDpressureBalance),1)**2+b0(i,j,k+8*int(TDpressureBalance),2)**2+b0(i,j,k+8*int(TDpressureBalance),3)**2 ) )  ) - ( (1.0-0*TDpressureBalance)*b0(i,j,k+8*int(TDpressureBalance),2)*b0(i,j,k+8*int(TDpressureBalance),2) / ( ( b0(i,j,k+8*int(TDpressureBalance),1)**2+b0(i,j,k+8*int(TDpressureBalance),2)**2+b0(i,j,k+8*int(TDpressureBalance),3)**2 ) ) ))
		fitdist= b0(i,j,k+0*int(TDpressureBalance),1)/sqrt( b0(i,j,k+0*int(TDpressureBalance),1)**2+b0(i,j,k+0*int(TDpressureBalance),2)**2+b0(i,j,k+0*int(TDpressureBalance),3)**2 )! - 1.0*b0(i,j,k+int(TDpressureBalance),2) *b0(i,j,k+int(TDpressureBalance),2) / ( ( b0(i,j,k+int(TDpressureBalance),1)**2+b0(i,j,k+int(TDpressureBalance),2)**2+b0(i,j,k+int(TDpressureBalance),3)**2 ) ) 
		!write(*,*) 'k,fitdist', k,fitdist
		!write(*,*) 'dist1,dist2,ratio', b0(i,j,k+int(TDpressureBalance),1)/ (  sqrt( b0(i,j,k+int(TDpressureBalance),1)**2+b0(i,j,k+int(TDpressureBalance),2)**2+b0(i,j,k+int(TDpressureBalance),3)**2 ) ) , (1.0-0*TDpressureBalance)*b0(i,j,k+int(TDpressureBalance),2)*b0(i,j,k+int(TDpressureBalance),2) / ( ( b0(i,j,k+int(TDpressureBalance),1)**2+b0(i,j,k+int(TDpressureBalance),2)**2+b0(i,j,k+int(TDpressureBalance),3)**2 ) ), (   ( b0(i,j,k+int(TDpressureBalance),1) / (  sqrt( b0(i,j,k+int(TDpressureBalance),1)**2+b0(i,j,k+int(TDpressureBalance),2)**2+b0(i,j,k+int(TDpressureBalance),3)**2 ) )  ) - ( (1.0-0*TDpressureBalance)*b0(i,j,k+int(TDpressureBalance),2)*b0(i,j,k+int(TDpressureBalance),2) / ( ( b0(i,j,k+int(TDpressureBalance),1)**2+b0(i,j,k+int(TDpressureBalance),2)**2+b0(i,j,k+int(TDpressureBalance),3)**2 ) ) ))
		!*fitdist = (   FSDriftSpeed*va_f*( b0(i,j,k,1) / (  sqrt( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) )  ) - ( (1.0-TDpressureBalance)*(va_f)*b0(i,j,k,2)*b0(i,j,k,2) / ( ( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) ) ) /(FSDriftSpeed*va_f))
		!fitdist = (FSDriftSpeed*va_f*( b0(i,j,k-0,1) / (  sqrt( b0(i,j,k-0,1)**2+b0(i,j,k-0,2)**2+b0(i,j,k-0,3)**2 ) )  ) - ( TDpressureBalance*(va_f)*b0(i,j,k-0,2) )*b0(i,j,k-0,2) / ( ( b0(i,j,k-0,1)**2+b0(i,j,k-0,2)**2+b0(i,j,k-0,3)**2 ) ))/(FSDriftSpeed*va_f) 

		!fitdist = b0(i,j,k,1) /  sqrt( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 )  + b0(i,j,k,2) / sqrt( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 )
                    !endif
!        	!write(*,*) 'k,fitdist',k,fitdist,b0(nx/2,j,k,1) ,sqrt( b0(nx/2,j,k,1)**2+b0(nx/2,j,k,2)**2+b0(nx/2,j,k,3)**2 )
         	if (pad_ranf() .le. fitdist) montecarlo = 1
!
!
!        	
!        
         22 continue
        endif
        
!        else if ((population .eq. 8) )then !TD for foreshock ions
!        !Fill in a bit more in solar wind.
!        montecarlo = 0
!        do 20 while (montecarlo .eq. 0) 
!        	xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz-1)-qz(1))
!        	if ( xp(l,3) .le. qz(nz/2) ) then
!        	fitdist = 1 - 0.5 - 0.5*tanh( (qz(nz/2) - xp(l,3)) / (ddthickness*delz) )
!        	else if ( xp(l,3) .gt. qz(nz/2) ) then
!        	fitdist = 1 - 0.5 - 0.5*tanh( (xp(l,3)- qz(nz/2)) / (ddthickness*delz) )
!        	endif
!        	!fitdist = cosh( (xp(l,3) - qz(nz/2)) / (0.5*ddthickness) )**(-2)
!        	if (pad_ranf() .le. fitdist) montecarlo = 1
!        20 continue
!        !write(*,*) 'z...........',xp(l,3)/101
        
        
        
        
        !Gradient within TD
!*        montecarlo = 0
!*        do 21 while (montecarlo .eq. 0)
!*        
!*        	xp(l,3) = qz(nz/2)-qz(1) - qz(1)*FSBeamWidth + (1.0-pad_ranf())*(qz(1)*FSBeamWidth)
!*        	if (xp(l,3) .le. qz(nz/2-2*ddthickness) .or. population .eq. 3) then
!*        		montecarlo=1
!*        	endif
!*        	if (xp(l,3) .gt. qz(nz/2-2*ddthickness)) then
!*			fitdist = (qz(nz/2+ddthickness*2)-xp(l,3))/(qz(1)*4*ddthickness)
!*			!call get_pindex(i,j,k,l)
!*			!write(*,*) 'k,fitdist',k,fitdist
!*       		if (pad_ranf() .le. fitdist) montecarlo = 1
!*        	endif

        	
        
!*        21 continue
        
        !montecarlo = 0
       ! do 21 while (montecarlo .eq. 0)
        !	xp(l,3) = qz(nz/2) - qz(1)*FSBeamWidth +(1.0-pad_ranf())*(qz(1)*FSBeamWidth)
       ! 	!write(*,*) 'xp(l,3),qz(nz/2),qz(1)*FSBeamWidth...........',qz(nz/2) ,- qz(1)*FSBeamWidth
       ! 	fitdist = 1.0 - (qz(nz/2)-xp(l,3))/(qz(1)*FSBeamWidth)
       ! 	!write(*,*) 'fitdist...........',fitdist
       ! 	if (pad_ranf() .le. fitdist) montecarlo = 1
       ! 
       ! 21 continue
        
        
        
        
        
        !Foreshock Bubble
        !xp(l,3) = qz(3*nz/4) - (1.0-pad_ranf())*(float(FSBeamWidth)*(qz(2)-qz(1)))
   ! else if (population .eq. 3) then !Middle Beam, with soft edges
    !	!!xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz)-qz(1))
    !    xp(l,3) = (qz(nz/2)-qz(1)-ddthickness*(qz(2)-qz(1))) - (1.0-pad_ranf())*(float(FSBeamWidth)*(qz(2)-qz(1)))
!
   !     !montecarlo = 0
    !    !do 20 while (montecarlo .eq. 0) !Monte Carlo to fit particles onto initial condition
    !    !    xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz)-qz(1))
    !    !    fitdist = exp( - ( ( xp(l,3) - ( qz(nz)-qz(1) )  /2.0 ) / ( 50*(qz(2)-qz(1)) ) )**10 )
    !    !    if (pad_ranf() .le. fitdist) montecarlo = 1
    !    !20 continue
!
    else if ((population .eq. 4) .or. (population .eq. 5))then !TD
        !Fill in a bit more in solar wind.
        montecarlo = 0
        do 20 while (montecarlo .eq. 0) 
        	xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz-1)-qz(1))
        	if ( xp(l,3) .le. qz(nz/2) ) then
        	fitdist = 1 - 0.5 - 0.5*tanh( (qz(nz/2) - xp(l,3)) / (ddthickness*delz) )
        	else if ( xp(l,3) .gt. qz(nz/2) ) then
        	fitdist = 1 - 0.5 - 0.5*tanh( (xp(l,3)- qz(nz/2)) / (ddthickness*delz) )
        	endif
        	!fitdist = cosh( (xp(l,3) - qz(nz/2)) / (0.5*ddthickness) )**(-2)
        	if (pad_ranf() .le. fitdist) montecarlo = 1
        20 continue
        !write(*,*) 'z...........',xp(l,3)/101
        
    else if (population .eq. 6) then !Top FS
            
            
            montecarlo = 0
        
         do 23 while (montecarlo .eq. 0)
        	xp(l,3) = (qz(nz/2)-qz(1)) + (1.0-pad_ranf())*(float(FSBeamWidth)*(qz(2)-qz(1)))
         	!xp(l,3) = qz((nz)/2) - qz(1)*FSBeamWidth +(1.0-pad_ranf())*(qz(1)*FSBeamWidth)
         	call get_pindex(i,j,k,l)

		!fitdist = (   FSDriftSpeed*va_f*( b0(i,j,k+int(TDpressureBalance),1) / (  sqrt( b0(i,j,k+int(TDpressureBalance),1)**2+b0(i,j,k+int(TDpressureBalance),2)**2+b0(i,j,k+int(TDpressureBalance),3)**2 ) )  ) - ( (1.0-0*TDpressureBalance)*(va_f)*b0(i,j,k+int(TDpressureBalance),2)*b0(i,j,k+int(TDpressureBalance),2) / ( ( b0(i,j,k+int(TDpressureBalance),1)**2+b0(i,j,k+int(TDpressureBalance),2)**2+b0(i,j,k+int(TDpressureBalance),3)**2 ) ) ) /(FSDriftSpeed*va_f))
		!fitdist = (   ( b0(i,j,k+8*int(TDpressureBalance),1) / (  sqrt( b0(i,j,k+8*int(TDpressureBalance),1)**2+b0(i,j,k+8*int(TDpressureBalance),2)**2+b0(i,j,k+8*int(TDpressureBalance),3)**2 ) )  ) - ( (1.0-0*TDpressureBalance)*b0(i,j,k+8*int(TDpressureBalance),2)*b0(i,j,k+8*int(TDpressureBalance),2) / ( ( b0(i,j,k+8*int(TDpressureBalance),1)**2+b0(i,j,k+8*int(TDpressureBalance),2)**2+b0(i,j,k+8*int(TDpressureBalance),3)**2 ) ) ))
		fitdist= b0(i,j,k+0*int(TDpressureBalance),2)/sqrt( b0(i,j,k+0*int(TDpressureBalance),1)**2+b0(i,j,k+0*int(TDpressureBalance),2)**2+b0(i,j,k+0*int(TDpressureBalance),3)**2 )!
         	if (pad_ranf() .le. fitdist) montecarlo = 1        
         23 continue
            
            
            
            
            
    else !All z
        xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz-1)-qz(1))
    endif

    m_arr(l) = mass
    mrat(l) = mratio

    
!if (xp(l,1) .le. qx(1)) then
!write(*,*) 'x...........',xp(l,1),ijkp(l,1)
	!endif


    va = b0_init/sqrt(mu0*mion*nf_init/1e9)/1e3
    vth2 = vth!


    vx = vth2*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
    vy = vth2*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
    vz = vth2*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())

 
	eoverm = q/mO
    if ((population .eq. 1) .or. (population .eq. 5) .or. (population .eq. 0) .or. (population .eq. 4)) then !Solar Wind, Left Edge
    !	if (xp(l,3) .lt. qz(nz/3)) then
        vp(l,1) = va_f*va+vx
        
      !  else if (xp(l,3) .gt. qz(2*nz/3)) then
      !  vp(l,1) = va_f*va+vx
        
       ! else
     !   vp(l,1) = 0.8*va_f*va+vx
     !   endif
        
        
        vp(l,2) = vy
        vp(l,3) = vz !- 1.0*va
        
        
        
        
        !ExB
        EvBx = -( ( 0*b0(i,j,k,3) - 0*b0(i,j,k,2) )  )
    	EvBy =  ( ( va_f*va*b0(i,j,k,3) - 0*b0(i,j,k,1) )  )
    	EvBz = -( ( va_f*va*b0(i,j,k,2) - 0*b0(i,j,k,1) )  )
        
        !vp(l,1) = vp(l,1) + ( EvBy*b0(i,j,k,3) - EvBz*b0(i,j,k,2) ) / ( ( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) )
        !vp(l,2) = vp(l,2) - ( EvBx*b0(i,j,k,3) - EvBz*b0(i,j,k,1) ) / ( ( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) )
        !vp(l,3) = vp(l,3) + ( EvBx*b0(i,j,k,2) - EvBy*b0(i,j,k,1) ) / ( ( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) )
        
        !Tilted SW Velocity
 !       vp(l,1) = vx + va_f*va*cos(15.0/180.0*PI)
!        vp(l,2) = vy + va_f*va*sin(15.0/180.0*PI)
!        vp(l,3) = vz 
	beta_p(l) = beta_particle
    else if ((population .eq. 2) .or. (population .eq. 3) .or. (population .eq. 7) .or. (population .eq. 6)) then !Foreshock Ions, Right Edge
!	vp(l,1) = -FSDriftSpeed*va_f*va + FSThermalRatio*vx 
      !  vp(l,2) = FSThermalRatio*vy                       
       ! vp(l,3) = FSThermalRatio*vz                      
    
    	call get_pindex(i,j,k,l)
    	!write(*,*) 'k', i
    	!i=int(nx/2)
    	i=i-1
    	EvBx = -( ( up(i,j,k,2)*b1(i,j,k,3) - up(i,j,k,3)*b1(i,j,k,2) )  )
    	EvBy =  ( ( up(i,j,k,1)*b1(i,j,k,3) - up(i,j,k,3)*b1(i,j,k,1) )  )
    	EvBz = -( ( up(i,j,k,1)*b1(i,j,k,2) - up(i,j,k,2)*b1(i,j,k,1) )  )
    	!ExB Drift using the local E/B
  !  	vp(l,1) = -FSDriftSpeed*va_f*va+FSThermalRatio*vx + ( up(i,j,k,1)*b1(i,j,k,2) - up(i,j,k,2)*b1(i,j,k,1) )*b1(i,j,k,2) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
   !     vp(l,2) = FSThermalRatio*vy                       - ( -up(i,j,k,1)*b1(i,j,k,2) + up(i,j,k,2)*b1(i,j,k,1) )*b1(i,j,k,1) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
  !      vp(l,3) = FSThermalRatio*vz                       + ( up(i,j,k,3)*( b1(i,j,k,2)*b1(i,j,k,2) + b1(i,j,k,1)*b1(i,j,k,1) ) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
    	
  
      	!ExB Drift using the initial ExB..constant distribution for injection
     ! 	vp(l,1) = -FSDriftSpeed*va_f*va*( b0(i,j,k,1) / (  sqrt( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) ) ) + FSThermalRatio*vx + ExB(i,j,k,1)
     !   vp(l,2) = -FSDriftSpeed*va_f*va*( b0(i,j,k,2) / (  sqrt( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) ) ) + FSThermalRatio*vy + ExB(i,j,k,2)
     !   vp(l,3) = 												    + FSThermalRatio*vz + ExB(i,j,k,3)                 
  	!write(*,*)
      	
      	

      	
      	
      	!Constant solar wind bulk flow for ExB
    	!*vp(l,1) = -FSDriftSpeed*va_f*va*( b0(i,j,k,1) / (  sqrt( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) ) ) + FSThermalRatio*vx + ( (va_f*va)*b0(i,j,k,2) )*b0(i,j,k,2) / ( ( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) )
        !*vp(l,2) = -FSDriftSpeed*va_f*va*( b0(i,j,k,2) / (  sqrt( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) ) ) + FSThermalRatio*vy - ( (va_f*va)*b0(i,j,k,2) )*b0(i,j,k,1) / ( ( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) )
        !*vp(l,3) = FSThermalRatio*vz    
        
               
        !Current B (local) for ExB using current (local) ion flow. (2022)
        vp(l,1) = -FSDriftSpeed*va_f*va*( b1(i,j,k,1) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + FSThermalRatio*vx*1.0 + ( EvBy*b1(i,j,k,3) - EvBz*b1(i,j,k,2) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
        vp(l,2) = -FSDriftSpeed*va_f*va*( b1(i,j,k,2) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + FSThermalRatio*vy - ( EvBx*b1(i,j,k,3) - EvBz*b1(i,j,k,1) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
        vp(l,3) = -FSDriftSpeed*va_f*va*( b1(i,j,k,3) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + FSThermalRatio*vz + ( EvBx*b1(i,j,k,2) - EvBy*b1(i,j,k,1) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
        
        !Temperature Anisotrophy
      !  vp(l,1) = vp(l,1) - 4*vx*  ( b1(i,j,k,1) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) 
        
       ! vp(l,2) = vp(l,2) - 4*vx*  ( b1(i,j,k,2) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) 
        
       ! vp(l,3) = vp(l,3) - 4*vx*  ( b1(i,j,k,3) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) 
 
     
        !Gyrophase Bunched with local B !Updated 12/22/2021 sqrt(8.0*8.0-4.0*4.0)
!        vp(l,1) = -FSDriftSpeed*va_f*va*( b1(i,j,k,1) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + FSThermalRatio*vx + ( EvBy*b1(i,j,k,3) - EvBz*b1(i,j,k,2) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
!        vp(l,2) = -FSDriftSpeed*va_f*va*( b1(i,j,k,2) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + FSThermalRatio*vy + 0*sqrt(8.0*8.0-4.0*4.0)*vth2 - ( EvBx*b1(i,j,k,3) - EvBz*b1(i,j,k,1) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
!        vp(l,3) = -FSDriftSpeed*va_f*va*( b1(i,j,k,3) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + FSThermalRatio*vz + 6.0*vth2		       + ( EvBx*b1(i,j,k,2) - EvBy*b1(i,j,k,1) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
            
     

        
        
        
        !Gyrophase Bunched with local B
!       vp(l,1) = -FSDriftSpeed*va_f*va*( b1(i,j,k,1) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + FSThermalRatio*vx 				    + ( EvBy*b1(i,j,k,3) - EvBz*b1(i,j,k,2) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
!        vp(l,2) = -FSDriftSpeed*va_f*va*( b1(i,j,k,2) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + FSThermalRatio*vy - sqrt(2.0)/2.0*6*vth2 - ( EvBx*b1(i,j,k,3) - EvBz*b1(i,j,k,1) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
!        vp(l,3) = -FSDriftSpeed*va_f*va*( b1(i,j,k,3) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + FSThermalRatio*vz + sqrt(2.0)/2.0*6*vth2 + ( EvBx*b1(i,j,k,2) - EvBy*b1(i,j,k,1) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
        
	!Gyrophase Bunched with constant B
!	EvBx = -( ( 0.0*b0(i,j,k,3) - 0.0*b0(i,j,k,2) )  )
!    	EvBy =  ( ( va_f*va*b0(i,j,k,3) - 0.0*b0(i,j,k,1) )  )
!    	EvBz = -( ( va_f*va*b0(i,j,k,2) - 0.0*b0(i,j,k,1) )  )
	
!        vp(l,1) = -FSDriftSpeed*va_f*va*( b0(i,j,k,1) / (  sqrt( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) ) ) + FSThermalRatio*vx                                     + ( EvBy*b0(i,j,k,3) - EvBz*b0(i,j,k,2) ) / ( ( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) )
!        vp(l,2) = -FSDriftSpeed*va_f*va*( b0(i,j,k,2) / (  sqrt( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) ) ) + FSThermalRatio*vy + sqrt(2.0)/2.0*FSThermalRatio*vth2 - ( EvBx*b0(i,j,k,3) - EvBz*b0(i,j,k,1) ) / ( ( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) )
!        vp(l,3) = -FSDriftSpeed*va_f*va*( b0(i,j,k,3) / (  sqrt( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) ) ) + FSThermalRatio*vz + sqrt(2.0)/2.0*FSThermalRatio*vth2 + ( EvBx*b0(i,j,k,2) - EvBy*b0(i,j,k,1) ) / ( ( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) )
        
        

             
        
        !Shell distribution beam +exB
        !theta2 = pad_ranf()*2*PI
        !phi2   = pad_ranf()*PI   
        !vp(l,1) = -FSDriftSpeed*va_f*va*( b1(i,j,k,1) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + FSThermalRatio*vth2*cos(theta2)*sin(phi2) + ( EvBy*b1(i,j,k,3) - EvBz*b1(i,j,k,2) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) 
        !vp(l,2) = -FSDriftSpeed*va_f*va*( b1(i,j,k,2) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + FSThermalRatio*vth2*sin(theta2)*sin(phi2) - ( EvBx*b1(i,j,k,3) - EvBz*b1(i,j,k,1) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
        !vp(l,3) = -FSDriftSpeed*va_f*va*( b1(i,j,k,3) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + FSThermalRatio*vth2*cos(phi2) + ( EvBx*b1(i,j,k,2) - EvBy*b1(i,j,k,1) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
        
        
        
      !!!!!!!!!  
        
        !Ring beam +exB Updated 3/19/2022 for vperp parameter scan, constant B
 !       i=i-1
 !   	EvBx = -( ( up(i,j,k,2)*b0(i,j,k,3) - up(i,j,k,3)*b0(i,j,k,2) )  )
 !   	EvBy =  ( ( up(i,j,k,1)*b0(i,j,k,3) - up(i,j,k,3)*b0(i,j,k,1) )  )
 !   	EvBz = -( ( up(i,j,k,1)*b0(i,j,k,2) - up(i,j,k,2)*b0(i,j,k,1) )  )
        
!	vring = 10.0
!	vr = sqrt(vy*vy + vz*vz)
!	vtheta = atan(vz/vy)
	
!	vyring = (FSThermalRatio*vr+vring*vth2)*cos(vtheta)
!	vzring = (FSThermalRatio*vr+vring*vth2)*sin(vtheta)
!	vy = sign(vyring,vy)
!	vz = vzring
	        
!	vp(l,1) = -FSDriftSpeed*va_f*va*( b1(i,j,k,1) / (  sqrt( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) ) ) + 0*FSThermalRatio*vx + ( EvBy*b0(i,j,k,3) - EvBz*b0(i,j,k,2) ) / ( ( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) )
!        vp(l,2) = -FSDriftSpeed*va_f*va*( b1(i,j,k,2) / (  sqrt( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) ) ) + vy                - ( EvBx*b0(i,j,k,3) - EvBz*b0(i,j,k,1) ) / ( ( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) )
!        vp(l,3) = -FSDriftSpeed*va_f*va*( b1(i,j,k,3) / (  sqrt( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) ) ) + vz                + ( EvBx*b0(i,j,k,2) - EvBy*b0(i,j,k,1) ) / ( ( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) )
        
     !!!!!!!!!1
        
	!Ring beam +exB Updated 12/22/2021
!	vring = sqrt(8.0*8.0-4.0*4.0)
!	vr = sqrt(vy*vy + vz*vz)
!	vtheta = atan(vz/vy)
	
!	vyring = (FSThermalRatio*vr+vring*vth2)*cos(vtheta)
!	vzring = (FSThermalRatio*vr+vring*vth2)*sin(vtheta)
!	vy = sign(vyring,vy)
!	vz = vzring
	        
!	vp(l,1) = -FSDriftSpeed*va_f*va*( b1(i,j,k,1) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + FSThermalRatio*vx + ( EvBy*b1(i,j,k,3) - EvBz*b1(i,j,k,2) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
 !       vp(l,2) = -FSDriftSpeed*va_f*va*( b1(i,j,k,2) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + vy - ( EvBx*b1(i,j,k,3) - EvBz*b1(i,j,k,1) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
!        vp(l,3) = -FSDriftSpeed*va_f*va*( b1(i,j,k,3) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + vz + ( EvBx*b1(i,j,k,2) - EvBy*b1(i,j,k,1) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
        
        
        !Ring Beam, Vth=0, Vperp = 8
!      	vring = 8.0
!	vr = sqrt(vy*vy + vz*vz)
!	vtheta = atan(vz/vy)
	
!	vyring = (vring*vth2)*cos(vtheta)
!	vzring = (vring*vth2)*sin(vtheta)
!	vy = sign(vyring,vy)
!	vz = vzring
	
!        vp(l,1) = -FSDriftSpeed*va_f*va*( b1(i,j,k,1) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + 0  + ( EvBy*b1(i,j,k,3) - EvBz*b1(i,j,k,2) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
!        vp(l,2) = -FSDriftSpeed*va_f*va*( b1(i,j,k,2) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + vy - ( EvBx*b1(i,j,k,3) - EvBz*b1(i,j,k,1) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
!        vp(l,3) = -FSDriftSpeed*va_f*va*( b1(i,j,k,3) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + vz + ( EvBx*b1(i,j,k,2) - EvBy*b1(i,j,k,1) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
        
	
	
	
	!Ring beam
!        vxring = 0*0.5*FSThermalRatio*vth2 *( b1(i,j,k,1) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) )
!        vyring = 0.5*1.0*FSThermalRatio*vth2 !1.5*1.0
!        vzring = 0.5*1.0*FSThermalRatio*vth2 !1.5*1.0	


!        montecarlo = 0
!        do 24 while (montecarlo .eq. 0)
!        	vx = vth2*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
!        	if (FSThermalRatio*vx .ge. vxring) montecarlo = 1 !cut-off is half of thermal.
!        24 continue
        
!        montecarlo = 0
!        do 25 while (montecarlo .eq. 0)
!        	vy = vth2*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
!        	if (abs(FSThermalRatio*vy) .ge. abs(vyring)) montecarlo = 1 !cut-off is half of thermal.
!        25 continue
        
!        montecarlo = 0
!        do 26 while (montecarlo .eq. 0)
!        	vz = vth2*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
!        	if (abs(FSThermalRatio*vz) .ge. abs(vzring)) montecarlo = 1 !cut-off is half of thermal.
!        26 continue
        
!        vp(l,1) = -FSDriftSpeed*va_f*va*( b1(i,j,k,1) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + FSThermalRatio*vx + ( EvBy*b1(i,j,k,3) - EvBz*b1(i,j,k,2) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
!        vp(l,2) = -FSDriftSpeed*va_f*va*( b1(i,j,k,2) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + FSThermalRatio*vy - ( EvBx*b1(i,j,k,3) - EvBz*b1(i,j,k,1) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
!        vp(l,3) = -FSDriftSpeed*va_f*va*( b1(i,j,k,3) / (  sqrt( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ) ) + FSThermalRatio*vz + ( EvBx*b1(i,j,k,2) - EvBy*b1(i,j,k,1) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
        
        
        
        
        
        
        
        
                         
  	!write(*,*) 'thermal speed, parallel speed, ExB speed...',up(i,j,k,3),EvBx*b1(i,j,k,2)/ ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ),- EvBy*b1(i,j,k,1)/ ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
  	
  	
    	!ExB Drift using the initial ExB..constant distribution for injection
!    	vp(l,1) = -FSDriftSpeed*va_f*va + FSThermalRatio*vx + ( (va_f*va+vth)*b0(nx/2,j,k,2) - vth*b0(nx/2,j,k,1) )*b0(nx/2,j,k,2) / ( ( b0(nx/2,j,k,1)**2+b0(nx/2,j,k,2)**2+b0(nx/2,j,k,3)**2 ) )
!        vp(l,2) = FSThermalRatio*vy                       - ( -(va_f*va+vth)*b0(nx/2,j,k,2) + vth*b0(nx/2,j,k,1) )*b0(nx/2,j,k,1) / ( ( b0(nx/2,j,k,1)**2+b0(nx/2,j,k,2)**2+b0(nx/2,j,k,3)**2 ) )
!        vp(l,3) = FSThermalRatio*vz                       + ( vth*( b0(nx/2,j,k,2)*b0(nx/2,j,k,2) + b0(nx/2,j,k,1)*b0(nx/2,j,k,1) ) ) / ( ( b0(nx/2,j,k,1)**2+b0(nx/2,j,k,2)**2+b0(nx/2,j,k,3)**2 ) )
    	
    	
    	!vp(l,1) = -FSDriftSpeed*va_f*va+FSThermalRatio*vx + ( EvBy*b1(i,j,k,3) - EvBz*b1(i,j,k,2) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
        !vp(l,2) = FSThermalRatio*vy                       - ( EvBx*b1(i,j,k,3) - EvBz*b1(i,j,k,1) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
        !vp(l,3) = FSThermalRatio*vz                       + ( EvBx*b1(i,j,k,2) - EvBy*b1(i,j,k,1) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
    	
        !vp(l,1) = -FSDriftSpeed*va_f*va+FSThermalRatio*vx !+ up(i,j,k,1) !+ ( E(i,j,k,2)*b1(i,j,k,3) - E(i,j,k,3)*b1(i,j,k,2) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
        !vp(l,2) = FSThermalRatio*vy                       !+ up(i,j,k,2) !- ( E(i,j,k,1)*b1(i,j,k,3) - E(i,j,k,3)*b1(i,j,k,1) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
        !vp(l,3) = FSThermalRatio*vz                       !+ up(i,j,k,3) !+ ( E(i,j,k,1)*b1(i,j,k,2) - E(i,j,k,2)*b1(i,j,k,1) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
        beta_p(l) = (ForeshockBeta)*beta_particle !2 for generating half of z side and 10 for 10% of FS ions.
	!write(*,*) 'vx,vy,vz', ( EvBy*b1(i,j,k,3) - EvBz*b1(i,j,k,2) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ), - ( EvBx*b1(i,j,k,3) - EvBz*b1(i,j,k,1) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) ),+ ( EvBx*b1(i,j,k,2) - EvBy*b1(i,j,k,1) ) / ( ( b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2 ) )
	!write(*,*) 'ux,uy,uz', up(i,j,k,1), up(i,j,k,2), up(i,j,k,3)
   ! else if ((population .eq. 0) .or. (population .eq. 4) ) then !Solar Wind Ions, Initially Everywhere
    !
    !    vp(l,1) = va_f*va+vx
    !    vp(l,2) = vy
    !    vp(l,3) = vz
    !    beta_p(l) = beta_particle

    !else if (population .eq. 3) then !Foreshock Ions, Initially Beam Middle
    !
    !    vp(l,1) = -FSDriftSpeed*va_f*va+FSThermalRatio*vx
    !    vp(l,2) = FSThermalRatio*vy
    !    vp(l,3) = FSThermalRatio*vz !1.0*va_f*va+va_f*vz
    !    beta_p(l) = 40*beta_particle !2 for generating half of z side and 10 for 10% of FS ions.
        
    !else if (population .eq. 16) then !Foreshock Ions, Top half above TD
    

        !vp(l,1) = -FSDriftSpeed*va_f*va*sqrt(2.0)/2.0 + FSThermalRatio*vx + (E(i,j,k,2)*b1(i,j,k,3) - E(i,j,k,3)*b1(i,j,k,2))/(  (b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2)  )
        !vp(l,2) = -FSDriftSpeed*va_f*va*sqrt(2.0)/2.0 + FSThermalRatio*vy - (E(i,j,k,1)*b1(i,j,k,3) - E(i,j,k,3)*b1(i,j,k,1))/((b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2))
        !vp(l,3) = FSThermalRatio*vz + (E(i,j,k,1)*b1(i,j,k,2) - E(i,j,k,2)*b1(i,j,k,1))/((b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2))!1.0*va_f*va+va_f*vz
        !beta_p(l) = (ForeshockBeta)*beta_particle !2 for generating half of z side and 10 for 10% of FS ions.
        
       ! write(*,*) 'vx,vy,vz,vexbx, vexb2, vexb3',vx, vy, vz,(E(i,j,k,2)*b1(i,j,k,3) - E(i,j,k,3)*b1(i,j,k,2)), - (E(i,j,k,1)*b1(i,j,k,3) - E(i,j,k,3)*b1(i,j,k,1)), (E(i,j,k,1)*b1(i,j,k,2) - E(i,j,k,2)*b1(i,j,k,1))
       !write(*,*) 'b0_init,b1...', sqrt(b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2), b1(i,j,k,1),b1(i,j,k,2),sqrt(b1(i,j,k,1)**2 + b1(i,j,k,2)**2 ),(E(i,j,k,2)*b1(i,j,k,3) - E(i,j,k,3)*b1(i,j,k,2))/((b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2)),(E(i,j,k,1)*b1(i,j,k,3) - E(i,j,k,3)*b1(i,j,k,1))/((b1(i,j,k,1)**2+b1(i,j,k,2)**2+b1(i,j,k,3)**2))
    endif

 

    do m=1,3
        vp1(l,m) = vp(l,m)
        input_E = input_E + 0.5*m_arr(l)*(vp(l,m)*km_to_m)**2/(beta * beta_p(l))
        input_p(m) = input_p(m) + m_arr(l) * vp(l,m) / (beta * beta_p(l))
    enddo

enddo




 call get_interp_weights()

 call update_np()

 call update_up(vp)





do i=1,nx
do j=1,ny
do k=1,nz
grav(i,j,k) = 0.0
enddo
enddo
enddo





end subroutine load_foreshock_Maxwellian
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calculateTotalPressure()
use inputs, only: mu0,q, mion,vth_bottom, mO
use dimensions, only: nx,ny,nz
use var_arrays, only: temp_p,np, btc, bt, Ptherm, PB, Ptotal, Ptotal_sum, Ptotal_average, Ni_tot
use mult_proc, only: my_rank
use gutsp, only: get_temperature
real:: moverq
!integer(4), intent(in)::
!integer(4), intent(out)::
integer:: i,j,k, cell_count

 cell_count = 0
moverq= mO/q
Ptotal_sum = 0
 call get_temperature()


do i=2,nx-1
do j=2,ny-1
do k=2,nz-1

Ptherm(i,j,k) = 2.0*abs(temp_p(i,j,k))*np(i,j,k)*1.0e-9

PB(i,j,k) = (moverq*moverq*bt(i,j,k,1)*bt(i,j,k,1) + moverq*moverq*bt(i,j,k,2)*bt(i,j,k,2) + moverq*moverq*bt(i,j,k,3)*bt(i,j,k,3) )/(2.0*mu0)

Ptotal(i,j,k) = Ptherm(i,j,k) + PB(i,j,k)
Ptotal_sum= Ptotal_sum + Ptotal(i,j,k)
 cell_count = cell_count + 1

enddo
enddo
enddo
Ptotal_average = Ptotal_sum/( cell_count )

end subroutine calculateTotalPressure

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine balanceTotalPressure(population)
use inputs, only: mu0,q, mion,vth_bottom, mO, ddthickness, dx
use dimensions, only: nx,ny,nz
use var_arrays, only: temp_p,np, b0, btc, bt, Ptherm, PB, Ptotal, additional_ions, beta, Ptotal_sum, Ptotal_average,Ni_tot, density_sum,density_average, sumAddedPerRow, avgAddedPerRow
use mult_proc, only: my_rank
use gutsp, only: get_temperature
real:: moverq, PthermPP,PSW_average,den_part,Ptherm_sum,PB_sum,Ptherm_average,PB_average
integer(4), intent(in):: population
!integer(4), intent(out)::
integer:: i,j,k, cell_count, addions, Ni_tot_1, Ni_tot_2


moverq= mO/q
den_part = 1/(beta*dx**3)







!if (my_rank .eq. 0) then
!write(*,*) 'PtotalAverage...PB...Ptherm...PtotalCellDD', Ptotal_average, PB(3,2,nz/2), Ptherm(3,2,nz/2), Ptotal(3,2,nz/2), i, j, k
!write(*,*) 'End'
!endif


!If the cell is below average, add particles //(vth,Ni_tot_1,Ni_tot_2,mass,mratio,population)







!!Check pressure balance
!call get_temperature()
!
!
!do i=1,nx
!do j=1,ny
!do k=1,nz
!
!Ptherm(i,j,k) = 2.0*abs(temp_p(i,j,k))*np(i,j,k)*1.0e-9
!
!PB(i,j,k) = (moverq*moverq*bt(i,j,k,1)*bt(i,j,k,1) + moverq*moverq*bt(i,j,k,2)*bt(i,j,k,2) + moverq*moverq*bt(i,j,k,3)*bt(i,j,k,3) )/(2.0*mu0)
!
!Ptotal(i,j,k) = Ptherm(i,j,k) + PB(i,j,k)
!Ptotal_average = Ptotal_average + Ptotal(i,j,k)
!cell_count = cell_count + 1
!
!
!
!
!!write(*,*) 'temp_p(i,j,k)...np(i,j,k)...bt(i,j,k,1)...bt(i,j,k,2)...bt(i,j,k,3)', temp_p(i,j,k), np(i,j,k), moverq*bt(i,j,k,1), moverq*bt(i,j,k,2), moverq*bt(i,j,k,3)
!!write(*,*) 'PtotalAverage...PB...Ptherm...PtotalCellDD...cell_count', Ptotal_average, PB(i,j,k), Ptherm(i,j,k), Ptotal(i,j,k), cell_count
!
!
!enddo
!enddo
!enddo
!
!!Calculate average pressure throughout simulation
!Ptotal_average = Ptotal_average/( cell_count )










if (population .eq. 0) then !Initialization Pressure Balance

density_sum = 0
density_average = 0
Ptotal_sum = 0
Ptotal_average = 0
 cell_count = 0
Ptherm_sum = 0
PB_sum = 0

do k=1,nz
sumAddedPerRow(k)=0
avgAddedPerRow(k)=0
enddo

!Check pressure balance
 call get_temperature()


do i=3,nx-2
do j=2,ny-1
do k=3,nz-2

!Ptherm(i,j,k) = 2.0*abs(temp_p(i,j,k))*np(i,j,k)*1.0e-9

Ptherm(i,j,k) = 2.0*abs(temp_p(i,j,k))*np(i,j,k)*1.0e-9!*1.602e-19
PB(i,j,k) = (moverq*moverq*b0(i,j,k,1)*b0(i,j,k,1) + moverq*moverq*b0(i,j,k,2)*b0(i,j,k,2) + moverq*moverq*b0(i,j,k,3)*b0(i,j,k,3) )/(2.0*mu0)

Ptotal(i,j,k) = Ptherm(i,j,k) + PB(i,j,k)

if ((k .lt. (nz/2 - 3*ddthickness) ) .or. (k .gt. (nz/2 + 3*ddthickness))) then

Ptotal_sum = Ptotal_sum + Ptotal(i,j,k)
Ptherm_sum = Ptherm_sum + Ptherm(i,j,k)
PB_sum = PB_sum + PB(i,j,k)
 cell_count = cell_count + 1
!density_sum = density_sum + np(i,j,k)

else
endif

!write(*,*) 'temp_p(i,j,k)...np(i,j,k)...bt(i,j,k,1)...bt(i,j,k,2)...bt(i,j,k,3)', temp_p(i,j,k), np(i,j,k), moverq*bt(i,j,k,1), moverq*bt(i,j,k,2), moverq*bt(i,j,k,3)
!write(*,*) 'PtotalAverage...PB...Ptherm...PtotalCellDD...cell_count', Ptotal_average, PB(i,j,k), Ptherm(i,j,k), Ptotal(i,j,k), cell_count


enddo
enddo
enddo

!Calculate average pressure throughout simulation
Ptotal_average = Ptotal_sum/( cell_count )
Ptherm_average = Ptherm_sum/cell_count
PB_average = PB_sum / cell_count
!density_average = density_sum/cell_count







do i=2,nx-1
do j=2,ny-1
do k=nz/2 - 3*ddthickness,nz/2 + 3*ddthickness
   additional_ions=0

!
do while (Ptotal(i,j,k) .lt. 1.0*Ptotal_average)
!do while (np(i,j,k) .lt. 0.99*density_average)
!            PthermPP = (Ptherm(i,j,k)/(np(i,j,k)/den_part))
!            addions = nint((Ptotal_average - Ptotal(i,j,k))/PthermPP)
!            write(*,*) 'temp_p(i,j,k),np(i,j,k)', temp_p(i,j,k),np(i,j,k)
!            write(*,*) 'Ptherm(i,j,k)', Ptherm(i,j,k)
!	     write(*,*) 'Ptherm(i,j,k), PB(i,j,k)', Ptherm(i,j,k),PB(i,j,k)
!            write(*,*) 'den_part,beta', den_part,beta
!            write(*,*) 'Pdiff', Ptotal_average,Ptotal(i,j,k)
!            write(*,*) 'PthermPP,addions', PthermPP,addions

            Ni_tot_1 = Ni_tot +1
            Ni_tot_2 = Ni_tot + 1!addions
            Ni_tot = Ni_tot_2
            call add_ion(vth_bottom,Ni_tot_1,Ni_tot_2,mion,1.0,0,i,j,k) !Population either 1 or 0
            additional_ions = additional_ions+1!addions
	    sumAddedPerRow(k)= sumAddedPerRow(k) + 1 
	    !write(*,*) 'sumAddedPerRow', sumAddedPerRow(k)

!            write(*,*) 'added 1 ion at...Ptotal at cell, PB, Ptherm, PtotalAverage', i,j,k, Ptotal(i,j,k), PB(i,j,k), Ptherm(i,j,k), Ptotal_average
!            write(*,*) '# of ions in cell', np(i,j,k),Ni_tot


            !Check Pressure again after adding particle
            call get_temperature()
            !Ptherm(i,j,k) = 2.0*abs(temp_p(i,j,k))*np(i,j,k)*1e-9
            Ptherm(i,j,k) = 2.0*abs(temp_p(i,j,k))*np(i,j,k)*1e-9!*1.602e-19

!            PB(i,j,k) = ( moverq*moverq*bt(i,j,k,1)*bt(i,j,k,1) + moverq*moverq*bt(i,j,k,2)*bt(i,j,k,2) + moverq*moverq*bt(i,j,k,3)*bt(i,j,k,3) )/(2.0*mu0)

            Ptotal(i,j,k) = Ptherm(i,j,k) + PB(i,j,k)
        enddo
if (my_rank .eq. 0) then
    !write(*,*) 'Ptotal(i,j,k)...Ptotal_average...', Ptotal(i,j,k),Ptotal_average,i,j,k,additional_ions
endif
        !write(*,*) 'total added ions in this cell...', additional_ions
enddo
enddo
enddo


do k=1,nz

avgAddedPerRow(k) =sumAddedPerRow(k)/(nx-2)
write(*,*) 'k,sum,avg', k,sumAddedPerRow(k),avgAddedPerRow(k)
enddo 

!write(*,*) 'Ptotal_average,Ptherm_average,PB_average', Ptotal_average,Ptherm_average,PB_average




else if (population .eq. 1) then !New Solar Wind Particles Pressure Balance
PSW_average = 0
 cell_count=0


!Check pressure balance
call get_temperature()


do i=1,2
do j=2,ny-1
do k=2,nz-1

!Ptherm(i,j,k) = 2.0*abs( temp_p(2,j,k)     )*np(i,j,k)*1.0e-9

Ptherm(i,j,k) = 2.0*abs(temp_p(i,j,k))*np(i,j,k)*1e-9!*1.602e-19
PB(i,j,k) = (moverq*moverq*b0(i,j,k,1)*b0(i,j,k,1) + moverq*moverq*b0(i,j,k,2)*b0(i,j,k,2) + moverq*moverq*b0(i,j,k,3)*b0(i,j,k,3) )/(2.0*mu0)
!Ptherm(i,j,k) = 2.0*abs(temp_p(1,j,k)+temp_p(2,j,k))*(np(1,j,k)+np(2,j,k))*1e-9/4.0!*1.602e-19
!PB(1,j,k) = (moverq*moverq*b0(1,j,k,1)*b0(1,j,k,1) + moverq*moverq*b0(1,j,k,2)*b0(1,j,k,2) + moverq*moverq*b0(1,j,k,3)*b0(1,j,k,3) )/(2.0*mu0)

Ptotal(i,j,k) = Ptherm(i,j,k) + PB(i,j,k)

if ((k .lt. (nz/2 - 3*ddthickness) ) .or. (k .gt. (nz/2 + 3*ddthickness))) then

PSW_average = PSW_average + Ptotal(i,j,k)
 cell_count = cell_count + 1
else

endif



!write(*,*) 'temp_p(i,j,k)...np(i,j,k)...bt(i,j,k,1)...bt(i,j,k,2)...bt(i,j,k,3)', temp_p(i,j,k), np(i,j,k), moverq*bt(i,j,k,1), moverq*bt(i,j,k,2), moverq*bt(i,j,k,3)
!write(*,*) 'PtotalAverage...PB...Ptherm...PtotalCellDD...cell_count', Ptotal_average, PB(i,j,k), Ptherm(i,j,k), Ptotal(i,j,k), cell_count


enddo
enddo
enddo

!Calculate average pressure throughout simulation
!!PSW_average = PSW_average/( cell_count )






!write(*,*) 'Ptotal_average', Ptotal_average




do i=2,2
do j=2,ny-1
do k=nz/2 - 3*ddthickness, nz/2 + 3*ddthickness

		Ni_tot_1 = Ni_tot + 1
            Ni_tot_2 = Ni_tot + avgAddedPerRow(k)/2! addions
            Ni_tot = Ni_tot_2!addions
            
            call add_ion(vth_bottom,Ni_tot_1,Ni_tot_2,mion,1.0,1,i,j,k) !Population either 1 or 0
            
            
        	!if (my_rank .eq. 0) then
         		!write(*,*) 'k,avgAddedPerRow...',k,avgAddedPerRow(k)
       		!endif




!!!!additional_ions=0
!!!!!do while (Ptotal(i,j,k) .lt. 1.0*Ptotal_average) !2.4
!do while (np(i,j,k) .lt. 0.99*density_average)
            !!!!!PthermPP = (Ptherm(i,j,k)/(np(i,j,k)/den_part))
!            addions = nint((Ptotal_average - Ptotal(i,j,k))/PthermPP)
            !addions = nint((PSW_average - (PB(1,j,k) + Ptherm(i,j,k)))/PthermPP)
            !!!!!!addions = nint((Ptotal_average - (PB(1,j,k) + Ptherm(i,j,k)))/PthermPP)
!            write(*,*) 'temp_p(i,j,k),np(i,j,k)', temp_p(i,j,k),np(i,j,k)
!            write(*,*) 'Ptherm(i,j,k)', Ptherm(i,j,k)
!            write(*,*) 'den_part,beta', den_part,beta
!            write(*,*) 'Pdiff', Ptotal_average,Ptotal(i,j,k)
!            write(*,*) 'PthermPP,addions', PthermPP,addions
 
            !!!!!!Ni_tot_1 = Ni_tot +1
            !!!!!!Ni_tot_2 = Ni_tot + 1! addions
            !!!!!!Ni_tot = Ni_tot_2!addions
            
            !!!!!!call add_ion(vth_bottom,Ni_tot_1,Ni_tot_2,mion,1.0,1,i,j,k) !Population either 1 or 0
            !!!!!!additional_ions = additional_ions+1!addions

!            if (my_rank .eq. 0)  then
!                write(*,*) 'PB, PTherm, Ptotal,PSW_average...', PB(3,j,k),Ptherm(i,j,k),Ptotal(i,j,k),PSW_average,i,j,k,additional_ions,Ni_tot_1,Ni_tot_2,Ni_tot
!                write(*,*) 'Ni_tot_1,Ni_tot_2,Ni_tot...', Ni_tot_1,Ni_tot_2,Ni_tot
!            endif

!            write(*,*) 'added 1 ion at...Ptotal at cell, PB, Ptherm, PtotalAverage', i,j,k, Ptotal(i,j,k), PB(i,j,k), Ptherm(i,j,k), Ptotal_average


            !Check Pressure again after adding particle
            !!!!!call get_temperature()
            !!!!!Ptherm(i,j,k) = 2.0*abs(temp_p(i,j,k))*np(i,j,k)*1e-9!abs(   (temp_p(1,j,k)+temp_p(2,j,k))/2.0      )*np(i,j,k)*1e-9*2.0
	     !Ptherm(i,j,k) = 2.0*abs(temp_p(1,j,k)+temp_p(2,j,k))*(np(1,j,k)+np(2,j,k))*1e-9/4.0!*1.602e-19
!            PB(i,j,k) = ( moverq*moverq*btc(i,j,k,1)*btc(i,j,k,1) + moverq*moverq*btc(i,j,k,2)*btc(i,j,k,2) + moverq*moverq*btc(i,j,k,3)*btc(i,j,k,3) )/(2.0*mu0)

            !!!!!Ptotal(i,j,k) = Ptherm(i,j,k) + PB(i,j,k)
!!!!enddo
        !write(*,*) 'total added ions in this cell...', additional_ions
!write(*,*) 'PB, PTherm,additional_ions',PB(i,j,k),Ptherm(i,j,k),additional_ions,Ptotal(i,j,k),Ptotal_average

!if (my_rank .eq. 0) then
!write(*,*) 'PB, PTherm, Ptotal,PSW_average...', PB(3,j,k),Ptherm(i,j,k),Ptotal(i,j,k),PSW_average,i,j,k,additional_ions
!endif

enddo
enddo
enddo
!Ni_tot = Ni_tot + additional_ions
!Ptotal(i,j,k),PSW_average,i,j,k,additional_ions,Ni_tot_1,Ni_tot_2,Ni_tot
endif












end subroutine balanceTotalPressure
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine add_ion(vth,Ni_tot_1,Ni_tot_2,mass,mratio,population,x,y,z)
use dimensions
use boundary
use inputs, only: PI, vsw, dx, dy, km_to_m, beta_particle, kboltz, mion, amp, grad, nf_init,b0_init,mu0,boundx, Lo, q, mO, va_f, delz
use grid, only: qx,qy,qz,dz_grid
use gutsp
use var_arrays, only: np,vp,vp1,xp,input_p,up,Ni_tot,input_E,ijkp,m_arr,mrat,beta,beta_p,wght,grav,temp_p,mix_ind,b0
use mult_proc, only: my_rank
implicit none
integer(4), intent(in):: Ni_tot_1, Ni_tot_2,population, x, y, z
real, intent(in):: mratio, mass, vth
real:: Lo_y

integer:: disp
real:: vth2, vx, vy, vz, va, va_x, Temp, Tempcalc, pl_beta(nx,ny,nz), sechdist
integer:: l,m,i,j,k

disp = 0

do i=1,nx
do j=1,ny
do k=1,nz
pl_beta(i,j,k) = 1.0
enddo
enddo
enddo

va = b0_init/sqrt(mu0*mion*nf_init/1e9)/1e3

do l = Ni_tot_1,Ni_tot_2
mix_ind(l) = 0
!if (my_rank .eq. 0) then
!write(*,*) 'mix_ind(l)...',l,mix_ind(l)
!endif
if (population .eq. 1) then

    if (x .eq. 1) then
        xp(l,1) = (1.0-pad_ranf())*(qx(1))!+qx(1)
    elseif (x .eq. 2) then
        xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(2)-qz(1))
    else
        xp(l,1) = ( qx(x) - qx(1) ) + (1.0-pad_ranf()) * ( qx(2) -qx(1) )
    endif

xp(l,2) = ( qy(y) - qy(1) ) + (1.0-pad_ranf()) * ( qy(2) -qy(1) )
xp(l,3) = ( qz(z) - qz(1) ) + (1.0-pad_ranf()) * ( qz(2) -qz(1) )

else if (population .eq. 0) then

xp(l,1) = ( qx(x) - qx(1) ) + (1.0-pad_ranf()) * ( qx(2) -qx(1) )
xp(l,2) = ( qy(y) - qy(1) ) + (1.0-pad_ranf()) * ( qy(2) -qy(1) )
xp(l,3) = ( qz(z) - qz(1) ) + (1.0-pad_ranf()) * ( qz(2) -qz(1) )

endif




vx = vth*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
vy = vth*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())
vz = vth*sqrt(-log(pad_ranf()))*cos(PI*pad_ranf())

if (population .eq. 1) then !Solar Wind Ions
vp(l,1) = va_f*va+vx
vp(l,2) = vy
vp(l,3) = vz

else if (population .eq. 0) then
vp(l,1) = va_f*va+vx
vp(l,2) = vy
vp(l,3) = vz

endif

m_arr(l) = mass
mrat(l) = mratio
beta_p(l) = beta_particle

 call get_pindex(i,j,k,l)

do m=1,3
vp1(l,m) = vp(l,m)
input_E = input_E + 0.5*m_arr(l)*(vp(l,m)*km_to_m)**2/(beta * beta_p(l))
input_p(m) = input_p(m) + m_arr(l) * vp(l,m) / (beta * beta_p(l))
enddo

enddo

 call get_interp_weights()
 call update_np()
 call update_up(vp)

do i=1,nx
do j=1,ny
do k=1,nz
grav(i,j,k) = 0.0
enddo
enddo
enddo

end subroutine add_ion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module part_init
