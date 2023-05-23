program hybrid
      use mpi
      use boundary
      use dimensions
      use Var_Arrays
      use mult_proc
      use inputs
      use misc
      use initial
      use part_init
      use gutsf
      use gutsp
      use grid_interp
      use grid
      use iso_fortran_env
      
      
      implicit none

      real:: Evp,Euf,EB1,EB1x,EB1y,EB1z,EE,EeP,input_EeP
      real:: input_chex, input_bill,dtsub
      real:: puf(3),peb(3)
      character(2):: filenum!(16) !max 16 processors
      integer:: ierr,t1,t2,cnt_rt,m,ndiag,seed,ndiag_part
      real(8):: time
      integer(4):: Ni_tot_sw, Ni_tot_1,Ni_tot_2
      integer:: n,ntf,l !looping indicies
      integer(4):: Nx_boundary, sw_delayTime, FS_boundary, FS_initial, TD_initial,TD_boundary
      real (real64) :: dp
      integer, parameter :: dp_kind = kind(dp)
      real::sw_speed
      

      call readInputs()
      call initparameters()
      
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, procnum, ierr)
      
      
      call system_clock(t1, cnt_rt)
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !Initialize all variables

      write(filenum, '(I2)') my_rank
      

      Ni_tot=int((nx-2)*(ny-2)*(nz-2)*(ppc/procnum)) !1D
      Ni_tot_0 = Ni_tot
      Ni_tot_sw = Ni_tot
      Ni_tot_sys = Ni_tot*procnum
      
      if (my_rank .eq. 0) then
            call check_inputs()
            write(*,*) 'Total particles, PPP, #pc', Ni_tot_sys,Ni_tot,procnum
            write(*,*) 'Particles per cell... ', Ni_tot_sys/((nz-2)*(ny-2)*(nx-2))
            write(*,*) ' '
      endif
      
      ndiag = 0
      ndiag_part = 0

      seed = t1 + my_rank*100
      call seed_mpi(my_rank)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      
     
      input_E = 0.0
      input_p = 0.0
      input_chex = 0.0
      input_bill = 0.0
      
      input_Eb = 0.0
      input_EeP= 0.0
      
      call grd7()
      call grd6_setup(b0,bt,b12,b1,b1p2,nu,input_Eb)
      call get_beta(Ni_tot_sys,beta)
      
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
      !Background Solar Wind
      call load_foreshock_Maxwellian(vth,1,Ni_tot,mion,0) !SW
  
      !Initial TD, TD center, turn off for real shock 1/28/22
      if (boundx .eq. 4) then
            Ni_tot_1 = Ni_tot + 1
            TD_initial = procnum*floor(float(int((nx-2)*(ny-2)*TDcellBalance*ppc/procnum))/procnum) !for 15dd, 5.06 for 15dd (30depletionwidth), then divide by two because B drops only 25%, 4.24664 for 15dd (15depletionwidith), 0.283 for 1dd, 2.83 for 20dd, 2.53 for 16dd
            Ni_tot_2 = Ni_tot_1 + TD_initial -1
            Ni_tot = Ni_tot_2  
            call load_foreshock_Maxwellian(vth,Ni_tot_1,Ni_tot_2,mion,4) !TD
      elseif (boundx .eq. 5) then
            !Real Shock, no TD balance
      endif
  
      if (my_rank .eq. 0) then
            write(*,*) 'Total Particles, TD Particles',Ni_tot, TD_initial
      end if
  


      if (my_rank .eq. 0) then
            call check_inputs()     
      endif

      call f_update_tlev(b1,b12,b1p2,bt,b0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write parameter file
      if (my_rank .eq. 0) then

            open(109, file=trim(out_dir)//'para.dat', &
                  status='unknown',form='unformatted')
      
            write(109) nx,ny,nz,dx,dy,delz
            write(109) nt,dtsub_init,ntsub,dt,nout,mion
            write(109) out_dir
            write(109) vtop,vbottom
            write(109) Ni_max
            write(109) mion
            write(109) nf_init
            write(109) b0
            write(109) vth
            write(109) alpha,beta
            close(109)
            
      endif
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Inititalize output files
      if ((my_rank .eq. 0)) then !For 0th thread
            open(110,file=trim(out_dir)//'c.np.dat',status='unknown',form='unformatted')
            open(120,file=trim(out_dir)//'c.mixed_0.dat',status='unknown',form='unformatted')
            open(130,file=trim(out_dir)//'c.b1.dat',status='unknown',form='unformatted')
            open(140,file=trim(out_dir)//'c.aj.dat',status='unknown',form='unformatted')
            open(150,file=trim(out_dir)//'c.E.dat',status='unknown',form='unformatted')
            open(160,file=trim(out_dir)//'c.energy.dat',status='unknown',form='unformatted')
            open(170,file=trim(out_dir)//'c.np_cold.dat',status='unknown',form='unformatted')
            open(175,file=trim(out_dir)//'c.temp_p_cold.dat',status='unknown',form='unformatted')
            open(176,file=trim(out_dir)//'c.np_mixed.dat',status='unknown',form='unformatted')
            open(177,file=trim(out_dir)//'c.temp_p_mixed.dat',status='unknown',form='unformatted')
            open(315,file=trim(out_dir)//'c.gradP.dat',status='unknown',form='unformatted')
            open(317,file=trim(out_dir)//'c.ue.dat',status='unknown',form='unformatted')
            
            open(180,file=trim(out_dir)//'c.up.dat',status='unknown',form='unformatted')
            open(181,file=trim(out_dir)//'c.up_cold.dat',status='unknown',form='unformatted')
            open(182,file=trim(out_dir)//'c.up_mixed.dat',status='unknown',form='unformatted')
            open(190,file=trim(out_dir)//'c.momentum.dat',status='unknown',form='unformatted')
            open(192,file=trim(out_dir)//'c.p_conserve.dat',status='unknown',form='unformatted')
            open(300,file=trim(out_dir)//'c.temp_p.dat',status='unknown',form='unformatted')
            open(301,file=trim(out_dir)//'c.tp.dat',status='unknown',form='unformatted')
            open(302,file=trim(out_dir)//'c.tp_cold.dat',status='unknown',form='unformatted')
            open(303,file=trim(out_dir)//'c.tp_mixed.dat',status='unknown',form='unformatted')
            open(305,file=trim(out_dir)//'c.xp_0.dat',status='unknown',form='unformatted')
            open(310,file=trim(out_dir)//'c.vp_0.dat',status='unknown',form='unformatted')
            !open(315,file=trim(out_dir)//'c.gradP_0.dat',status='unknown',form='unformatted')
            !open(317,file=trim(out_dir)//'c.ue_0.dat',status='unknown',form='unformatted')
            open(320,file=trim(out_dir)//'c.np_wake.dat',status='unknown',form='unformatted')
            !open(330,file=trim(out_dir)//'c.vSC_0.dat',status='unknown',form='unformatted')
            open(340,file=trim(out_dir)//'c.xp_mixed_0_.dat',status='unknown',form='unformatted')
            open(350,file=trim(out_dir)//'c.mnp.dat',status='unknown',form='unformatted')
       
            open(501,file=trim(out_dir)//'c.energy_p.dat',status='unknown',form='unformatted')

            open(601,file=trim(out_dir)//'c.np_cold_foreshock.dat',status='unknown',form='unformatted')
            open(602,file=trim(out_dir)//'c.up_cold_foreshock.dat',status='unknown',form='unformatted')
            open(603,file=trim(out_dir)//'c.temp_p_cold_foreshock.dat',status='unknown',form='unformatted')
            open(604,file=trim(out_dir)//'c.tp_cold_foreshock.dat',status='unknown',form='unformatted')


            open(611,file=trim(out_dir)//'c.np_mixed_foreshock.dat',status='unknown',form='unformatted')
            open(612,file=trim(out_dir)//'c.up_mixed_foreshock.dat',status='unknown',form='unformatted')
            open(613,file=trim(out_dir)//'c.temp_p_mixed_foreshock.dat',status='unknown',form='unformatted')
            open(614,file=trim(out_dir)//'c.tp_mixed_foreshock.dat',status='unknown',form='unformatted')
      endif
       
      if (my_rank .gt. 0) then !for Subsequent threads
            open(120,file=trim(out_dir)//'c.mixed_'//trim(filenum)//'.dat',status='unknown',form='unformatted')
            open(305,file=trim(out_dir)//'c.xp_'//trim(filenum)//'.dat',status='unknown',form='unformatted')
            open(310,file=trim(out_dir)//'c.vp_'//trim(filenum)//'.dat',status='unknown',form='unformatted')
            !open(315,file=trim(out_dir)//'c.gradP_'//trim(filenum)//'_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            !open(317,file=trim(out_dir)//'c.ue_'//trim(filenum)//'_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            !open(330,file=trim(out_dir)//'c.vSC_'//trim(filenum)//'_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
      endif
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Initialize Injection Parameters
      Nx_boundary = floor((nz-2)*(ny-2)*(1)*ppc/procnum)
      FS_boundary = (((cos(ByConeAngle/180.0*pi)*float(FSBeamWidth))/&
            (nz-2+1*TDcellBalance))*float(Nx_boundary)*(FSDensityRatio/(100.0/ForeshockBeta)))
      

      sw_speed = va_f*(b0_init/sqrt(mu0*mion*nf_init/1e9)/1e3)
      sw_delayTime = int((qx(2)-qx(1))/(sw_speed*dt))

      Nx_boundary = floor(Nx_boundary/ ((qx(2)-qx(1))/(sw_speed*dt)))
      FS_boundary = floor(FSDriftSpeed*FS_boundary/ ((qx(2)-qx(1))/(sw_speed*dt)))
      TD_boundary = floor(TD_initial/float(nx-2)/ ((qx(2)-qx(1))/(sw_speed*dt)))
      
      if (my_rank .eq. 0) then
            write(*,*) 'Nx_boundary per Time Step, Nx_boundary particles in column', &
                  Nx_boundary,4*int((nz-2)*(ny-2)*(1)*ppc/procnum)/4
            write(*,*) 'FS_boundary, FS_boundary particles in column', &
                  FS_boundary, &
                  (float(FSBeamWidth)/nz)*float(4*int((nz-2)*(ny-2)*(1)*ppc/procnum)/4)*(FSDensityRatio/(100.0/ForeshockBeta))
            write(*,*) 'TD_boundary,TD_initial', TD_boundary, TD_initial
            write(*,*) 'sw displacement in one timestep, qxdiff', sw_speed*dt, qx(2)-qx(1)
            write(*,*) 'mpertimestep, rounded', (qx(2)-qx(1))/(sw_speed*dt), sw_delayTime
            write(*,*) 'Processor Count', procnum
            write(*,*) '  '
      endif


      !MAIN LOOP
      do m = 1, nt
            if (my_rank .eq. 0) then
                  write(*,*) 'time...', m, m*dt, Ni_tot
            endif

            call get_interp_weights()

            call update_np() !np at n+1/2

            call update_up(vp) !up at n+1/2

            call get_gradP()

            call curlB(bt,np,aj)

            call edge_to_center(bt,btc)
            call extrapol_up()

            call get_Ep()


            !Injection when previously injected ions move out of cell.
            if (m .gt. 0.5*sw_delayTime ) then
                 
                  Ni_tot_1 = Ni_tot+1
                  Ni_tot_2 = Ni_tot + Nx_boundary - 1
                  Ni_tot = Ni_tot_2
                  call load_foreshock_Maxwellian(vth,Ni_tot_1,Ni_tot_2,mion,1) !SW at left boundary
            
                  if (boundx .ne. 5) then
                        Ni_tot_1 = Ni_tot+1
                        Ni_tot_2 = Ni_tot + TD_boundary -1
                        Ni_tot = Ni_tot_2
                        call load_foreshock_Maxwellian(vth,Ni_tot_1,Ni_tot_2,mion,5) !TD at left boundary
                  
                        if (quasiparallel .eq. 1) then 
                              Ni_tot_1 = Ni_tot+1
                              Ni_tot_2 = Ni_tot + FS_boundary - 1
                              Ni_tot = Ni_tot_2
                              call load_foreshock_Maxwellian(vth,Ni_tot_1,Ni_tot_2,mion,2) !FS Bot Right boundary
                        endif
                  
                        if (quasiparallel .eq. 2) then 
                              Ni_tot_1 = Ni_tot+1
                              Ni_tot_2 = Ni_tot + FS_boundary - 1
                              Ni_tot = Ni_tot_2
                              call load_foreshock_Maxwellian(vth,Ni_tot_1,Ni_tot_2,mion,2) !FS Bot Right boundary
                  
                              Ni_tot_1 = Ni_tot+1
                              Ni_tot_2 = Ni_tot + FS_boundary - 1
                              Ni_tot = Ni_tot_2
                              call load_foreshock_Maxwellian(vth,Ni_tot_1,Ni_tot_2,mion,6) !FS Top Right Boundary
                        endif
                  endif

                  call f_update_tlev(b1,b12,b1p2,bt,b0)  
                  
            endif
            


            call get_interp_weights()


            call update_np()                  !np at n+1/2

            call update_up(vp)            !up at n+1/2

            call get_gradP()
 
            call curlB(bt,np,aj)



            call edge_to_center(bt,btc)


            call extrapol_up()

            call get_Ep()

            
            call get_vplus_vminus()
            call improve_up()

            

            call get_Ep()
   


            call get_vplus_vminus()


            call get_vp_final()


            
            call move_ion_half() !1/2 step ion move to n+1/2

            call get_interp_weights()

            call update_np() !np at n+1/2

            call update_up(vp) !up at n+1/2
            
            call get_gradP()
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Subcycling loop

            dtsub = dtsub_init
            ntf=ntsub
            call check_time_step(bt,np,dtsub,ntf)
            
            do n=1,ntf
                  call curlB(bt,np,aj)
                  
                  call predict_B(b12,b1p2,bt,E,aj,up,nu,dtsub)
                  
                  call correct_B(b1,b1p2,E,aj,up,np,nu,dtsub)

                  call f_update_tlev(b1,b12,b1p2,bt,b0)
                  
            enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           
 
            call move_ion_half()       !final ion move to n+1
            call update_foreshock()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       diagnositc output
            
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            
            if (my_rank .eq. 0) then
                  write(160) m
                  write(160) input_E,input_EeP,Evp,Euf,EB1,EB1x,EB1y,EB1z,EE, &
                        EeP,input_chex,input_bill
                  write(190) m
                  write(190) puf,peb,input_p
            endif
            
            ndiag = ndiag+1
            ndiag_part = ndiag_part + 1
            if (ndiag .eq. nout) then
                  call get_temperature() 

                  call update_rho()

                  call update_mixed() !Update Foreshock Ion Labels
                  
                  call update_np_cold() !Update solar wind ion densities                
                  
                  call get_temperature_cold() !Update solar wind ion temperatures    

                  call update_np_cold_foreshock() !Update solar wind ion densities                
                  
                  call get_temperature_cold_foreshock() !Update solar wind ion temperatures                            
                  
                  call update_np_mixed() !Update foreshock ion densities           
                  
                  call get_temperature_mixed() !Update foreshock ion temperatures

                  call update_np_mixed_foreshock() !Update foreshock ion densities           
                  
                  call get_temperature_mixed_foreshock() !Update foreshock ion temperatures                  
                  
                  call update_up_cold(vp) !Update solar wind flow
                  
                  call update_up_mixed(vp) !Update foreshock flow

                  call update_up_cold_foreshock(vp) !Update solar wind flow
                  
                  call update_up_mixed_foreshock(vp) !Update foreshock flow                  

 
                  if (my_rank .eq. 0) then
                        call face_to_center(E,Ec)

                        write(110) m
                        write(110) np
                        !write(120) m
                        !write(120) mix_ind
                        write(130) m
                        write(130) bt
                        write(140) m
                        write(140) curlBcurrent!aj!*alpha*np3
                        write(150) m
                        write(150) Ec                        
                        write(170) m
                        write(170) np_cold
                        write(175) m
                        write(175) temp_p_cold/1.6e-19
                        write(176) m
                        write(176) np_mixed
                        write(177) m
                        write(177) temp_p_mixed/1.6e-19
                        write(180) m
                        write(180) up
                        write(181) m
                        write(181) up_cold
                        write(182) m
                        write(182) up_mixed
                        write(300) m
                        write(300) temp_p/1.6e-19  !output in eV
                        write(301) m
                        write(301) tp/1.6e-19
                        write(302) m
                        write(302) tp_cold/1.6e-19
                        write(303) m
                        write(303) tp_mixed/1.6e-19
                        !write(305) m
                        !write(305) xp
                        !write(310) m
                        !write(310) vp
                        write(315) m
                        write(315) gradP
                        write(317) m
                        write(317) ue
                        !write(330) m
                        !write(330) vSC
                        !write(340) m
                        !write(340) testPartPos
                        write(350) m
                        write(350) mnp

                        write(601) m
                        write(601) np_cold_foreshock
                        write(602) m
                        write(602) up_cold_foreshock
                        write(603) m
                        write(603) temp_p_cold_foreshock/1.6e-19
                        write(604) m
                        write(604) tp_cold_foreshock/1.6e-19                   

                        write(611) m
                        write(611) np_mixed_foreshock
                        write(612) m
                        write(612) up_mixed_foreshock
                        write(613) m
                        write(613) temp_p_mixed_foreshock/1.6e-19
                        write(614) m
                        write(614) tp_mixed_foreshock/1.6e-19           

                   endif
                   ndiag = 0
            endif

                  !Output Dist data
              if (ndiag_part .eq. nout) then
                   if (my_rank .ge. 0) then
                        write(120) m
                        write(120) mix_ind
                        write(305) m
                        write(305) xp
                        write(310) m
                        write(310) vp
                    endif    
                    ndiag_part = 0
            endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
      enddo     !End Main Loop
      
      close(110)
      close(120)
      close(130)
      close(140)
      close(150)
      close(160)
      close(170)
      close(172)
      close(175)
      close(176)
      close(177)
      close(180)
      close(181)
      close(182)
      close(190)
      close(192)
      close(210)
      close(220)
      close(221)
      close(300)
      close(301)
      close(302)
      close(303)
      close(305)
      close(310)
      close(315)
      close(317)
      close(320)
!      close(330)
!      close(340)
      close(342)
      close(350)
!      close(401)
!      close(402)
!      close(403)
      close(501)
      
      close(601)
      close(602)
      close(603)
      close(604)
      close(611)
      close(612)
      close(613)
      close(614)

      call system_clock(t2,cnt_rt)
      time=(real(t2,dp_kind) - real(t1,dp_kind))/real(cnt_rt,dp_kind)
      if (my_rank .eq. 0) then
            write(*,*)
            write(*,*) 'Elapsed time .....', time, ' sec'
            write(*,*)
      endif
      
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      
      call MPI_FINALIZE(ierr)

end program hybrid

      
