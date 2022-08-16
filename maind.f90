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
      use chem_rates
      use grid
      use iso_fortran_env
      
      implicit none

      real:: Evp,Euf,EB1,EB1x,EB1y,EB1z,EE,EeP,input_EeP,vB0x,vB0y,vB0z, eoverm, mO_q,b0eoverm,phi
      real:: input_chex, input_bill,dtsub
      real:: pup(3),puf(3),peb(3)
      character(2):: filenum!(16) !max 16 processors
      character(1):: mstart
      integer:: ierr,t1,t2,cnt_rt,m,mstart_n,ndiag,seed,ndiag_part
      real(8):: time
      logical:: restart = .false.
      integer(4):: Ni_tot_sw, Ni_tot_1,Ni_tot_2
      integer:: i,j,k,n,ntf,mm,tl,l !looping indicies
      integer:: N_1, N_2, Ntot_initial, Nx_boundary, sw_delayTime, FS_boundary, testParticlesAssigned,FSstartupTime, FS_initial, TD_initial,TD_boundary, TD_procCount, TD_procRand
      real (real64) :: dp
      integer, parameter :: dp_kind = kind(dp)
      real::va_x, sw_speed,swfsRatio
      
!      filenum = (/'1 ','2 ','3 ','4 ','5 ','6 ','7 ','8 ','9 ', &
!            '10','11','12','13','14','15','16'/)
                 
      call readInputs()
      call initparameters()
      
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, procnum, ierr)
      
      
      call system_clock(t1, cnt_rt)
      
!Initialize all variables

      write(filenum, '(I2)') my_rank
      

      Ni_tot=int((nx-2)*(ny-2)*(nz-2)*(ppc/procnum)) !1D
      Ni_tot_0 = Ni_tot
      Ni_tot_sw = Ni_tot
      Ni_tot_sys = Ni_tot*procnum
      
      if (my_rank .eq. 0) then
            call check_inputs()
            write(*,*) 'Total particles, PPP, #pc', Ni_tot_sys,Ni_tot,procnum
            write(*,*) 'Partilces per cell... ', Ni_tot_sys/((nz-2)*(ny-2)*(nx-2))
            write(*,*) ' '
      endif
      
      mstart_n = 0 !number of times restarted
      write(mstart, '(I1)') mstart_n
      
      ndiag = 0
      ndiag_part = 0
      prev_Etot = 1.0
!      nuei = 0.0

      seed = t1 + my_rank*100
!      call random_initialize(seed)
      call seed_mpi(my_rank)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      
      if (.not. restart) then
!            do i=1,nx
!                  do j=1,ny
!                        do k=1,nz
                              input_E = 0.0
                              input_p = 0.0
                              input_chex = 0.0
                              input_bill = 0.0
                             
                              input_Eb = 0.0
                              input_EeP= 0.0
!                        enddo
!                  enddo
!            enddo
      endif
      

      
      call grd7()
!      call grid_gaussian()
      call grd6_setup(b0,bt,b12,b1,b1p2,nu,input_Eb)
      call get_beta(Ni_tot_sys,beta)
   
      input_E = 0.0
      bndry_Eflux = 0.0
      
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
      !Initialize particles: use load Maxwellian, or sw_part_setup, etc.
!      call load_Maxwellian(vth,1,mion,1.0)
!      call load_const_ppc(vth,1,mion,1.0)
   !write(*,*) 'np,np_cold', np(50,2,50),np_cold(50,2,50)
  !Initialize particles: use load Maxwellian, or sw_part_setup, etc.
  call load_foreshock_Maxwellian(vth_bottom,1,Ni_tot,mion,1.0,0,1.0) !SW
  !write(*,*) 'Total Particles',Ni_tot
  
  
  !Initial TD
 
  !Ni_tot_2 = Ni_tot_1 + int(Ni_tot*((1.0/procnum)*ddthickness/nz))
  !TD_initial = int((float(Ni_tot)*(float(ddthickness)/float(nz))))
  !TD_initial = ceiling(float(Ni_tot)/float(nz))
!*   TD_initial = int((nx-2)*(ny-2)*2*0.3375*ppc/procnum) !for 1dd

  !Integral of 1 - 0.5 -0.5*tanh(x/ddthickness). For ddthickness = 12, 4.

  !write(*,*) 'Ni_tot_1,Ni_tot_2 after TD,TD_intial', Ni_tot_1,Ni_tot_2, TD_initial

  !write(*,*) 'Ni_tot_1,Ni_tot_2,Ni_tot_',Ni_tot_1,Ni_tot_2,Ni_tot,TD_initial
  
  
  
  !TD center, turn off for real shock 1/28/22
   Ni_tot_1 = Ni_tot + 1
     TD_initial = procnum*floor(float(int((nx-2)*(ny-2)*2*1.0*0.5*TDcellBalance*TDpressureBalance*ppc/procnum))/procnum) !for 15dd, 5.06 for 15dd (30depletionwidth), then divide by two because B drops only 25%, 4.24664 for 15dd (15depletionwidith), 0.283 for 1dd, 2.83 for 20dd, 2.53 for 16dd
   Ni_tot_2 = Ni_tot_1 + TD_initial 
  Ni_tot = Ni_tot_2  
  call load_foreshock_Maxwellian(vth_bottom,Ni_tot_1,Ni_tot_2,mion,1.0,4,1.0) !TD
  
  
  
  
  !write(*,*) 'TD Particles',int(Ni_tot*((1.0/procnum)*ddthickness/nz))
    !write(*,*) 'TD Particles',int(Ni_tot*((0.5*ddthickness/nz)))
    	!!Ni_tot_2 = Ni_tot_1 + float(Nx_boundary)*(float(ddthickness)/(nz))*float(sw_delayTime)/procnum
    	
    	
    	
  !Load Foreshock 
!*  Ni_tot_1 = Ni_tot + 1
  FS_initial = (float(FSBeamWidth)/nz)*int((nz-2)*(ny-2)*(nx-2)*ppc/procnum)*int(FSDensityRatio/(100.0/ForeshockBeta))!2.5!
!*  Ni_tot_2 = Ni_tot_1 + FS_initial
!*  Ni_tot = Ni_tot_2
!*  call load_foreshock_Maxwellian(vth_bottom,Ni_tot_1,Ni_tot_2,mion,1.0,3,1.0) !Initial FS Beam
  !write(*,*) 'FS Particles',FS_initial
  
  if (my_rank .eq. 0) then
  write(*,*) 'Total Particles, TD Particles, Foreshock Particles',Ni_tot, TD_initial, FS_initial
  end if
  
  
  !write(*,*) 'Ni_tot_1,Ni_tot_2,Ni_tot...',Ni_tot_1,Ni_tot_2,Ni_tot
 ! call balanceTotalPressure(0)
  !Ni_tot_1 = Ni_tot+1
  !Ni_tot_2 = Ni_tot + 16*int((ny-2)*(nx-2)*ppc/procnum)*(float(FSBeamWidth))/4
 ! Ni_tot_2 = Ni_tot + int(((nx-2)*(ny-2)*(FSBeamWidth)*ppc/procnum))
  !Ni_tot = Ni_tot_2

  !call load_foreshock_Maxwellian(vth_bottom,Ni_tot_1,Ni_tot_2,mion,1.0,3,0.0) !FS
  
    
   !write(*,*) 'FSbeam count', 8*int(((ny-2)*(nx-2)*ppc/procnum)*(float(FSBeamWidth)))
  
  
  
      if (my_rank .eq. 0) then
            call check_inputs()     
      endif
!      call load_ring_beam(57.0,int(Ni_tot*0.1),mion,1.0)
      
!      write(*,*) 'Particles per cell... (Ni_tot_sys/(nz-2)', Ni_tot_sys/(nz-2)
      
      call f_update_tlev(b1,b12,b1p2,bt,b0)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!  Check for restart flag

      !write(*,*) 'restart status...', restart
      if ((restart) .and. (mstart_n .gt. 0)) then
!            if (my_rank .eq. 0) then
            write(*,*) 'opening restar vars....'
            open(210,file=trim(out_dir)//'restart.vars',status='unknown',form='unformatted')
            write(*,*) 'reading restart vars......'
!            read(210) b1,b12,b1p2,bt,btmf,btc,np,np3,vp,vp1,vplus,vminus, &
!                  up,xp,aj,nu,Ep,E,temp_p,mnp,beta,beta_p,Evp,Euf, &
!                  EB1,EB1x,EB1y,EB1z,EE,EeP,input_E,Ni_tot, &
!                  ijkp,input_p,input_EeP,input_Eb,prev_Etot,bndry_Eflux,grav, &
!                  input_chex,input_bill,mrat,m_arr
                  
            read(210) b1,b12,b1p2,bt,btmf,btc,np,np3, &
                        up,aj,nu,E,temp_p,mnp,beta,Evp,Euf, &
                        EB1,EB1x,EB1y,EB1z,EE,EeP, &
                        input_EeP,input_Eb,prev_Etot,bndry_Eflux,grav, &
                        input_chex,input_bill
            write(*,*) 'restarting hybrid ....'
            
!            endif
            
!            if (my_rank .gt. 0) then
                  open(211,file=trim(out_dir)//'restart.part'//trim(filenum),status='unknown',form='unformatted')
                  read(211) vp,vp1,vplus,vminus,xp,Ep,input_E,Ni_tot,ijkp,input_p,mrat,m_arr,beta_p
!            endif
            close(210)
            close(211)
      endif

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
                  write(109) mion,m_pu,m_heavy
                  write(109) np_top,np_bottom
                  write(109) b0_top,b0_bottom
                  write(109) vth_top,vth_bottom
                  write(109) alpha,beta
                  close(109)
                  
! Write fft parameter file
!                  open(401, file = trim(out_dir)//'fft_11400.dat',status='unknown',form='unformatted')
!                  write(401) dt,nt,omega_p
                  
!                  open(402, file = trim(out_dir)//'fft_14000.dat',status='unknown',form='unformatted')
!                  write(402) dt,nt,omega_p
                  
!                  open(403, file = trim(out_dir)//'fft_17000.dat',status='unknown',form='unformatted')
!                 write(403) dt,nt,omega_p

            endif
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Inititalize diagnositc output files

      if ((my_rank .eq. 0) .and. restart) then
            open(110,file=trim(out_dir)//'c.np_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(115,file=trim(out_dir)//'c.np_b_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(120,file=trim(out_dir)//'c.mixed_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(130,file=trim(out_dir)//'c.b1_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(140,file=trim(out_dir)//'c.aj_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(150,file=trim(out_dir)//'c.E_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(160,file=trim(out_dir)//'c_.energy_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(170,file=trim(out_dir)//'c.np_cold_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(175,file=trim(out_dir)//'c.temp_p_cold_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(176,file=trim(out_dir)//'c.np_mixed_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(177,file=trim(out_dir)//'c.temp_p_mixed_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
	    open(315,file=trim(out_dir)//'c.gradP_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(317,file=trim(out_dir)//'c.ue_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            
            !diagnostics chex,bill,satnp
            
            open(180,file=trim(out_dir)//'c.up_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(181,file=trim(out_dir)//'c.up_cold_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(182,file=trim(out_dir)//'c.up_mixed_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(190,file=trim(out_dir)//'c.momentum_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(192,file=trim(out_dir)//'c.p_conserve_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(300,file=trim(out_dir)//'c.temp_p_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(301,file=trim(out_dir)//'c.tp_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(302,file=trim(out_dir)//'c.tp_cold_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(303,file=trim(out_dir)//'c.tp_mixed_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(305,file=trim(out_dir)//'c.xp_0_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(310,file=trim(out_dir)//'c.vp_0_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            !open(315,file=trim(out_dir)//'c.gradP_0_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(317,file=trim(out_dir)//'c.ue_0)'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(320,file=trim(out_dir)//'c.np_wake_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            !open(330,file=trim(out_dir)//'c.vSC_0)'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(340,file=trim(out_dir)//'c.xp_mixed_0_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(342,file=trim(out_dir)//'c.test_part_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            open(350,file=trim(out_dir)//'c.mnp_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            
       endif
       if ((my_rank .eq. 0) .and. (.not. restart)) then
            open(110,file=trim(out_dir)//'c.np.dat',status='unknown',form='unformatted')
            open(115,file=trim(out_dir)//'c.np_b.dat',status='unknown',form='unformatted')
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
            !diagnostics chex,bill,satnp
            
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
            open(342,file=trim(out_dir)//'c.test_part.dat',status='unknown',form='unformatted')
            open(350,file=trim(out_dir)//'c.mnp.dat',status='unknown',form='unformatted')
       
            open(501,file=trim(out_dir)//'c.energy_p.dat',status='unknown',form='unformatted')
       endif
       
       if (my_rank .gt. 0) then
       	    open(120,file=trim(out_dir)//'c.mixed_'//trim(filenum)//'.dat',status='unknown',form='unformatted')
            open(305,file=trim(out_dir)//'c.xp_'//trim(filenum)//'.dat',status='unknown',form='unformatted')
            open(310,file=trim(out_dir)//'c.vp_'//trim(filenum)//'.dat',status='unknown',form='unformatted')
            !open(315,file=trim(out_dir)//'c.gradP_'//trim(filenum)//'_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            !open(317,file=trim(out_dir)//'c.ue_'//trim(filenum)//'_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
            !open(330,file=trim(out_dir)//'c.vSC_'//trim(filenum)//'_'//trim(mstart)//'.dat',status='unknown',form='unformatted')
       endif
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Ntot_initial = Ni_tot
Nx_boundary = floor((nz-2)*(ny-2)*(1)*ppc/procnum)
!FS_boundary = (float((FSBeamWidth-int(ddthickness/4)))/(nz-2))*float(Nx_boundary)*int(FSDensityRatio/(100.0/ForeshockBeta))!2.5! 5.0 for 1/20 !2.50 for 1/40
!FS_boundary = (float((FSBeamWidth))/(nz-2+1*TDcellBalance*TDpressureBalance))*float(Nx_boundary)*int(FSDensityRatio/(100.0/ForeshockBeta))!2.5! 5.0 for 1/20 !2.50 for 1/40, 2.83 for 20dd
FS_boundary = ((cos(ByConeAngle/180.0*pi)*float(FSBeamWidth))/(nz-2+1*TDcellBalance*TDpressureBalance))*float(Nx_boundary)*int(FSDensityRatio/(100.0/ForeshockBeta))!2.5! 5.0 for 1/20 !2.50 for 1/40, 2.83 for 20dd
  !FS_initial = (float(FSBeamWidth)/nz)*int((nz-2)*(ny-2)*(nx-2)*ppc/procnum)*int(FSDensityRatio/(100.0/ForeshockBeta))!2.5!
  !FS_boundary = FS_initial/(nx-2)
!FS_boundary = 2*4*(float(FSBeamWidth)/nz)*float(Nx_boundary)!2*8*FSDriftSpeed*((FSBeamWidth)/nz)*Nx_boundary!8*FSDriftSpeed!2* for 20% 4* is for 40beta (2.5*4). 10%
!FS_boundary = 2*2*(float(FSBeamWidth)/nz)*float(Nx_boundary)!2*8*FSDriftSpeed*((FSBeamWidth)/nz)*Nx_boundary!8*FSDriftSpeed!2* for 20% 4* is for 40beta (2.5*4). 5%
!FS_boundary = 2*(float(FSBeamWidth)/nz)*float(Nx_boundary)!2*8*FSDriftSpeed*((FSBeamWidth)/nz)*Nx_boundary!8*FSDriftSpeed!2* for 20% 4* is for 40beta (2.5*4). 2.5%
!FS_boundary = 2*(float(FSBeamWidth)/nz)*float(Nx_boundary) !2*2.5%
nTestParticles =  FS_boundary
va = b0_init/sqrt(mu0*mion*nf_init/1e9)/1e3
!sw_speed = va_f*va*cos(15.0/180.0*PI)!+vth
sw_speed = va_f*va
sw_delayTime = ((qx(2)-qx(1))/(sw_speed*dt))!/2.0;
testParticlesAssigned = 0
FSstartupTime=int(FS_boundary/sw_delayTime)
        Nx_boundary = floor(Nx_boundary/ ((qx(2)-qx(1))/(sw_speed*dt)))
        FS_boundary = floor(FSDriftSpeed*FS_boundary/ ((qx(2)-qx(1))/(sw_speed*dt)))
        TD_boundary = floor(TD_initial/float(nx-2)/ ((qx(2)-qx(1))/(sw_speed*dt)))!!int(0.05*float(Nx_boundary)*(float(ddthickness)/(nz)))
        swfsRatio= float(Nx_boundary)/ float(( Nx_boundary + FS_boundary))
if (my_rank .eq. 0) then
!write(*,*) 'b1(x,:,:),b0', b1(1,51,2,2),b0(1,51,2,2)
    !write(*,*) 'va', va
    write(*,*) 'Nx_boundary per Time Step, Nx_boundary particles in column', Nx_boundary,4*int((nz-2)*(ny-2)*(1)*ppc/procnum)/4
    write(*,*) 'FS_boundary, FS_boundary particles in column', FS_boundary,(float(FSBeamWidth)/nz)*float(4*int((nz-2)*(ny-2)*(1)*ppc/procnum)/4)*(FSDensityRatio/(100.0/ForeshockBeta))
    write(*,*) 'TD_boundary,TD_initial', TD_boundary, TD_initial
    write(*,*) 'sw displacement in one timestep, qxdiff', sw_speed*dt, qx(2)-qx(1)
    write(*,*) 'mpertimestep, rounded', (qx(2)-qx(1))/(sw_speed*dt), sw_delayTime
    write(*,*) 'Processor Count', procnum
    write(*,*) 'FS_initial', FS_initial,(nz-2)*(nx-2)*(ny-2)
    write(*,*) '  '
endif


	!do k=1, nz
        !        write(*,*) 'b0(nx,2,k), bmag',k,b0(nx,2,k,1),b0(nx,2,k,1)**2,b0(nx,2,k,2),sqrt( b0(nx,2,k,1)**2 + b0(nx,2,k,2)**2 + b0(nx,2,k,3)**2 )
        !enddo
        
!write(*,*) 'Proc number', my_rank
	!Ni_tot_1 = Ni_tot+1
        !Ni_tot_2 = Ni_tot + Nx_boundary + FS_boundary
        !Ni_tot = Ni_tot_2

        !call load_foreshock_Maxwellian(vth,Ni_tot_1,Ni_tot_2,mion,1.0,1,swfsRatio) !SW
        do k=nz/2-2*ddthickness,nz/2+2*ddthickness
        	if (my_rank .eq. 0) then
        		!write(*,*) 'k,Bx,Btot,cos(theta)...',k,b0(nx,2,k,1),b0(nx,2,k,2),(  sqrt( b0(nx,2,k,1)**2+b0(nx,2,k,2)**2+b0(nx,2,k,3)**2 )),-FSDriftSpeed*va_f*( b0(nx,2,k,1) / (  sqrt( b0(nx,2,k,1)**2+b0(nx,2,k,2)**2+b0(nx,2,k,3)**2 ) )  ) + ( (va_f)*b0(nx,2,k,2) )*b0(nx,2,k,2) / ( ( b0(nx,2,k,1)**2+b0(nx,2,k,2)**2+b0(nx,2,k,3)**2 ) )
                endif
                
        enddo

!       MAIN LOOP
      do m = 1, nt !mstart_n+1, nt
      !write(*,*) 'Start of Main Loop.....',m
            if (my_rank .eq. 0) then
                  write(*,*) 'time...', m, m*dt, Ni_tot
            endif

do k=nz/2-2*ddthickness,nz/2+2*ddthickness
		if (my_rank .eq. 0) then
        		!write(*,*) 'gradP...',gradP(nx/2,2,k,1),gradP(nx/2,2,k,2),gradP(nx/2,2,k,3)
                endif
enddo

            call get_interp_weights()

            call update_np()                  !np at n+1/2

            call update_up(vp)            !up at n+1/2
		!if (m .eq. 1 ) then
		!	do i=1,nx
                !  	do j=1,ny
                !        do k=1,nz 
		!    		vB0x = -( ( up(i,j,k,2)*b0(i,j,k,3) - up(i,j,k,3)*b0(i,j,k,2) )  )
    		!		vB0y =  ( ( up(i,j,k,1)*b0(i,j,k,3) - up(i,j,k,3)*b0(i,j,k,1) )  )
    		!		vB0z = -( ( up(i,j,k,1)*b0(i,j,k,2) - up(i,j,k,2)*b0(i,j,k,1) )  )
    		!	enddo
    		!	enddo
    		!	enddo
    		!	do i=1,nx
                !  	do j=1,ny
                !        do k=1,nz 
    		!	    	ExB(i,j,k,1) =  ( vB0y*b0(i,j,k,3) - vB0z*b0(i,j,k,2) ) / ( ( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) )
       		!		ExB(i,j,k,2) = -( vB0x*b0(i,j,k,3) - vB0z*b0(i,j,k,1) ) / ( ( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) )
       		!		ExB(i,j,k,3) =  ( vB0x*b0(i,j,k,2) - vB0y*b0(i,j,k,1) ) / ( ( b0(i,j,k,1)**2+b0(i,j,k,2)**2+b0(i,j,k,3)**2 ) )
       		!	enddo
    		!	enddo
    		!	enddo
    	!
		!endif

            !energy diagnostics
            call get_bndry_Eflux(b1,E,bndry_Eflux)

            call Energy_diag(Evp,Euf,EB1,EB1x,EB1y,EB1z,EE,EeP)
            call get_gradP()
 
            call curlB(b0,bt,np,aj)
            call edge_to_center(bt,btc)
            call extrapol_up()
            call get_Ep()
            
            
            


if (m .gt. 0.5*sw_delayTime ) then

        !write(*,*) 'Ni_tot_1,Ni_tot_2,Ni_tot...',Ni_tot_1,Ni_tot_2,Ni_tot
            !    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       ! !do TD_procCount=1,TD_boundary
        !	TD_procRand= nint(pad_ranf()*procnum)
       ! 	write(*,*) 'TD_procRand',TD_procRand,my_rank
        ! 		if (my_rank .eq. TD_procRand) then
        ! 			Ni_tot_1 = Ni_tot + 1
        ! 			Ni_tot_2 = Ni_tot_1 
        !			Ni_tot = Ni_tot_2
        !			call load_foreshock_Maxwellian(vth_bottom,Ni_tot_1,Ni_tot_2,mion,1.0,5,1.0) !TD
        !			write(*,*) 'TD proc Count, TD_Boundary ',TD_procCount,TD_boundary,TD_procRand
        ! 		endif
!
        !!enddo
        !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        
        !Everysolar wind delay time, generate TD as a whole
        !!!!!!if (mod(m-10,sw_delayTime) == 0) then
        !write(*,*) 'right before TD injection', float(nx-2)*((qx(2)-qx(1))/(sw_speed*dt)),float(TD_initial),float(nx-2)*((qx(2)-qx(1))/(sw_speed*dt))/float(TD_initial)
       !if ((mod(m,int(TD_initial/float(nx-2)*((qx(2)-qx(1))/(sw_speed*dt)))) .eq. 0) .and. TDpressureBalance .eq. 1.0) then
!*        if ((mod(m,int(float(nx-2)*((qx(2)-qx(1))/(sw_speed*dt))/TD_initial)) .eq. 0) .and. TDpressureBalance .eq. 1.0) then
!*        Ni_tot_1 = Ni_tot + 1
  	!Ni_tot_2 = Ni_tot_1 + float(Nx_boundary)*(float(sw_delayTime)/ddthickness*float(ddthickness)/(nz))
  	!!!Ni_tot_2 = Ni_tot_1 + float(Nx_boundary)*(float(ddthickness)/(nz))*float(sw_delayTime)/procnum
!*  	Ni_tot_2 = Ni_tot_1! + TD_boundary !! +TD_initial/ int(float(nx-2)*((qx(2)-qx(1))/(sw_speed*dt))/TD_initial) !TD_boundary
!*        Ni_tot = Ni_tot_2
        !write(*,*) 'TD',Ni_tot_2 - Ni_tot_1, Nx_boundary,(0.1*ddthickness/nz)
        !write(*,*) 'Ni_tot_1,Ni_tot_2 after TD',Ni_tot_1,Ni_tot_2
        !write(*,*) 'Producing 1 particles per timestep',int(float(nx-2)*((qx(2)-qx(1))/(sw_speed*dt))/TD_initial)
!*        call load_foreshock_Maxwellian(vth_bottom,Ni_tot_1,Ni_tot_2,mion,1.0,5,1.0) !TD
        
        !!!!!!!endif
!*        endif
        
        
        
        !call balanceTotalPressure(1)
        !do i=2,2
	!do j=2,ny-1
	!do k=nz/2 - 3*ddthickness, nz/2 + 3*ddthickness
	!            Ni_tot_1 = Ni_tot + 1
       !     Ni_tot_2 = Ni_tot + avgAddedPerRow(k)/2! addions
        !    Ni_tot = Ni_tot_2!addions
        !    
        !    call add_ion(vth_bottom,Ni_tot_1,Ni_tot_2,mion,1.0,1,i,j,k) !Population either 1 or 0
        !    
        !    
        !	if (my_rank .eq. 0) then
        ! 		write(*,*) 'k,avgAddedPerRow...',k,avgAddedPerRow(k)
       	!	endif
        !            
	!enddo
	!enddo
	!enddo
	
	
        Ni_tot_1 = Ni_tot+1
        Ni_tot_2 = Ni_tot + Nx_boundary-1
        Ni_tot = Ni_tot_2
        call load_foreshock_Maxwellian(vth,Ni_tot_1,Ni_tot_2,mion,1.0,1,1.0) !SW
        
        !if (m .ge. 1*Nx*sw_delayTime) then 
    !    Ni_tot_1 = Ni_tot+1
    !    Ni_tot_2 = Ni_tot + FS_boundary
    !    Ni_tot = Ni_tot_2
   !     call load_foreshock_Maxwellian(vth,Ni_tot_1,Ni_tot_2,mion,1.0,7,0.0) !FS Bot Left
        
        Ni_tot_1 = Ni_tot+1
        Ni_tot_2 = Ni_tot + TD_boundary
        Ni_tot = Ni_tot_2
        call load_foreshock_Maxwellian(vth,Ni_tot_1,Ni_tot_2,mion,1.0,5,1.0) !TD
        
        Ni_tot_1 = Ni_tot+1
        Ni_tot_2 = Ni_tot + FS_boundary
        Ni_tot = Ni_tot_2
        call load_foreshock_Maxwellian(vth,Ni_tot_1,Ni_tot_2,mion,1.0,2,0.0) !FS Bot Right
        
        
        if (quasiparallel .eq. 2) then
        	!write(*,*) 'start parallel injection'
        	Ni_tot_1 = Ni_tot+1
        	Ni_tot_2 = Ni_tot + int(cos(magneticShear/180*pi)*float(FS_boundary))
        	!write(*,*) 'FS Top...', int(cos(magneticShear/180*pi)*float(FS_boundary))
        	Ni_tot = Ni_tot_2
!
        	call load_foreshock_Maxwellian(vth,Ni_tot_1,Ni_tot_2,mion,1.0,6,0.0) !FS Top
        endif
      

      
      
	call f_update_tlev(b1,b12,b1p2,bt,b0)
            
        
        
        
        
        
        
        
endif



!    if (m .ge. 0.5*Nx*sw_delayTime) then 
 !   	if (testParticlesAssigned .eq. 0) then
 !   		l=Ni_tot
 !   		tl = 1
!		do while (tl .le. nTestParticles_max) !tl= 1,nTestParticles
!			
!			do while (mix_ind(l) .eq. 0)
!				l=l-1
!			end do
!			testPartIndex(tl) = l
			
			
!			if (my_rank .eq. 0) then
  !  				write(*,*) 'tl', tl,tl,mix_ind(l),l
!			endif
!			l=l-1
!			tl = tl+1
!		enddo
!		do tl= 1,nTestParticles_max
!			testPartPos(tl,1) = xp(testPartIndex(tl),1)
!			testPartPos(tl,2) = xp(testPartIndex(tl),2)
!			testPartPos(tl,3) = xp(testPartIndex(tl),3)
!		enddo
!		testParticlesAssigned=1
 !   	endif!testparticle assignment
 !   endif!End delay for test particles

!update testparticle positions
!if (m .lt. 0.5*Nx*sw_delayTime) then 
!	do tl= 1,nTestParticles_max
!			testPartPos(tl,1) = 0.0!qx(nx)
!			testPartPos(tl,2) = 0.0!qy(2)
!			testPartPos(tl,3) = 0.0!qz(nz/2.0-ddthickness)
			
			!if (my_rank .eq. 0) then
    			!	write(*,*) 'tl,x,y,z', tl,testPartPos(tl,1),testPartPos(tl,2),testPartPos(tl,3)
			!endif
!	enddo
!else
!	do tl= 1,nTestParticles_max
!		if (testPartIndex(tl) .eq. 0) then

!		else
!			testPartPos(tl,1) = xp(testPartIndex(tl),1)
!			testPartPos(tl,2) = xp(testPartIndex(tl),2)
!			testPartPos(tl,3) = xp(testPartIndex(tl),3)
!		endif
			!if (my_rank .eq. 0) then
    			!	write(*,*) 'tl,x,y,z', tl,testPartPos(tl,1),testPartPos(tl,2),testPartPos(tl,3)
			!endif
!	enddo
!endif

       !if (my_rank .eq. 0) then
       !		write(*,*) 'Ni_tot, SolarWind Ions, Foreshock Ions,Diff,Ratio...', 				Ni_tot,Nx_boundary,FS_boundary,Ni_tot_2-Ni_tot_1,swfsRatio
       !endif
                        
       
            call get_interp_weights()

            call update_np()                  !np at n+1/2

            call update_up(vp)            !up at n+1/2

            !energy diagnostics
            call get_bndry_Eflux(b1,E,bndry_Eflux)

            call Energy_diag(Evp,Euf,EB1,EB1x,EB1y,EB1z,EE,EeP)
            call get_gradP()
 
            call curlB(b0,bt,np,aj)
            call edge_to_center(bt,btc)
            call extrapol_up()
            call get_Ep()

            
            call get_vplus_vminus()
            call improve_up()

            call get_Ep()
         
            call get_vplus_vminus()
            call get_vp_final()

            call move_ion_half() !1/2 step ion move to n+1/2



!Generate Particles on Left and Right Side, 1st half time step
if (m .gt. 0*sw_delayTime) then


    !if (mod(m,sw_delayTime) == 0) then

        !Ni_tot_1 = Ni_tot+1
       ! Ni_tot_2 = Ni_tot + Nx_boundary + FS_boundary
       ! Ni_tot = Ni_tot_2

       ! call load_foreshock_Maxwellian(vth,Ni_tot_1,Ni_tot_2,mion,1.0,1,swfsRatio) !SW
       if (my_rank .eq. 0) then
       		!write(*,*) '(m-1*Nx*sw_delayTime)/sw_delayTime,m-1*Nx*sw_delayTime, sw_delayTime, FSstartupTime...',(m-1*Nx*sw_delayTime)/sw_delayTime, m-1*Nx*sw_delayTime,sw_delayTime,FSstartupTime
       		!!write(*,*) 'Ni_tot, SolarWind Ions, Foreshock Ions,Diff,Ratio...', Ni_tot,Nx_boundary,FS_boundary,Ni_tot_2-Ni_tot_1,swfsRatio
       endif
    !endif
           !if (my_rank .eq. 0) then
       		!write(*,*) '(m-1*Nx*sw_delayTime)/sw_delayTime,m-1*Nx*sw_delayTime, sw_delayTime, FSstartupTime...',(m-1*Nx*sw_delayTime)/sw_delayTime, m-1*Nx*sw_delayTime,sw_delayTime,FSstartupTime
       !endif
    
if (m .ge. 1*Nx*sw_delayTime) then !wait for peaceful SW 
    !if (mod(m,sw_delayTime) == 0) then  
       !if ((m-1*Nx*sw_delayTime)/sw_delayTime .lt. FSstartupTime) then
       !Ni_tot_1 = Ni_tot+1
       !Ni_tot_2 = Ni_tot + (m-1*Nx*sw_delayTime)/sw_delayTime
       !Ni_tot = Ni_tot_2
       !else
!       Ni_tot_1 = Ni_tot+1
!       Ni_tot_2 = Ni_tot + FS_boundary
!       Ni_tot = Ni_tot_2
       !endif

 !      call load_foreshock_Maxwellian(vth,Ni_tot_1,Ni_tot_2,mion,1.0,2)  !Foreshock
       
    !endif
endif








endif









!write(*,*) 'before subcycle,b1,b0', b1(1,2,51,2),b0(1,2,51,2)

            call get_interp_weights()

            call update_np()                  !np at n+1/2

            call update_up(vp)            !up at n+1/2
            
            call get_gradP()
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Subcycling loop

            dtsub = dtsub_init
            ntf=ntsub
            call check_time_step(bt,np,dtsub,ntf)
            
            do n=1,ntf
                  call curlB(b0,bt,np,aj)

                  
                  call predict_B(b0,b12,b1p2,bt,E,aj,up,nu,dtsub)
                  

                  
                  call correct_B(b0,b1,b1p2,E,aj,up,np,nu,dtsub)

                  
                  call f_update_tlev(b1,b12,b1p2,bt,b0)

                  
            enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           
 
            call move_ion_half()       !final ion move to n+1
    




 !call fieldMove_Boundary(b1,m)

!save velocity of particles at SC position
 
 !call get_v_sc(150,2,250)
                   

 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       diagnositc output
            
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            
            if (my_rank .eq. 0) then
                  write(160) m
                  write(160) input_E,input_EeP,Evp,Euf,EB1,EB1x,EB1y,EB1z,EE, &
                        EeP,input_chex,input_bill
                  write(190) m
                  write(190) pup,puf,peb,input_p
                   !write(*,*) 'n...', nSC
                   !write(*,*) 'x...', xSC
                   !write(*,*) 'v...', vSC
                   
                  !fft output
                  !write(401) b1(2,2,11400,1), b1(2,2,11400,2), b1(2,2,11400,3)
                  !write(402) b1(2,2,14000,1), b1(2,2,14000,2), b1(2,2,14000,3)
                  !write(403) b1(2,2,17000,1), b1(2,2,17000,2), b1(2,2,17000,3)
                  
            endif
            
            ndiag = ndiag+1
            ndiag_part = ndiag_part + 1
            if (ndiag .eq. nout) then
                  call get_temperature()
                  call update_rho()
                  !call get_v_dist()
                  call update_mixed()
                  !write(*,*) 'outside update_np_cold', Ni_tot
                  !Update solar wind ion densities
		  call update_np_cold()                
		  !Update solar wind ion temperatures
		  call get_temperature_cold()               
		  !Update foreshock ion densities
		  call update_np_mixed()                
		  !Update foreshock ion temperatures
		  call get_temperature_mixed()
		  !Update solar wind flow
		  call update_up_cold(vp)
		  !Update foreshock flow
		  call update_up_mixed(vp)

 
                  if (my_rank .eq. 0) then
                  	call face_to_center(E,Ec)
                        write(110) m
                        write(110) np
                        write(115) m
!                        write(115) np_b
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
                        !write(305) testPartPos!ionPos(1:5,:)!xp(1:100,:)
                        !write(310) m
                        !write(310) vp
                        write(315) m
                        write(315) gradP
                        write(317) m
                        write(317) ue
                        !write(330) m
                        !write(330) vSC
                        write(340) m
                        write(340) testPartPos
                        write(350) m
                        write(350) mnp

                        
!                        write(342) m
!                        write(342) vp(299985:Ni_max,:)
                        
                        
                   endif
                   ndiag = 0
              endif
                   !Output Dist data
              if (ndiag_part .eq. nout*10) then
                   if (my_rank .ge. 0) then
                        write(120) m
                        write(120) mix_ind
                        write(305) m
                        write(305) xp
                        write(310) m
                        write(310) vp
                       ! write(315) m
                       ! write(315) nSC
                       ! write(317) m
                       ! write(317) xSC
                       ! write(330) m
                       ! write(330) vSC
                    endif    
                        !ndiag_part = 0  !Only print the first one. 
            endif
!            write(*,*) 'minimum density.....', minval(np(:,:,:))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Write restart file

            
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
      enddo     !End Main Loop
      
!      if (my_rank .eq. 0) then
            if (restart) then
                  if (my_rank .eq. 0) then
                        write(*,*) 'Writing restart file...'
                        
                        !Write grid data
                        
                        open(220,file=trim(out_dir)//'restart.vars',status='unknown',form='unformatted')
                        write(220) b1,b12,b1p2,bt,btmf,btc,np,np3, &
                        up,aj,nu,E,temp_p,mnp,beta,Evp,Euf, &
                        EB1,EB1x,EB1y,EB1z,EE,EeP, &
                        input_EeP,input_Eb,prev_Etot,bndry_Eflux,grav, &
                        input_chex,input_bill
                  endif
                              
                              
                  !write individual processor data
                  open(212,file=trim(out_dir)//'restart.part'//trim(filenum),status='unknown',form='unformatted')
                  write(212) vp,vp1,vplus,vminus,xp,Ep,input_E,Ni_tot,ijkp,input_p,mrat,m_arr,beta_p
                  
                  close(220)
                  close(212)
            endif
!      endif
      
      close(110)
      close(115)
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

      
