program hybrid

      use mpi
      use dimensions
      use Var_Arrays
      use mult_proc, only: my_rank, procnum
      use inputs
      use misc
      
      implicit none
      
      real:: Evp,Euf,EB1,EB1x,EB1y,EB1z,EE,input_EeP,prev_Etot
      real:: input_chex, input_bill
      real:: pup(3),puf(3),peb(3)
      character(2):: filenum(16) !max 16 processors
      integer:: ierr,t1,t2,cnt_rt,m,mstart,ndiag,seed
      real(8):: time
      logical:: restart = .false.
      integer(4):: Ni_tot_sw,Ni_tot_sys
      integer:: i,j,k !looping indicies
      
      filenum = (/'1 ','2 ','3 ','4 ','5 ','6 ','7 ','8 ','9 ', &&
            '10','11','12','13','14','15','16'/)
            
      call readInputs()
      call initparameters()
      
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, procnum, ierr)
      
      call system_clock(t1, cnt_rt)
      
!Initialize all variables

      Ni_tot=nz*(ppc/procnum) !1D
      Ni_tot_0 = Ni_tot
      Ni_tot_sw = Ni_tot
      Ni_tot_sys = Ni_tot*procnum
      
      if (my_rank .eq. 0) then
            call check_inputs()
            write(*,*) 'Partilces per cell... ', Ni_tot_sys/nz
            write(*,*) ' '
      endif
      
      mstart = 0
      ndiag = 0
      prev_Etot = 1.0
!      nuei = 0.0
      
      seed = t1 + my_rank*100
      call random_initialize(seed)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      
      if (.not. restart) then
            do i=1,nx
                  do j=1,ny
                        do k=1,nz
                              input_E = 0.0
                              input_p = 0.0
                              input_chex = 0.0
                              input_bill = 0.0
                        enddo
                  enddo
            enddo
      endif
      
      call grd7()
      call grd6_setup(b0,bt,b12,b1,b1p2,nu,input_Eb)
      call get_beta(Ni_tot_sys,beta)
      
      input_E = 0.0
      
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      
!       Initialize particles
!      call sw_part_setup_maxwl_mix(np,vp,vp1,xp,input_p,up,np_t_flg,np_b_flg)
      
      call load_Maxwellian !Add variables here
!      call load_aniso_Maxwellian

      Ni_tot_sys = Ni_tot*procnum
      write(*,*) Ni_tot_sys, Ni_tot, procnum
      write(*,*) 'Particles per cell...' Ni_tot_sys/nz
      
      call f_update_tlev(b1,b12,b1p2,bt,b0)
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!  Check for restart flag

      write(*,*) 'restart status...', restart
      if (restart) then
            write(*,*) 'opening restar vars....'
            open(210,file='restart.vars',status='unknown',form='unformatted')
            write(*,*) 'reading restart vars......'
            read(210) b1,b12,b1p2,bt,btmp,nn,n,nf,vp,vp1,vplus,vminus, &
                  up,xp,uf,uf2,ufp2,aj,Ep,Ef,E,uplus,uminus,Evp,Euf, &
                  EB1,EB1x,EB1y,EB1z,EE,EeP,input_E,Ni_tot, &
                  ijkp,mstart,input_p,input_EeP,prev_Etot,nf1,nf3,nfp1, &
                  input_chex,input_bill,pf,pf1,mrat,m_arr
            write(*,*) 'restarting hybrid ....'
            
            if (my_rank .gt. 0) then
                  open(211,file='restart.part'//trim(filenum(my_rank)),status='unknown',form='unformatted')
                  read(211) vp,vp1,vplus,vminus,xp,Ep,input_E,Ni_tot,ijkp,input_p,mrat,m_arr
            endif
      endif
      
      close(211)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!       Inititalize diagnositc output files

      if (my_rank .eq. 0) then
            open(110,file=trim(out_dir)//'c.np.dat',status='unknown',form='unformatted')
            open(115,file=trim(out_dir)//'c.np_b.dat',status='unknown',form='unformatted')
            open(120,file=trim(out_dir)//'c.mixed.dat',status='unknown',form='unformatted')
            open(130,file=trim(out_dir)//'c.b1.dat',status='unknown',form='unformatted')
            open(140,file=trim(out_dir)//'c.aj.dat',status='unknown',form='unformatted')
            open(150,file=trim(out_dir)//'c.E.dat',status='unknown',form='unformatted')
            open(160,file=trim(out_dir)//'c.energy.dat',status='unknown',form='unformatted')
            
            !diagnostics chex,bill,satnp
            
            open(180,file=trim(out_dir)//'c.up.dat',status='unknown',form='unformatted')
            open(190,file=trim(out_dir)//'c.momentum.dat',status='unknown',form='unformatted')
            open(192,file=trim(out_dir)//'c.p_conserve.dat',status='unknown',form='unformatted')
            open(300,file=trim(out_dir)//'c.temp_p.dat',status='unknown',form='unformatted')
            open(305,file=trim(out_dir)//'c.xp_0.dat',status='unknown',form='unformatted')
            open(310,file=trim(out_dir)//'c.vp_0.dat',status='unknown',form='unformatted')
            open(315,file=trim(out_dir)//'c.mrat_0.dat',status='unknown',form='unformatted')
            open(317,file=trim(out_dir)//'c.beta_p_0.dat',status='unknown',form='unformatted')
            open(320,file=trim(out_dir)//'c.np_wake.dat',status='unknown',form='unformatted')
            open(330,file=trim(out_dir)//'c.up_t.dat',status='unknown',form='unformatted')
            open(340,file=trim(out_dir)//'c.up_b.dat',status='unknown',form='unformatted')
            open(342,file=trim(out_dir)//'c.test_part.dat',status='unknown',form='unformatted')
            open(350,file=trim(out_dir)//'c.mnp.dat',status='unknown',form='unformatted')
       endif
       
       if (my_rank .gt. 0) then
            open(305,file=trim(out_dir)//'c.xp_'//trim(filenum(my_rank))//'.dat',status='unknown',form='unformatted')
            open(310,file=trim(out_dir)//'c.vp_'//trim(filenum(my_rank))//'.dat',status='unknown',form='unformatted')
            open(315,file=trim(out_dir)//'c.mrat_'//trim(filenum(my_rank))//'.dat',status='unknown',form='unformatted')
            open(317,file=trim(out_dir)//'c.beta_p__'//trim(filenum(my_rank))//'.dat',status='unknown',form='unformatted')
       endif
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!       MAIN LOOP

      do m = mstart+1, nt
            if (my_rank .eq. 0) then
                  write(*,*) 'time...' m, m*dt,my_rank
            endif
            
            if (Ni_tot .lt. Ni_max) then
                  call Mass_load_Io(np,vp,vp1,xp,m,input_p,up)  !call ionizing subroutine (ionize_pluto, etc.)
            endif
            
            call get_interp_weights(xp)
            call update_np(np)                  !np at n+1/2
            call update_up(vp,np,up)            !up at n+1/2
            
            !energy diagnostics
            call get_bndry_Eflux(b1,E)
            call Energy_diag(vp,b0,b1,E,Evp.Euf,EB1,EB1x,EB1y,EB1z,EE,EeP,nu,up,np)
            
            call curlB(bt,np,aj)
            call edge_to_center(bt,btc)
            call extrapol_up(up,vp,vp1,np)
            call get_Ep(Ep,aj,np,up,btc,nu)
            call get_vplus_vminus(Ep,btc,vp,vplus,vminus)
            call improve_up(vp1,vplus,vminus,up,np)
            
            call get_Ep(Ep,aj,np,up,btc,nu)
            call get_vplus_vminus(Ep,btc,vp,vplus,vminus)
            call get_vp_final(Ep,vp,vp1,vplus)
            
            call move_ion_half(xp,vp,vp1,input_p) !1/2 step ion move to n+1/2
            
            call get_interp_weights(xp)
            call update_np(np)                   !np at n+1/2
            call update_up(vp,np,up)             !up at n+1/2
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Subcycling loop

            dtsub = dtsub_init
            ntf=ntsub
            
            call check_time_step(bt,np)
            
            do n=1,ntf
                  call curlB(bt,np,aj)
                  
                  call predict_B(b0,b1,b12,b1p2,bt,E,aj,up,np,nu)
                  
                  call correct_B(b0,b1,b1p2,E,aj,up,np,nu)
                  
                  call f_update_tlev(b1,b12,b1p2,bt,b0)
                  
            enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1           
 
            call move_ion_half(xp,vp,vp1,input_p)       !final ion move to n+1
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       diagnositc output
            
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            
            if (my_rank .eq. 0) then
                  write(160) m
                  write(160) input_E,input_EeP,Evp,Euf,EB1,EB1x,EB1y,EB1z,EE, &
                        EeP,input_chex,input_bill
                  write(190) m
                  write(190) pup,puf,peb,input_p
                  
            endif
            
            ndiag = ndiag+1
            if (ndiag .eq. nout) then
                  call get_temperature(xp,vp,np,temp_p)
                  call update_rho(mnp)
                  call update_mixed()
                  if (my_rank .eq. 0) then
                        write(110) m
                        write(110) np
                        write(115) m
                        write(115) np_b
                        write(120) m
                        write(120) mixed
                        write(130) m
                        write(130) bt
                        write(140) m
                        write(140) aj*alpha*np3
                        write(150) m
                        write(150) E
                        write(180) m
                        write(180) up
                        write(300) m
                        write(300) temp_p/1.6e-19
                        write(305) m
                        write(305) xp
                        write(310) m
                        write(310) vp
                        write(315) m
                        write(315) mrat
                        write(317) m
                        write(317) beta_p
                        write(330) m
                        write(330) up_t
                        write(340) m
                        write(340) up_b
                        write(350) m
                        write(350) mn
                        
                        ndiag = 0
                   endif
                   
                   if (my_rank .gt. 0) then
                        write(305) m
                        write(305) xp
                        write(310) m
                        write(310) vp
                        write(315) m
                        write(315) mrat
                        write(317) m
                        write(317) beta_p
                        
                        ndiag = 0
                  endif
                   
            endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Write restart file

            if (my_rank .eq. 0) then
                  if (m .eq. mrestart) then
                        open(220),file='restart.vars.new',status='unknown',form='unformatted')
                        write(220) b1,b12,b1p2,bt,btmp,nn,n,nf,vp,vp1,vplus,vminus, &
                              up,xp,uf,uf2,ufp2,aj,Ep,Ef,E,uplus,uminus,Evp,Euf, &
                              EB1,EB1x,EB1y,EB1z,EE,EeP,input_E,Ni_tot, &
                              ijkp,mstart,input_p,input_EeP,prev_Etot,nf1,nf3,nfp1, &
                              input_chex,input_bill,pf,pf1,mrat,m_arr
                              
                  endif
            endif
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            
      enddo             !End Main Loop
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
      close(180)
      close(190)
      close(192)
      close(210)
      close(220)
      close(221)
      close(300)
      close(305)
      close(310)
      close(315)
      close(317)
      close(320)
      close(330)
      close(340)
      close(342)
      close(350)
      
      call system_clock(t2,cnt_rt)
      time=(real(t2) - real(t1))/real(cnt_rt)
      if (my_rank .eq. 0) then
            write(*,*)
            write(*,*) 'Elapsed time .....', time, ' sec'
            write(*,*)
      endif
      
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      
      call MPI_FINALIZE(ierr)
      
end program hybrid
      
