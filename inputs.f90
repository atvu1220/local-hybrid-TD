module inputs
      use dimensions
      use mpi
      implicit none
      save
      
      real:: b0_init, nf_init,dt_frac, vsw, vth, plasma_beta, dx_frac, &
            nu_init_frac,lambda_i,mO, ppc, nu_init, &
            etemp0, mion, va, va_f, &
            FSDriftSpeed, FSDensityRatio, ForeshockBeta, &
            TDcellBalance, magneticShear, &
            ByConeAngle, ByConeAngleDelta, FSThermalRatio
            
      real, parameter:: amu=1.6605e-27!, mion = 3.841e-26
      integer:: mp, nt, nout, boundx,ddthickness,FSBeamWidth, quasiparallel
      integer:: Ni_tot_0

      real, parameter:: q=1.6e-19         !electron charge

!       Grid paramters
      
      real:: dx,dy,delz,dt,dtsub_init

!       time stepping paramters

      integer, parameter:: ntsub = 10    !number of subcycle timesteps
      
!       Output directory

      character(50):: out_dir
      
!       Neutral cloud expansion characteristics
      real:: vtop, vbottom
      
!       Max number of ions to be produced
!      integer(4):: Ni_tot_0
      
!       Misc. constants

      real, parameter:: pi = 3.14159
      real, parameter:: mu0 = pi*4.0e-7, &
                        epsilon = 8.85e-12, &
                        rtod=180.0/pi, &
                        km_to_m = 1.0e3, &
                        kboltz = 1.38e-29, & !This is in units of km already...1.38e-23 * 10^-9 km^-3 * 1000
                        tempf0 = 50*11600.0
          
      real, parameter:: beta_particle = 1.0     !beta value of intial particles       
      real, parameter:: beta_pu = 10.0            !beta value of pickup ions
      real:: omega_p                            !ion gyrofrequency
      
!       Electron ion collision frequency
      real, parameter:: lww2 = 0.985!0.995!0.995!0.985!0.995!0.999!0.999!9!1.0!0.9990!0.9995!0.9995!0.9994!0.995          !artificial diffusion for the magnetic field update
      real, parameter:: lww1 = (1-lww2)/6.0     !divide by six for nearest neighbor
      
!       Density scaling paramter, alpha, and ion particle array dims

      real:: alpha
      
      
      contains
      
            subroutine readInputs()
                 implicit none
                  
                 open(unit=100, file= 'inputs.dat', status='old')
                  
                 read(100,*) b0_init
                 write(*,*) 'b0_init...........',b0_init
                 read(100,*) nf_init
                 write(*,*) 'nf_init...........',nf_init
                 read(100,*) dt_frac
                 write(*,*) 'dt_frac...........',dt_frac
                 read(100,*) nt
                 write(*,*) 'nt................',nt
                 read(100,*) nout
                 write(*,*) 'nout..............',nout
                 read(100,*) va_f
                 write(*,*) 'va_f...............',va_f
                 read(100,*) plasma_beta
                 write(*,*) 'plasma_beta...............',plasma_beta
                 read(100,*) dx_frac
                 write(*,*) 'dx_frac...........',dx_frac
                 read(100,*) nu_init_frac
                 write(*,*) 'nu_init_frac......',nu_init_frac
                 read(100,*) ppc
                 write(*,*) 'part per cell.....',ppc
                 read(100,*) etemp0
                 write(*,*) 'electon temperature (eV)...', etemp0
                 read(100,*) boundx
                 write(*,*) 'boundary condition......', boundx
                 read(100,*) ddthickness
                 write(*,*) 'discontinuity thickness........',ddthickness
                 read(100,*) FSBeamWidth
                 write(*,*) 'Foreshock Beam Width........',FSBeamWidth
                 read(100,*) FSDensityRatio
                 write(*,*) 'Foreshock Density Ratio........',FSDensityRatio
                 read(100,*) FSThermalRatio
                 write(*,*) 'Foreshock Thermal Ratio........',FSThermalRatio
                 read(100,*) FSDriftSpeed
                 write(*,*) 'Foreshock Drift Speed........',FSDriftSpeed
                 read(100,*) magneticShear
                 write(*,*) 'Magnetic Shear........',magneticShear
                 read(100,*) quasiparallel
                 write(*,*) 'Number of quasiparallel sides........',quasiparallel
                 read(100,*) ByConeAngle
                 write(*,*) 'By Cone Angle........', ByConeAngle
                 read(100,*) ByConeAngleDelta
                 write(*,*) 'By Cone Angle Delta........', ByConeAngleDelta
                 read(100,*) out_dir
                 write(*,*) 'output dir........',out_dir
                 
                 close(100)
                 
                
            end subroutine readInputs
            
            
            
            subroutine initparameters()
                  implicit none
                  mion = amu*1.0!3.841e-26
                  write(*,*) 'mion...',mion

                !vth is input as the plasma beta
                  vth = sqrt(plasma_beta*B0_init**2/(mu0*mion*nf_init/1e9))/1e3
                  write(*,*) 'vth...',vth
                  
                  omega_p = q*b0_init/mion
                  vsw=va_f*vth
                  lambda_i = (3e8/sqrt((nf_init/1e9)*q*q/(8.85e-12*mion)))/1e3
                  write(*,*) 'lambda_i',lambda_i    
                  dx= lambda_i*dx_frac
                  dy= lambda_i*dx_frac           !units in km
                  delz = lambda_i*dx_frac       !dz at release coordinates
                  ddthickness = int(float(ddthickness)/dx_frac) !Ddthickness in the number of cells...ion inertial lengths.
                  dt= dt_frac*mion/(q*b0_init)  !main time step
                  dtsub_init = dt/ntsub         !subcycle time step for magnetic field updates
                  vtop = vsw
                  vbottom = -vsw
                  
                  Ni_tot_0 = int(Ni_max*1.0)
                  write(*,*) 'Ni_tot_0...',Ni_tot_0, Ni_max
                  
                  ForeshockBeta = 50!20.0!0.5; !Weighting of Foreshock Ions. 0.5% = 1/200, 1 Foreshock Ion to 1/200 of Solar Wind Ions, FSBeta = 20  means that one foreshock particle is 5%
                  TDcellBalance = -0.25*float(ddthickness)*(-1.0 + log(cosh(1.0))) !Calculates how many additional cells of particles to fill in TD for pressure balance 0-a
                  write(*,*) 'TDcellBalance...',TDcellBalance
                  write(*,*) 'Foreshock Beta...',ForeshockBeta
                  
                  mO = mion     !mass of H (kg)
                  
                  nu_init = nu_init_frac*omega_p
                  
                  alpha = (mu0/1.0e3)*q*(q/mion)  !mH...determines particle scaling
                  write(*,*) 'alpha...', alpha
                  
            end subroutine initparameters
            
            
            subroutine check_inputs()
                  use mult_proc, only: my_rank
                  use grid, only: dz_grid
                  use var_arrays, only: np
                  implicit none
                  integer:: i,j,k
                  
                  real:: ak, btot, a1, a2, womega, phi, deltat, cwpi
                  
      ! Check input paramters
      
                  if (my_rank .eq. 0) then
                        write(*,*) mu0,q,mion
                        write(*,*) 'alpha...', alpha
                        write(*,*) 'c/wpi...', lambda_i,dx,dy,delz
                        write(*,*) 'dt......', dt
                        write(*,*) 'dt_sub...', dtsub_init
                        
                        
                        do i=1,nx
                              do j=1,ny
                                    do k=1,nz
                        
                                          ak = 2.0/dz_grid(k)
                                          btot = b0_init*q/mion
                                          a1 = ak**2*btot/(alpha*np(i,j,k))
                                          a2 = (ak*btot)**2/(alpha*np(i,j,k))
                                          womega = 0.5*(a1+ sqrt(a1**2 + 4.0*a2))
                                          phi = womega/ak
                                          deltat = dz_grid(k)/phi
                                          if (deltat/dtsub_init .le. 2.0) then
                                                write(*,*) 'deltat/dtsub....', deltat/dtsub_init
                                                write(*,*) 'Field time stp too long...', i,j,k
                                                write(*,*) 'ak,btot...', ak,btot
                                                write(*,*) 'womega,b0_init', womega,b0_init
                                                write(*,*) 'phi,dz_grid(k)...', phi,dz_grid(k)
                                                write(*,*) 'np(i,j,k)', np(1,1,k)
                                                
                                                stop
                                          endif
                                    enddo
                              enddo
                        enddo
                        
                        write(*,*) 'Courant check (>1)...', deltat/dtsub_init
                        write(*,*) 'Max Number of Particles...', Ni_max
      
                        write(*,*) '  '
                        write(*,*) '  '
                        
                        va = b0_init/sqrt(mu0*mion*nf_init/1e9)/1e3
                        
                        write(*,*) 'Alfven velocity......', va
                        write(*,*) 'Thermal velocity.....', vth
                        write(*,*) 'Mach number..........', va/(va)
                        
                        write(*,*) 'Thermal gyroradius...', mion*vth/(q*b0_init), mion*vth/(q*b0_init)/dx
                        cwpi = 3.0e8/sqrt((nf_init/1e9)*q*q/(epsilon*mion))
                        write(*,*) 'Ion inertial length..', cwpi/1e3,cwpi/1e3/dx
                        write(*,*) '  '
                        
                  endif
                  
                  
            end subroutine check_inputs
            
end module inputs
                        
