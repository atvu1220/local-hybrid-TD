module Var_Arrays
      use dimensions
      implicit none
      real::      b0(nx,ny,nz,3), &     !ambient mag field
                  b1(nx,ny,nz,3), &     !1st order mag field   
                  b12(nx,ny,nz,3), &    !b1 at previous time step
                  b1p2(nx,ny,nz,3), &   !temp b1 at time level m+1
                  bt(nx,ny,nz,3), &     !total mag field, mc covarient
                  btmf(nx,ny,nz,3), &   !main cell contravarient bt field
                  btc(nx,ny,nz,3), &    !btmf at cell center for particle move
                  np(nx,ny,nz), &       !particle ion density at level n, n+1/2
                  np3(nx,ny,nz,3), &    
                  vp(Ni_max,3), &       !particle velocity at t level n+1/2
                  vp1(Ni_max,3), &      !particle velocity at t level n
                  vplus(Ni_max,3), &    !v+ used in velocity update
                  vminus(Ni_max,3), &   !v- used in velocity update
                  up(nx,ny,nz,3), &     !particle flow at n, n+1/2
                  xp(Ni_max,3), &       !coordinates of ion particles
                  aj(nx,ny,nz,3), &     !curlB/(alpha*n)
                  nu(nx,ny,nz), &       !collision frequency
                  Ep(Ni_max,3), &       !Ion particle electric field
                  E(nx,ny,nz,3), &      !E field from electron mom. eq.
                  Ec(nx,ny,nz,3),&      !E field at center of cell
                  ExB(nx,ny,nz,3),&     !ExB Drift Velocity from Initial Solar Wind Setup
                  temp_p(nx,ny,nz), &   !temperature
                  mnp(nx,ny,nz), &      !mass density
                  beta, beta_p(Ni_max), &  !variable for particle scaling
                  m_arr(Ni_max), &      !Mass of macroparticles
                  np_t(nx,ny,nz), &     !Ion density of Top
                  np_b(nx,ny,nz), &     !Ion Density of Bot
                  up_t(nx,ny,nz,3), &   !Ion bulk flow of Top
                  up_b(nx,nz,nz,3), &   !Ion bulk flow of Bot
                  input_p(3), &
                  input_E, input_Eb, prev_Etot, &
                  gradP(nx,ny,nz,3),&            !electron pressure gradient
                  mixed(nx,ny,nz),&     !0/1 Marker for either solar wind or foreshock
                  Ptherm(nx,ny,nz), &
                  PB(nx,ny,nz), &
                  Ptotal(nx,ny,nz), &
                  np_cold(nx,ny,nz), & 
                  temp_p_cold(nx,ny,nz), &
                  np_cold_foreshock(nx,ny,nz), & 
                  temp_p_cold_foreshock(nx,ny,nz), &
                  np_mixed(nx,ny,nz), & 
                  temp_p_mixed(nx,ny,nz), &
                  np_mixed_foreshock(nx,ny,nz), & 
                  temp_p_mixed_foreshock(nx,ny,nz), &                 
                  up_cold(nx,ny,nz,3), &
                  up_mixed(nx,ny,nz,3), &
                  up_cold_foreshock(nx,ny,nz,3), &
                  up_mixed_foreshock(nx,ny,nz,3), &
                  tp_cold(nx,ny,nz,3), &
                  tp_mixed(nx,ny,nz,3), &
                  tp_cold_foreshock(nx,ny,nz,3), &
                  tp_mixed_foreshock(nx,ny,nz,3), &
                  tp(nx,ny,nz,3),&
                  curlBcurrent(nx,ny,nz,3), & !mu0*J current.
                  ue(nx,ny,nz,3)
      
      integer(4):: Ni_tot, Ni_tot_sys, Ni_init
      
      !Location (indices) of particles in the grid
      integer:: ijkp(Ni_max,3)
      logical:: in_bounds(Ni_max)
      real:: mix_ind(Ni_max)
      
      !Weight variables for trilinear interpolation
      real:: wght(Ni_max,8)
      
end module Var_Arrays