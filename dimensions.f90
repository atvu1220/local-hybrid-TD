module dimensions
      implicit none
      save
      integer, parameter:: nx = 202, ny = 3, nz = 400 !1001,1001 !Normal simulations, more normal
      !integer, parameter:: nx = 302, ny = 3, nz = 600 !1001,1001 !Normal simulations
      
     ! integer, parameter:: nx = 1002, ny = 3, nz = 900 !1001,1001 !Bow shock reflected, big 2022
      !integer, parameter:: nx = 600, ny = 3, nz = 150 !1001,1001 !Bow shock reflected
      !integer, parameter:: nx = 202, ny = 3, nz = 600 !1001,1001
      !integer, parameter:: nx = 802, ny = 3, nz = 600 !1001,1001 ! zero shear, reflect, 2000000 particles
      !integer, parameter:: nx = 602, ny = 3, nz = 900 !1001,1001 !15 degreee shear, reflect, 3000000 particles
      !!integer, parameter:: nx = 1002, ny = 3, nz = 900 !1001,1001 !15 degreee shear, reflect, 3000000 particles
      !!integer*4, parameter:: Ni_max =8000000!2500000!500000  
      integer*4, parameter:: Ni_max =2000000!2500000!500000  
     
      integer, parameter:: nTestParticles_max = 1000 !Number of FS Test Particles
      
end module dimensions
      
