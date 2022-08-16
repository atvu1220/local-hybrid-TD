module grid
      use dimensions
      implicit none
      save
      real:: qx(nx), qy(ny), qz(nz) !cartesian coordinates of grid points
      real:: lambda(nz)             !latitude of given grid point
      real:: dz_grid(nz), dz_cell(nz) !dz for grid/cell, variable in z
      real:: dx_grid(nx), dx_cell(nx), dy_grid(ny), dy_cell(ny)
      real:: xrat(nx), yrat(ny), zrat(nz)
      integer:: ri,rj,rk
      
end module grid
      