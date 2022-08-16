module grid_interp
      implicit none
      contains

      subroutine edge_to_center(bt,btc)
            use dimensions
            use boundary
            use grid, only: xrat,yrat,zrat
            implicit none
            real, intent(inout):: bt(nx,ny,nz,3)
            real, intent(out):: btc(nx,ny,nz,3)
            real:: b1,b2, btmf(nx,ny,nz,3)
            integer:: i,j,k,im,jm,km, ip, jp, kp
            
            call boundary_vector(bt)
!            call periodic(bt)
            
            do i=2,nx-1
                  do j=2,ny-1
                        do k=2,nz-1
                              im=i-1
                              jm=j-1
                              km=k-1
                              
                              b2 = bt(i,jm,k,1) + yrat(j)*(bt(i,j,k,1)-bt(i,jm,k,1))
                              b1 = bt(i,jm,km,1)+ yrat(j)*(bt(i,j,km,1)-bt(i,jm,km,1))
                              
                              btc(i,j,k,1) = b1 + zrat(k)*(b2-b1)
                              
                              b2 = bt(im,j,k,2) + xrat(i)*(bt(i,j,k,2)-bt(im,j,k,2))
                              b1 = bt(im,j,km,2)+ xrat(i)*(bt(i,j,km,2)-bt(im,j,km,2))
                              
                              btc(i,j,k,2) = b1 + zrat(k)*(b2-b1)
                              
                              b2 = bt(i,jm,k,3) + yrat(j)*(bt(i,j,k,3)-bt(i,jm,k,3))
                              b1 = bt(im,jm,k,3)+ yrat(j)*(bt(im,j,k,3)-bt(im,jm,k,3))
                              
                              btc(i,j,k,3) = b1 + xrat(i)*(b2-b1)
                              
                        enddo
                  enddo
            enddo
          
          do i=2,nx-1
                  do j=2,ny-1
                        do k=2,nz-1
                              ip=i+1
                              jp=j+1
                              kp=k+1
                              
                              btmf(i,j,k,1) = 0.5* ( btc(i,j,k,1) + btc(ip,j,k,1) )
                              btmf(i,j,k,2) = 0.5* ( btc(i,j,k,2) + btc(i,jp,k,2) )  
                              btmf(i,j,k,3) = 0.5* ( btc(i,j,k,3) + btc(i,j,kp,3) )    
                              
                        enddo
                  enddo
            enddo
          
          
          
            
            call edge_to_face(bt,btmf)
            call face_to_center(btmf,btc)
            
            call boundary_vector(btc)
!            call periodic(btc)
            
      end subroutine edge_to_center
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine edge_to_face(bt,btmf)
            use dimensions
            use boundary
            use grid, only: xrat,yrat,zrat
            implicit none
            real, intent(inout):: bt(nx,ny,nz,3)
            real, intent(out):: btmf(nx,ny,nz,3)
            real:: btc(nx,ny,nz,3),b1,b2
            integer:: i,j,k,im,jm,km
            
            call boundary_vector(bt)
!            call periodic(bt)
            
            do i=2,nx
                  do j =2,ny
                        do k=2,nz
                              im=i-1
                              jm=j-1
                              km=k-1
                              
                              b2 = bt(i,jm,k,1) + yrat(j)*(bt(i,j,k,1)-bt(i,jm,k,1))
                              b1 = bt(i,jm,km,1)+ yrat(j)*(bt(i,j,km,1)-bt(i,jm,km,1))
                              
                              btc(i,j,k,1) = b1 + zrat(k)*(b2-b1)
                              
                              b2 = bt(im,j,k,2) + xrat(i)*(bt(i,j,k,2)-bt(im,j,k,2))
                              b1 = bt(im,j,km,2)+ xrat(i)*(bt(i,j,km,2)-bt(im,j,km,2))
                              
                              btc(i,j,k,2) = b1 + zrat(k)*(b2-b1)
                              
                              b2 = bt(i,jm,k,3) + yrat(j)*(bt(i,j,k,3)-bt(i,jm,k,3))
                              b1 = bt(im,jm,k,3)+ yrat(j)*(bt(im,j,k,3)-bt(im,jm,k,3))
                              
                              btc(i,j,k,3) = b1 + xrat(i)*(b2-b1)
                              
                        enddo
                  enddo
            enddo
            
            call boundary_vector(btc)
!            call periodic(btc)
            
            do i=2,nx-1
                  do j=2,ny-1
                        do k=2,nz-1
                              btmf(i,j,k,1) = 0.5*(btc(i,j,k,1)+btc(i+1,j,k,1))
                              btmf(i,j,k,2) = 0.5*(btc(i,j,k,2)+btc(i,j+1,k,2))
                              btmf(i,j,k,3) = 0.5*(btc(i,j,k,3)+btc(i,j,k+1,3))
                        enddo
                  enddo
            enddo
!            call boundary_vector(btmf)            
!            call periodic(btmf)
            
      end subroutine edge_to_face
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine face_to_center(v,vc)
            use dimensions
            use boundary
            use grid, only: xrat,yrat,zrat
            implicit none
            real, intent(inout):: v(nx,ny,nz,3)
            real, intent(out):: vc(nx,ny,nz,3)
            integer:: i,j,k,im,jm,km
            
            call boundary_vector(v)
!            call periodic(v)
            
            do i=2,nx
                  do j=2,ny
                        do k=2,nz
                              im=i-1
                              jm=j-1
                              km=k-1
                              
                              vc(i,j,k,1) = xrat(i)*(v(i,j,k,1)-v(im,j,k,1)) + v(im,j,k,1)
                              vc(i,j,k,2) = yrat(j)*(v(i,j,k,2)-v(i,jm,k,2)) + v(i,jm,k,2)
                              vc(i,j,k,3) = zrat(k)*(v(i,j,k,3)-v(i,j,km,3)) + v(i,j,km,3)
                              
                        enddo
                  enddo
            enddo
            
            call boundary_vector(vc)
!            call periodic(vc)
            
      end subroutine face_to_center
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine grav_to_center(grav,gravc)
            use dimensions
            use boundary
            use grid, only: zrat
            implicit none
            real, intent(inout):: grav(nx,ny,nz)
            real, intent(out):: gravc(nx,ny,nz)
            integer:: i,j,k,km
            
            call boundary_scalar(grav)
!            call periodic_scalar(grav)
            do i=2,nx
                  do j=2,ny
                        do k=2,nz
                              
                              km=k-1
                              
                              gravc(i,j,k) = zrat(k)*(grav(i,j,k)-grav(i,j,km)) + grav(i,j,km)
                        enddo
                  enddo
            enddo
            
            call boundary_scalar(gravc)
!            call periodic_scalar(gravc)
            
      end subroutine grav_to_center
      
end module grid_interp

            
