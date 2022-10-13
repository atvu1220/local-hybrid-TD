module initial
use dimensions
implicit none
save
contains
      
subroutine grd6_setup(b0,bt,b12,b1,b1p2,nu,input_Eb)
use inputs, only: q, mO, PI, b0_init, nu_init, km_to_m, mu0,delz, pi, &
ddthickness, magneticShear,ByConeAngle, boundx, quasiparallel
use grid, only: dx_cell, dy_cell, dz_cell,qz
implicit none
real, intent(out):: b0(nx,ny,nz,3), &
bt(nx,ny,nz,3), &
b12(nx,ny,nz,3), &
b1(nx,ny,nz,3), &
b1p2(nx,ny,nz,3), &
nu(nx,ny,nz)
real, intent(inout):: input_Eb
					
real:: eoverm, mO_q, vol, b0eoverm
real:: phi, dtheta, ByConeAngleDelta
integer:: i,j,k,m,Bsetup

!if (quasiparallel .eq. 1) then
!if (ByConeAngle .gt. 1.0) then
!Bsetup= 15
!if (boundx .eq. 5) then
!Bsetup = 17 !2 TD, reflecting, periodic
!endif
!else
!Bsetup= 8
!if (boundx .eq. 5) then
!Bsetup = 16 !2 TD, reflecting, periodic
!endif
!endif
!elseif (quasiparallel .eq. 2) then
!Bsetup = 21
!endif

Bsetup = 24

dtheta = 2*pi / 4 /(2*ddthickness)

eoverm = q/mO
mO_q = mO/q
b0eoverm=b0_init*eoverm
phi = 2.0*PI/180.0
do i=1,nx
do j=1,ny
do k=1,nz   


if (Bsetup .eq. 8) then !BLMN Coordinates, with variable magneticShear from input.dat
!Constant B everywhere, BM
b0(i,j,k,1) = b0_init*eoverm*0.5*cos(magneticShear/180*pi/2)
b0(i,j,k,2) = -(b0_init*eoverm*0.5*sin(magneticShear/180*pi/2))
b0(i,j,k,3) = 0.0
if (k .le. nz/2) then !BL bottom
b0(i,j,k,1) = b0(i,j,k,1) + &
(cos(0.0) - cos(magneticShear/180.0*pi/2.0) ) * b0_init*eoverm*0.5*tanh( ( qz(nz/2)-qz(k))/(ddthickness*delz))
b0(i,j,k,2) = (b0(i,j,k,2) + & 
(sin(magneticShear/180.0*pi/2.0) ) * b0_init*eoverm*0.5*tanh( ( qz(nz/20)-qz(k))/(ddthickness*delz)))
b0(i,j,k,3) = 0.0
endif
if (k .gt. nz/2) then !BL top
b0(i,j,k,1) = b0(i,j,k,1) - &
(cos(magneticShear/180.0*pi/2) - cos(magneticShear/180.0*pi)   ) * b0_init*eoverm*0.5* &
tanh( (qz(k)-qz(nz/2) )/(ddthickness*delz))
b0(i,j,k,2) = (b0(i,j,k,2) - &
(sin(magneticShear/180.0*pi) - sin(magneticShear/180.0*pi/2.0) ) * b0_init*eoverm*0.5* &
tanh( (qz(k)-qz(nz/2) )/(ddthickness*delz)))
b0(i,j,k,3) = 0.0
endif
endif


if (Bsetup .eq. 15) then !BLMN Coordinates,with variable coneAngle By
if (k .le. nz/2) then !BL bottom 
b0(i,j,k,1) =  (cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+90)/180.0*pi))*0.25*b0_init*eoverm + &
(cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+90)/180.0*pi))*&
b0_init*eoverm*0.25*tanh( ( qz(nz/2)-qz(k))/(ddthickness*delz))
b0(i,j,k,2) =  (sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*0.25*b0_init*eoverm + &
(sin(ByConeAngle/180.0*pi)-sin((ByConeAngle+90)/180.0*pi))*&
b0_init*eoverm*0.25*tanh( ( qz(nz/2)-qz(k))/(ddthickness*delz))
b0(i,j,k,3) = 0.0
endif
if (k .gt. nz/2) then !BL top
b0(i,j,k,1) = -(cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+90)/180.0*pi))*0.25*b0_init*eoverm + &
(cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+90)/180.0*pi))*&
b0_init*eoverm*0.25*tanh( (qz(k)-qz(nz/2) )/(ddthickness*delz))
b0(i,j,k,2) =  (sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*0.25*b0_init*eoverm + &
(-sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*&
b0_init*eoverm*0.25*tanh( (qz(k)-qz(nz/2) )/(ddthickness*delz))
b0(i,j,k,3) = 0.0
endif
endif



if (Bsetup .eq. 16) then !BLMN Coordinates, with variable magneticShear from input.dat
!Periodic boundary along z, reflecting x boundary. 2 TDs
b0(i,j,k,1) = b0_init*eoverm*0.5*cos(magneticShear/180*pi/2)
b0(i,j,k,2) = b0_init*eoverm*0.5*sin(magneticShear/180*pi/2)
b0(i,j,k,3) = 0.0

if (k .le. nz/2) then
if (k .le. nz/3) then !BL bottom
b0(i,j,k,1) = b0(i,j,k,1) + &
(cos(0.0) - cos(magneticShear/180.0*pi/2.0) ) * &
b0_init*eoverm*0.5*tanh( ( qz(nz/3)-qz(k))/(ddthickness*delz))
b0(i,j,k,2) = b0(i,j,k,2) - &
(sin(magneticShear/180.0*pi/2.0) ) * &
b0_init*eoverm*0.5*tanh( ( qz(nz/3)-qz(k))/(ddthickness*delz))
b0(i,j,k,3) = 0.0
endif
if (k .gt. nz/3) then !BL top
b0(i,j,k,1) = b0(i,j,k,1) - &
(cos(magneticShear/180.0*pi/2) - cos(magneticShear/180.0*pi) ) * &
b0_init*eoverm*0.5*tanh( (qz(k)-qz(nz/3) )/(ddthickness*delz))
b0(i,j,k,2) = b0(i,j,k,2) + &
(sin(magneticShear/180.0*pi) - sin(magneticShear/180.0*pi/2.0) ) * &
b0_init*eoverm*0.5*tanh( (qz(k)-qz(nz/3) )/(ddthickness*delz))
b0(i,j,k,3) = 0.0
endif
endif

if (k .gt. nz/2) then
if (k .ge. 2*nz/3) then !BL top
b0(i,j,k,1) = b0(i,j,k,1) + &
(sin(magneticShear/180.0*pi) - sin(magneticShear/180.0*pi/2.0) ) * &
b0_init*eoverm*0.5*tanh( (qz(k)-qz(2*nz/3) )/(ddthickness*delz))
b0(i,j,k,2) = b0(i,j,k,2) - &
(cos(magneticShear/180.0*pi/2) - cos(magneticShear/180.0*pi)   ) * &
b0_init*eoverm*0.5*tanh( (qz(k)-qz(2*nz/3) )/(ddthickness*delz))
b0(i,j,k,3) = 0.0
endif
if (k .lt. 2*nz/3) then !BL bottom
b0(i,j,k,1) = b0(i,j,k,1) - &
( sin(magneticShear/180.0*pi/2.0) ) * &
b0_init*eoverm*0.5* &
tanh( ( qz(2*nz/3)-qz(k))/(ddthickness*delz))
b0(i,j,k,2) = b0(i,j,k,2) + &
(cos(0.0) - cos(magneticShear/180.0*pi/2.0) ) * &
b0_init*eoverm*0.5* &
tanh( ( qz(2*nz/3)-qz(k))/(ddthickness*delz))
b0(i,j,k,3) = 0.0
endif
endif

endif

if (Bsetup .eq. 17) then !BLMN Coordinates, with variable coneAngle By
!peridoic z, 2 TD, reflecting boundary
if (k .le. nz/2) then
if (k .le. nz/3) then !BL bottom 
b0(i,j,k,1) = &
(cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+90)/180.0*pi))*&
0.5*b0_init*eoverm + &
(cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+90)/180.0*pi))*&
b0_init*eoverm*0.25*tanh( ( qz(nz/3)-qz(k))/(ddthickness*delz))
b0(i,j,k,2) = &
(sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*&
0.5*b0_init*eoverm + &
(sin(ByConeAngle/180.0*pi)-sin((ByConeAngle+90)/180.0*pi))*&
b0_init*eoverm*0.25*tanh( ( qz(nz/3)-qz(k))/(ddthickness*delz))
b0(i,j,k,3) = 0.0
endif
if (k .gt. nz/3) then !BL top
b0(i,j,k,1) = &
+(cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+90)/180.0*pi))*&
0.25*b0_init*eoverm - &
(cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+90)/180.0*pi))*&
b0_init*eoverm*0.25*tanh( (qz(k)-qz(nz/3) )/(ddthickness*delz))
b0(i,j,k,2) = &
&(sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*&
0.25*b0_init*eoverm + &
(-sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*&
b0_init*eoverm*0.25*tanh( (qz(k)-qz(nz/3) )/(ddthickness*delz))
b0(i,j,k,3) = 0.0
endif
endif

if (k .gt. nz/2) then
if (k .gt. 2*nz/3) then !BL top
b0(i,j,k,1) = &
(-sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*&
0.25*b0_init*eoverm + &
(sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*&
b0_init*eoverm*0.25*tanh( (qz(k)-qz(2*nz/3) )/(ddthickness*delz))
b0(i,j,k,2) = &
(cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+90)/180.0*pi))*&
0.25*b0_init*eoverm - &
(cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+90)/180.0*pi))*&
b0_init*eoverm*0.25*tanh( (qz(k)-qz(2*nz/3) )/(ddthickness*delz))
b0(i,j,k,3) = 0.0
endif
if (k .le. 2*nz/3) then !BL bottom 
b0(i,j,k,1) = &
(-sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*&
0.25*b0_init*eoverm - &
(sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+90)/180.0*pi))*&
b0_init*eoverm*0.25*tanh( ( qz(2*nz/3)-qz(k))/(ddthickness*delz))
b0(i,j,k,2) = &
(cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+90)/180.0*pi))*&
0.25*b0_init*eoverm + &
(cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+90)/180.0*pi))*&
b0_init*eoverm*0.25*tanh( ( qz(2*nz/3)-qz(k))/(ddthickness*delz))
b0(i,j,k,3) = 0.0
endif
endif 

endif




if (Bsetup .eq. 21) then !Bot -By, Top +By, double Qpara
if (k .le. nz/2.0) then !BL bottom 
b0(i,j,k,1) =  2*(cos((ByConeAngle-10)/180.0*pi))*0.25*b0_init*eoverm
b0(i,j,k,2) =  -(sin((ByConeAngle-10)/180.0*pi)+&
sin(((ByConeAngle-10)-2*(ByConeAngle-10))/180.0*pi))*0.25*b0_init*eoverm + &
(-sin((ByConeAngle-10)/180.0*pi)+sin(((ByConeAngle-10)-2*(ByConeAngle-10))/180.0*pi))*b0_init*eoverm*&
0.25*tanh( ( qz(nz/2)-qz(k))/(ddthickness*delz))
b0(i,j,k,3) = 0.0
endif
if (k .gt. nz/2.0) then !BL top
b0(i,j,k,1) =  2*(cos((ByConeAngle-0)/180.0*pi))*0.25*b0_init*eoverm
b0(i,j,k,2) =   (sin((ByConeAngle-0)/180.0*pi)+ &
sin(((ByConeAngle-0)-2*(ByConeAngle-0))/180.0*pi))*0.25*b0_init*eoverm + &
(sin((ByConeAngle-0)/180.0*pi)-sin(((ByConeAngle-0)-2*(ByConeAngle-0))/180.0*pi))*b0_init*eoverm*&
0.25*tanh( (qz(k)-qz(nz/2) )/(ddthickness*delz))
b0(i,j,k,3) = 0.0
endif
endif


!2TD periodic for Bow shock
if (Bsetup .eq. 22) then !BLMN Coordinates, with variable coneAngle By
    !peridoic z, 2 TD, reflecting boundary
    if (k .le. nz/2) then
        if (k .le. nz/3) then !BL bottom 
            b0(i,j,k,1) = 2*(cos((ByConeAngle-0)/180.0*pi))*0.25*b0_init*eoverm
            b0(i,j,k,2) = -(sin((ByConeAngle-0)/180.0*pi) + &
            sin(((ByConeAngle-0)-2*(ByConeAngle-0))/180.0*pi))*0.25*b0_init*eoverm + &
            (-sin((ByConeAngle-0)/180.0*pi)+sin(((ByConeAngle-0)-2*(ByConeAngle-0))/180.0*pi))*&
            b0_init*eoverm*0.25*tanh( ( qz(nz/3)-qz(k) )/(ddthickness*delz))
            b0(i,j,k,3) = 0.0
        endif
        if (k .gt. nz/3) then !BL top
            b0(i,j,k,1) = 2*(cos((ByConeAngle-0)/180.0*pi))*0.25*b0_init*eoverm
            b0(i,j,k,2) = +(sin((ByConeAngle-0)/180.0*pi) + &
            sin(((ByConeAngle-0)-2*(ByConeAngle-0))/180.0*pi))*0.25*b0_init*eoverm + &
            (+sin((ByConeAngle-0)/180.0*pi)-sin(((ByConeAngle-0)-2*(ByConeAngle-0))/180.0*pi))*&
            b0_init*eoverm*0.25*tanh( (qz(k)-qz(nz/3) )/(ddthickness*delz))
            b0(i,j,k,3) = 0.0
        endif
    endif
        
    if (k .gt. nz/2) then
        if (k .gt. 2*nz/3) then !BL top
            b0(i,j,k,1) = 2*(cos((ByConeAngle-0)/180.0*pi))*0.25*b0_init*eoverm
            b0(i,j,k,2) = -(sin((ByConeAngle-0)/180.0*pi) + &
            sin(((ByConeAngle-0)-2*(ByConeAngle-0))/180.0*pi))*0.25*b0_init*eoverm + &
            (-sin((ByConeAngle-0)/180.0*pi)+sin(((ByConeAngle-0)-2*(ByConeAngle-0))/180.0*pi))*&
            b0_init*eoverm*0.25*tanh( (qz(k)-qz(2*nz/3) )/(ddthickness*delz))
            b0(i,j,k,3) = 0.0
        endif
        if (k .le. 2*nz/3) then !BL bottom 
            b0(i,j,k,1) = 2*(cos((ByConeAngle-0)/180.0*pi))*0.25*b0_init*eoverm
            b0(i,j,k,2) = +(sin((ByConeAngle-0)/180.0*pi) + &
            sin(((ByConeAngle-0)-2*(ByConeAngle-0))/180.0*pi))*0.25*b0_init*eoverm + &
            (+sin((ByConeAngle-0)/180.0*pi)-sin(((ByConeAngle-0)-2*(ByConeAngle-0))/180.0*pi))*&
            b0_init*eoverm*0.25*tanh( ( qz(2*nz/3)-qz(k))/(ddthickness*delz))
            b0(i,j,k,3) = 0.0
        endif

    endif 
    
endif



if (Bsetup .eq. 23) then !BLMN Coordinates, with variable coneAngle By
    !peridoic z, 2 TD, reflecting boundary
    ByConeAngle = -15.0
    ByConeAngleDelta = +30.0
    if (k .le. nz/2) then
    if (k .lt. nz/3) then !BL bottom 
    b0(i,j,k,1) = &
    (cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+ByConeAngleDelta)/180.0*pi))*&
    0.25*b0_init*eoverm + &
    (cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+ByConeAngleDelta)/180.0*pi))*&
    b0_init*eoverm*0.25*tanh( ( qz(nz/3)-qz(k))/(ddthickness*delz))

    b0(i,j,k,2) = &
    (sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+ByConeAngleDelta)/180.0*pi))*&
    0.25*b0_init*eoverm + &
    (sin(ByConeAngle/180.0*pi)-sin((ByConeAngle+ByConeAngleDelta)/180.0*pi))*&
    b0_init*eoverm*0.25*tanh( ( qz(nz/3)-qz(k)-qz(1))/(ddthickness*delz))
    
    !if (b0(i,j,k,2) > b0_init*eoverm) then
    !    b0(i,j,k,2) = b0_init*eoverm
    !endif

    !if (b0(i,j,k,2) < -b0_init*eoverm) then
    !    b0(i,j,k,2) = -b0_init*eoverm
    !endif

    b0(i,j,k,3) = 0.0
    endif
    if (k .ge. nz/3) then !BL top
    b0(i,j,k,1) = &
    +(cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+ByConeAngleDelta)/180.0*pi))*&
    0.25*b0_init*eoverm - &
    (cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+ByConeAngleDelta)/180.0*pi))*&
    b0_init*eoverm*0.25*tanh( (qz(k)-qz(nz/3) )/(ddthickness*delz))

    b0(i,j,k,2) = &
    (sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+ByConeAngleDelta)/180.0*pi))*&
    0.25*b0_init*eoverm + &
    (-sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+ByConeAngleDelta)/180.0*pi))*&
    b0_init*eoverm*0.25*tanh( (qz(k)-qz(nz/3)-qz(1) )/(ddthickness*delz))

    !if (b0(i,j,k,2) > b0_init*eoverm) then
    !    b0(i,j,k,2) = b0_init*eoverm
    !endif

    !if (b0(i,j,k,2) < -b0_init*eoverm) then
    !    b0(i,j,k,2) = -b0_init*eoverm
    !endif

    b0(i,j,k,3) = 0.0
    endif
    endif
    
    if (k .gt. nz/2) then
    if (k .gt. 2*nz/3) then !BL top
    b0(i,j,k,1) = &
    (cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+ByConeAngleDelta)/180.0*pi))*&
    0.25*b0_init*eoverm + &
    (cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+ByConeAngleDelta)/180.0*pi))*&
    b0_init*eoverm*0.25*tanh( (qz(k)-qz(2*nz/3) )/(ddthickness*delz))
    
    b0(i,j,k,2) = &
    (sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+ByConeAngleDelta)/180.0*pi))*&
    0.25*b0_init*eoverm + &
    (sin(ByConeAngle/180.0*pi)-sin((ByConeAngle+ByConeAngleDelta)/180.0*pi))*&
    b0_init*eoverm*0.25*tanh( (qz(k)-qz(2*nz/3)-qz(1) )/(ddthickness*delz))
    
    !if (b0(i,j,k,2) > b0_init*eoverm) then
    !    b0(i,j,k,2) = b0_init*eoverm
    !endif

    !if (b0(i,j,k,2) < -b0_init*eoverm) then
    !    b0(i,j,k,2) = -b0_init*eoverm
    !endif

    b0(i,j,k,3) = 0.0
    endif
    
    if (k .le. 2*nz/3) then !BL bottom 
    b0(i,j,k,1) = &
    +(cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+ByConeAngleDelta)/180.0*pi))*&
    0.25*b0_init*eoverm - &
    (cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+ByConeAngleDelta)/180.0*pi))*&
    b0_init*eoverm*0.25*tanh( ( qz(2*nz/3)-qz(k))/(ddthickness*delz))
    
    b0(i,j,k,2) = &
    (sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+ByConeAngleDelta)/180.0*pi))*&
    0.25*b0_init*eoverm + &
    (-sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+ByConeAngleDelta)/180.0*pi))*&
    b0_init*eoverm*0.25*tanh( ( qz(2*nz/3)-qz(k)-qz(1))/(ddthickness*delz))
    
    !if (b0(i,j,k,2) > b0_init*eoverm) then
    !    b0(i,j,k,2) = b0_init*eoverm
    !endif

    !if (b0(i,j,k,2) < -b0_init*eoverm) then
    !    b0(i,j,k,2) = -b0_init*eoverm
    !endif

    b0(i,j,k,3) = 0.0
    endif
    endif 
    
    endif


    if (Bsetup .eq. 24) then !BLMN Coordinates, with variable coneAngle By
        if (i .gt. nx-150) then
            b0(i,j,k,1) = 0.0!b0_init*eoverm*0.5
            b0(i,j,k,2) = b0_init*eoverm*0.5
            b0(i,j,k,3) = 0.0
        endif

        if (i .le. nx-150) then
        !peridoic z, 2 TD, reflecting boundary
        ByConeAngle = -15.0
        ByConeAngleDelta = +30.0
            if (k .le. nz/2) then
                if (k .lt. nz/3) then !BL bottom 
                    b0(i,j,k,1) = &
                    (cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+ByConeAngleDelta)/180.0*pi))*&
                    0.25*b0_init*eoverm + &
                    (cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+ByConeAngleDelta)/180.0*pi))*&
                    b0_init*eoverm*0.25*tanh( ( qz(nz/3)-qz(k))/(ddthickness*delz))
                
                    b0(i,j,k,2) = &
                    (sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+ByConeAngleDelta)/180.0*pi))*&
                    0.25*b0_init*eoverm + &
                    (sin(ByConeAngle/180.0*pi)-sin((ByConeAngle+ByConeAngleDelta)/180.0*pi))*&
                    b0_init*eoverm*0.25*tanh( ( qz(nz/3)-qz(k)-qz(1))/(ddthickness*delz))
                    
                    !if (b0(i,j,k,2) > b0_init*eoverm) then
                    !    b0(i,j,k,2) = b0_init*eoverm
                    !endif
                
                    !if (b0(i,j,k,2) < -b0_init*eoverm) then
                    !    b0(i,j,k,2) = -b0_init*eoverm
                    !endif
                
                    b0(i,j,k,3) = 0.0
                endif
                if (k .ge. nz/3) then !BL top
                    b0(i,j,k,1) = &
                    +(cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+ByConeAngleDelta)/180.0*pi))*&
                    0.25*b0_init*eoverm - &
                    (cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+ByConeAngleDelta)/180.0*pi))*&
                    b0_init*eoverm*0.25*tanh( (qz(k)-qz(nz/3) )/(ddthickness*delz))
                
                    b0(i,j,k,2) = &
                    (sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+ByConeAngleDelta)/180.0*pi))*&
                    0.25*b0_init*eoverm + &
                    (-sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+ByConeAngleDelta)/180.0*pi))*&
                    b0_init*eoverm*0.25*tanh( (qz(k)-qz(nz/3)-qz(1) )/(ddthickness*delz))
                
                    !if (b0(i,j,k,2) > b0_init*eoverm) then
                    !    b0(i,j,k,2) = b0_init*eoverm
                    !endif
                
                    !if (b0(i,j,k,2) < -b0_init*eoverm) then
                    !    b0(i,j,k,2) = -b0_init*eoverm
                    !endif
    
                    b0(i,j,k,3) = 0.0
                endif
            endif
        
            if (k .gt. nz/2) then
                if (k .gt. 2*nz/3) then !BL top
                    b0(i,j,k,1) = &
                    (cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+ByConeAngleDelta)/180.0*pi))*&
                    0.25*b0_init*eoverm + &
                    (cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+ByConeAngleDelta)/180.0*pi))*&
                    b0_init*eoverm*0.25*tanh( (qz(k)-qz(2*nz/3) )/(ddthickness*delz))
                    
                    b0(i,j,k,2) = &
                    (sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+ByConeAngleDelta)/180.0*pi))*&
                    0.25*b0_init*eoverm + &
                    (sin(ByConeAngle/180.0*pi)-sin((ByConeAngle+ByConeAngleDelta)/180.0*pi))*&
                    b0_init*eoverm*0.25*tanh( (qz(k)-qz(2*nz/3)-qz(1) )/(ddthickness*delz))
                    
                    !if (b0(i,j,k,2) > b0_init*eoverm) then
                    !    b0(i,j,k,2) = b0_init*eoverm
                    !endif
                
                    !if (b0(i,j,k,2) < -b0_init*eoverm) then
                    !    b0(i,j,k,2) = -b0_init*eoverm
                    !endif
                
                    b0(i,j,k,3) = 0.0
                endif
        
                if (k .le. 2*nz/3) then !BL bottom 
                    b0(i,j,k,1) = &
                    +(cos(ByConeAngle/180.0*pi)+cos((ByConeAngle+ByConeAngleDelta)/180.0*pi))*&
                    0.25*b0_init*eoverm - &
                    (cos(ByConeAngle/180.0*pi)-cos((ByConeAngle+ByConeAngleDelta)/180.0*pi))*&
                    b0_init*eoverm*0.25*tanh( ( qz(2*nz/3)-qz(k))/(ddthickness*delz))
                    
                    b0(i,j,k,2) = &
                    (sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+ByConeAngleDelta)/180.0*pi))*&
                    0.25*b0_init*eoverm + &
                    (-sin(ByConeAngle/180.0*pi)+sin((ByConeAngle+ByConeAngleDelta)/180.0*pi))*&
                    b0_init*eoverm*0.25*tanh( ( qz(2*nz/3)-qz(k)-qz(1))/(ddthickness*delz))
                    
                    !if (b0(i,j,k,2) > b0_init*eoverm) then
                    !    b0(i,j,k,2) = b0_init*eoverm
                    !endif
                
                    !if (b0(i,j,k,2) < -b0_init*eoverm) then
                    !    b0(i,j,k,2) = -b0_init*eoverm
                    !endif
                
                    b0(i,j,k,3) = 0.0
                endif
            endif 
        endif
    endif






enddo
enddo
enddo

do i=1,nx
do j=1,ny
do k= 1,nz
nu(i,j,k) = nu_init
do m = 1,3
bt(i,j,k,m) = b0(i,j,k,m)
b12(i,j,k,m) = b0(i,j,k,m)
b1(i,j,k,m) = b0(i,j,k,m)
b1p2(i,j,k,m) = b0(i,j,k,m)
vol = dx_cell(i)*dy_cell(j)*dz_cell(k)*km_to_m**3
input_Eb = input_Eb + (vol/(2.0*mu0))*(mO_q*b0(i,j,k,m))**2
enddo
enddo
enddo
enddo

end subroutine grd6_setup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

subroutine grd7()
use grid
use mult_proc, only: my_rank
use inputs, only: dx,dy,delz,out_dir
implicit none
integer:: i,j,k
real:: zplus,zminus,xplus,xminus,yplus,yminus

!!!!!!!!!Unstretched grids!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

do j=1,ny
qy(j) = j*dy
dy_grid(j) = dy
enddo

do i = 1,nx
qx(i) = i*dx
dx_grid(i) = dx
enddo

do k = 1,nz
qz(k) = k*delz
dz_grid(k) = delz
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

dz_cell(1) = dz_grid(1)
dz_cell(nz) = dz_grid(nz)
zrat(1) = 0.5
zrat(nz) = 0.5
do k=2, nz-1
dz_cell(k) = ((qz(k+1) + qz(k))/2.0) - ((qz(k) + qz(k-1))/2.0)
zplus = (qz(k+1) + qz(k))/2.0
zminus = (qz(k) + qz(k-1))/2.0
zrat(k) = (qz(k) - zminus)/(zplus - zminus)
enddo

dx_cell(1) = dx_grid(1)
dx_cell(nx) = dx_grid(nx)
xrat(1) = 0.5
xrat(nx) = 0.5
do i=2, nx-1
dx_cell(i) = ((qx(i+1) + qx(i))/2.0) - ((qx(i) + qx(i-1))/2.0)
xplus = (qx(i+1) + qx(i))/2.0
xminus = (qx(i) + qx(i-1))/2.0
xrat(i) = (qx(i) - xminus)/(xplus - xminus)
enddo

dy_cell(1) = dy_grid(1)
dy_cell(ny) = dy_grid(ny)
yrat(1) = 0.5
yrat(ny) = 0.5
do j=2, ny-1
dy_cell(j) = ((qy(j+1) + qy(j))/2.0) - ((qy(j) + qy(j-1))/2.0)
yplus = (qy(j+1) + qy(j))/2.0
yminus = (qy(j) + qy(j-1))/2.0
yrat(j) = (qy(j) - yminus)/(yplus - yminus)
enddo

!!!!!!!!!Print Coordinates!!!!!!!!!!!!!!!!!!!!!!

if (my_rank .eq. 0) then
open(40,file=trim(out_dir)//'c.coord.dat',status='unknown',form='unformatted')

write(40) nx
write(40) ny
write(40) nz
write(40) qx
write(40) qy
write(40) qz
write(40) dz_grid
write(40) dz_cell
close(40)
		
endif
		
end subroutine grd7
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
end module initial
            
