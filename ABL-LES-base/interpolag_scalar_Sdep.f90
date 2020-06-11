subroutine interpolag_scalar_Sdep()
! Subroutine written by Stimit Shah on 27th Sept 2010
! This subroutine computes the values of G_LM and G_MM 
! at positions (x-u*dt) for use in the subroutine lagrangian
use types,only:rprec
use param
use sgsmodule

implicit none
integer::jx,jy,jz 
integer::jjx,jjy,jjz
integer::jz_min

$if ($MPI)
  $define $lbz 0
  integer::recvcounts(nproc)
  integer::displs(nproc)
$else
  $define $lbz 1
$endif

real(kind=rprec)::frac_x,frac_y,frac_z
real(kind=rprec)::comp_x,comp_y,comp_z
real(kind=rprec),save,dimension(nx+2,ny+2,nz+2)::GG_LM,GG_MM
real(kind=rprec),save,dimension(nx+2,ny+2,nz+2)::GG_QN,GG_NN

integer::addx,addy,addz
integer::jxaddx,jyaddy,jzaddz
real(kind=rprec),save,dimension(nx+2,ny+2,nz+2)::Beta_t

do jz=1,nz
do jy=1,ny
do jx=1,nx
   GG_LM(jx+1,jy+1,jz+1) = G_LM(jx,jy,jz)
   GG_MM(jx+1,jy+1,jz+1) = G_MM(jx,jy,jz)
   GG_QN(jx+1,jy+1,jz+1) = G_QN(jx,jy,jz)
   GG_NN(jx+1,jy+1,jz+1) = G_NN(jx,jy,jz)
end do
end do
end do

! This is a bit like witch craft but the following lines 
! do take care of all the edges including the corners
if ( read_inflow_file ) then
   GG_LM(1,2:ny+1,2:nz+1) = GLM_hold
   GG_LM(nx+2,:,:) = GG_LM(nx+1,:,:)
else
   GG_LM(1,:,:) = GG_LM(nx+1,:,:)
   GG_LM(nx+2,:,:) = GG_LM(2,:,:)
end if

GG_LM(:,1,:) = GG_LM(:,ny+1,:) 
GG_LM(:,ny+2,:) = GG_LM(:,2,:) 

$if ($MPI)
  !--send G_LM @ jz=nz-1 to G_LM @ jz=0'
  !  i.e. GG_LM @ jz=nz to GG_LM @ jz=1'
  call mpi_sendrecv (GG_LM(1, 1, nz), (nx+2)*(ny+2), MPI_RPREC, up, 1,   &
                     GG_LM(1, 1, 1), (nx+2)*(ny+2), MPI_RPREC, down, 1,  &
                     comm, status, ierr)
  !--G_LM @ jz=nz and G_LM @ jz=1' should already be in sync (test?)
  !  i.e. GG_LM @ jz=nz+1 and GG_LM @ jz=2'
$endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  GG_LM(:,:,1) = GG_LM(:,:,2)
end if

$if ($MPI)
  !--send G_LM @ jz=2 to G_LM @ jz=nz+1'
  !  i.e. GG_LM @ jz=3 to GG_LM @ jz=nz+2'
  call mpi_sendrecv (GG_LM(1, 1, 3), (nx+2)*(ny+2), MPI_RPREC, down, 2,   &
                     GG_LM(1, 1, nz+2), (nx+2)*(ny+2), MPI_RPREC, up, 2,  &
                     comm, status, ierr)
$endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
  GG_LM(:,:,nz+2) = GG_LM(:,:,nz+1)
end if

if (read_inflow_file) then
  GG_MM(1,2:ny+1,2:nz+1)=GMM_hold
  GG_MM(nx+2,:,:) = GG_MM(nx+1,:,:) 
else
  GG_MM(1,:,:) = GG_MM(nx+1,:,:)
  GG_MM(nx+2,:,:) = GG_MM(2,:,:)
end if
GG_MM(:,1,:) = GG_MM(:,ny+1,:) 
GG_MM(:,ny+2,:) = GG_MM(:,2,:) 

$if ($MPI)
  !--send G_MM @ jz=nz-1 to G_MM @ jz=0'
  !  i.e. GG_MM @ jz=nz to GG_MM @ jz=1'
  call mpi_sendrecv (GG_MM(1, 1, nz), (nx+2)*(ny+2), MPI_RPREC, up, 3,   &
                     GG_MM(1, 1, 1), (nx+2)*(ny+2), MPI_RPREC, down, 3,  &
                     comm, status, ierr)
  !--G_MM @ jz=nz and G_MM @ jz=1' should already be in sync (test?)
  !  i.e. GG_MM @ jz=nz+1 and GG_MM @ jz=2'
$endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  GG_MM(:,:,1) = GG_MM(:,:,2)
end if

$if ($MPI)
  !--send G_MM @ jz=2 to G_MM @ jz=nz+1'
  !  i.e. GG_MM @ jz=3 to GG_MM @ jz=nz+2'
  call mpi_sendrecv (GG_MM(1, 1, 3), (nx+2)*(ny+2), MPI_RPREC, down, 4,   &
                     GG_MM(1, 1, nz+2), (nx+2)*(ny+2), MPI_RPREC, up, 4,  &
                     comm, status, ierr)
$endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
  GG_MM(:,:,nz+2) = GG_MM(:,:,nz+1) 
end if

if (read_inflow_file) then
  GG_QN(1,2:ny+1,2:nz+1)= GQN_hold
  GG_QN(nx+2,:,:) = GG_QN(nx+1,:,:)
else
  GG_QN(1,:,:) = GG_QN(nx+1,:,:)
  GG_QN(nx+2,:,:) = GG_QN(2,:,:)
end if
GG_QN(:,1,:) = GG_QN(:,ny+1,:) 
GG_QN(:,ny+2,:) = GG_QN(:,2,:) 

$if ($MPI)
  !--send G_QN @ jz=nz-1 to G_QN @ jz=0'
  !  i.e. GG_QN @ jz=nz to GG_QN @ jz=1'
  call mpi_sendrecv (GG_QN(1, 1, nz), (nx+2)*(ny+2), MPI_RPREC, up, 5,   &
                     GG_QN(1, 1, 1), (nx+2)*(ny+2), MPI_RPREC, down, 5,  &
                     comm, status, ierr)
  !--G_QN @ jz=nz and G_QN @ jz=1' should already be in sync (test?)
  !  i.e. GG_QN @ jz=nz+1 and GG_QN @ jz=2'
$endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  GG_QN(:,:,1) = GG_QN(:,:,2)
end if

$if ($MPI)
  !--send G_QN @ jz=2 to G_QN @ jz=nz+1'
  !  i.e. GG_QN @ jz=3 to GG_QN @ jz=nz+2'
  call mpi_sendrecv (GG_QN(1, 1, 3), (nx+2)*(ny+2), MPI_RPREC, down, 6,   &
                     GG_QN(1, 1, nz+2), (nx+2)*(ny+2), MPI_RPREC, up, 6,  &
                     comm, status, ierr)
$endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
  GG_QN(:,:,nz+2) = GG_QN(:,:,nz+1) 
end if

if ( read_inflow_file ) then
  GG_NN(1,2:ny+1,2:nz+1) = GNN_hold
  GG_NN(nx+2,:,:) = GG_NN(nx+1,:,:) 
else
  GG_NN(1,:,:) = GG_NN(nx+1,:,:)
  GG_NN(nx+2,:,:) = GG_NN(2,:,:)
end if
GG_NN(:,1,:) = GG_NN(:,ny+1,:) 
GG_NN(:,ny+2,:) = GG_NN(:,2,:) 

$if ($MPI)
  !--send G_NN @ jz=nz-1 to G_NN @ jz=0'
  !  i.e. GG_NN @ jz=nz to GG_NN @ jz=1'
  call mpi_sendrecv (GG_NN(1, 1, nz), (nx+2)*(ny+2), MPI_RPREC, up, 7,   &
                     GG_NN(1, 1, 1), (nx+2)*(ny+2), MPI_RPREC, down, 7,  &
                     comm, status, ierr)
  !--G_NN @ jz=nz and G_NN @ jz=1' should already be in sync (test?)
  !  i.e. GG_NN @ jz=nz+1 and GG_NN @ jz=2'
$endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  GG_NN(:,:,1) = GG_NN(:,:,2)
end if

$if ($MPI)
  !--send G_NN @ jz=2 to G_NN @ jz=nz+1'
  !  i.e. GG_NN @ jz=3 to GG_NN @ jz=nz+2'
  call mpi_sendrecv (GG_NN(1, 1, 3), (nx+2)*(ny+2), MPI_RPREC, down, 8,   &
                     GG_NN(1, 1, nz+2), (nx+2)*(ny+2), MPI_RPREC, up, 8,  &
                     comm, status, ierr)
$endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
  GG_NN(:,:,nz+2) = GG_NN(:,:,nz+1) 
end if
! end of witch craft

u_lag_s = u_lag_s/real(cs_count,kind=rprec)
v_lag_s = v_lag_s/real(cs_count,kind=rprec)

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  jz_min = 2
else
  jz_min = 1
end if

! must do backwards due to data depenency on plane below current one
do jz = nz,jz_min, -1
  u_lag_s(:,:,jz) = 0.5_rprec*(u_lag_s(:,:,jz) + u_lag_s(:,:,jz-1))
  v_lag_s(:,:,jz) = 0.5_rprec*(v_lag_s(:,:,jz) + v_lag_s(:,:,jz-1))
end do

w_lag_s = w_lag_s/real(cs_count,kind=rprec)

$if ($MPI)
  w_lag_s(:, :, $lbz) = BOGUS
$endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  w_lag_s (:,:,1) = 0.25_rprec*w_lag_s (:,:,2)
end if

! computes the 3-D inverse displacement arrays that describe
! the location where the point was at the previous step
u_lag_s = -u_lag_s*dt*real(cs_count,kind=rprec)/dx
v_lag_s = -v_lag_s*dt*real(cs_count,kind=rprec)/dy   
w_lag_s = -w_lag_s*dt*real(cs_count,kind=rprec)/dz

$if ($MPI)
  u_lag_s(:,:,$lbz) = BOGUS
  v_lag_s(:,:,$lbz) = BOGUS
  w_lag_s(:,:,$lbz) = BOGUS
$endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  ! because first plane is on u,v,p nodes
  ! this corrects for the fact that the first cell in the 
  ! z-direction has height dz/2, it doubles the zp fraction
  ! if this fraction relates to the cell beneath it
  do jy=1,ny
  do jx=1,nx
     w_lag_s(jx,jy,2) = w_lag_s(jx,jy,2) + min (w_lag_s(jx,jy,2), 0._rprec)
     w_lag_s(jx,jy,1) = w_lag_s(jx,jy,1) + max (w_lag_s(jx,jy,1), 0._rprec)
  end do
  end do

  w_lag_s(:,:,2) = min (1._rprec, w_lag_s(:,:,2))
  w_lag_s(:,:,1) = min (1._rprec, w_lag_s(:,:,2))
  w_lag_s(:,:,2) = max (-1._rprec, w_lag_s(:,:,2))
  w_lag_s(:,:,1) = max (-1._rprec, w_lag_s(:,:,2))
end if

! if(mod(jt,c_count).eq.0) print*,'Lagrangian CFL condition= ', &
!                          maxval(abs(u_lag_s(1:nx,:,1:nz)))

do jz=1,nz
 jjz = jz+1
 do jy=1,ny
  jjy = jy+1
   do jx=1,nx
      jjx = jx+1

      addx = int(sign(1._rprec,u_lag_s(jx,jy,jz)))
      addy = int(sign(1._rprec,v_lag_s(jx,jy,jz)))
      addz = int(sign(1._rprec,w_lag_s(jx,jy,jz)))
      jxaddx = jjx + addx 
      jyaddy = jjy + addy
      jzaddz = jjz + addz

      comp_x = abs(u_lag_s(jx,jy,jz))
      comp_y = abs(v_lag_s(jx,jy,jz))
      comp_z = abs(w_lag_s(jx,jy,jz)) 
      frac_x = 1._rprec - comp_x
      frac_y = 1._rprec - comp_y
      frac_z = 1._rprec - comp_z

      ! computes interpolated G_LM
      G_LM(jx,jy,jz)=frac_x*frac_y*&
          (GG_LM(jjx,jjy,jjz)*frac_z+GG_LM(jjx,jjy,jzaddz)*comp_z)&
          + frac_x*comp_y*&
          (GG_LM(jjx,jyaddy,jjz)*frac_z+GG_LM(jjx,jyaddy,jzaddz)*comp_z)&
          + comp_x*frac_y*&
          (GG_LM(jxaddx,jjy,jjz)*frac_z+GG_LM(jxaddx,jjy,jzaddz)*comp_z)&
          + comp_x*comp_y*&
          (GG_LM(jxaddx,jyaddy,jjz)*frac_z+GG_LM(jxaddx,jyaddy,jzaddz)*comp_z)
 
      ! computes interpolated G_MM
      G_MM(jx,jy,jz)=frac_x*frac_y*&
          (GG_MM(jjx,jjy,jjz)*frac_z+GG_MM(jjx,jjy,jzaddz)*comp_z)&
          + frac_x*comp_y*&
          (GG_MM(jjx,jyaddy,jjz)*frac_z+GG_MM(jjx,jyaddy,jzaddz)*comp_z)&
          + comp_x*frac_y*&
          (GG_MM(jxaddx,jjy,jjz)*frac_z+GG_MM(jxaddx,jjy,jzaddz)*comp_z)&
          + comp_x*comp_y*&
          (GG_MM(jxaddx,jyaddy,jjz)*frac_z+GG_MM(jxaddx,jyaddy,jzaddz)*comp_z)

      ! computes interpolated G_QN
      G_QN(jx,jy,jz)=frac_x*frac_y*&
          (GG_QN(jjx,jjy,jjz)*frac_z+GG_QN(jjx,jjy,jzaddz)*comp_z)&
          + frac_x*comp_y*&
          (GG_QN(jjx,jyaddy,jjz)*frac_z+GG_QN(jjx,jyaddy,jzaddz)*comp_z)&
          + comp_x*frac_y*&
          (GG_QN(jxaddx,jjy,jjz)*frac_z+GG_QN(jxaddx,jjy,jzaddz)*comp_z)&
          + comp_x*comp_y*&
          (GG_QN(jxaddx,jyaddy,jjz)*frac_z+GG_QN(jxaddx,jyaddy,jzaddz)*comp_z)

      ! computes interpolated G_NN
      G_NN(jx,jy,jz)=frac_x*frac_y*&
          (GG_NN(jjx,jjy,jjz)*frac_z+GG_NN(jjx,jjy,jzaddz)*comp_z)&
          + frac_x*comp_y*&
          (GG_NN(jjx,jyaddy,jjz)*frac_z+GG_NN(jjx,jyaddy,jzaddz)*comp_z)&
          + comp_x*frac_y*&
          (GG_NN(jxaddx,jjy,jjz)*frac_z+GG_NN(jxaddx,jjy,jzaddz)*comp_z)&
          + comp_x*comp_y*&
          (GG_NN(jxaddx,jyaddy,jjz)*frac_z+GG_NN(jxaddx,jyaddy,jzaddz)*comp_z)
end do
end do
end do

u_lag_s = 0._rprec
v_lag_s = 0._rprec
w_lag_s = 0._rprec

end subroutine interpolag_scalar_Sdep
