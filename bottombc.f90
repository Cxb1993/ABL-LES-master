module bottombc
use types,only:rprec
use param,only:nx,ny,ld
implicit none
! num_patch = # of patches, zo=surface roughness for the patches types
! ptypes=number of patches types to be used, usually we use 2

! for slices of patches
integer,parameter::num_patch=2,ptypes=2
real(kind=rprec),parameter::zo1=0.0002_rprec,zo2=0.0002_rprec
real(kind=rprec),parameter::theta_s1=274.15_rprec,theta_s2=274.15_rprec

!JJF this is for reading sfc parameters from a .txt file
integer,parameter::remote_to_ice=1
integer::row,col,max_rows,max_cols

! if remote_to_ice is 1, then bottom BC is set using outside file

! this is for a square patch
integer,parameter::square=0
real(kind=rprec),parameter::zo_out=0.000010_rprec, zo_square=0.000010_rprec
integer,parameter::Nx_sq=5,Ny_sq=3

! This is the common local ustar computed in obukhov and used in wallstress.f90
real(kind=rprec),dimension(nx,ny)::ustar_avg
real(kind=rprec),dimension(nx,ny)::zo,T_s,q_s
real(kind=rprec),dimension(nx,ny)::z_os

!TS add for non-neutral case
real(kind=rprec),dimension(nx,ny)::phi_m,psi_m,phi_h,psi_h

!VK The obukhov similarity functions are computed using obukhov(scalars_module.f90) for non-neutral scenario
integer,dimension(nx,ny)::patch
integer,dimension(num_patch)::patchnum


contains
subroutine patches()
!VK
! This assigns momentum roughness, temperature and wetness for the different patches
! and fills the lookup tables patch and patchnum. This is called from routines
! patch_or_remote.f90 depending whether to use remotely-sensed data or patches

use param
implicit none
integer::i, j, patchtype, begini, endi, type

! sc:
! without this, some compilers may give total junk
! this was cause the u_avg to be WAY off, and so txz was WAY off, etc.
patchnum=1

if (square==1) then !for square patch
   zo(:,:)=zo_out/z_i
   patch(:,:)= 1
   zo(Nx/2-(Nx_sq-1)/2:Nx/2+(Nx_sq-1)/2,&
        Ny/2-(Ny_sq-1)/2:Ny/2+(Ny_sq-1)/2) = zo_square/z_i
   patch(Nx/2-(Nx_sq-1)/2:Nx/2+(Nx_sq-1)/2,&
        Ny/2-(Ny_sq-1)/2:Ny/2+(Ny_sq-1)/2) = 2
   patchnum=Nx*Ny

! JJF - using a .txt file for surface files
else if (remote_to_ice==1) then
  max_rows=nx
  max_cols=ny
  ! opening/reading the sfc temp file
  open (11, file = 'T_s_remote_ice.txt')
  do j = 1,ny
    read(11,*) (T_s(i,j),i=1,nx)
  end do
  close(11)
  T_s = T_s*(1.0/T_scale)
  ! opening/reading the roughness file
  open (12, file = 'zo_remote_ice.txt')
  do j = 1,ny
    read(12,*) (zo(i,j),i=1,nx)
  end do
  close(12)
  zo = zo*(1.0/z_i)
  ! opening/reading the sfc humidity temp file
!  open (13, file = 'q_s_remote_ice.txt')
!  do j = 1,ny
!    read(13,*) (q_s(i,j),i=1,nx)
!  end do
!  close(13)

! for regular patches
else
   do j=1,ny
      endi = 0
      type=1
      do patchtype=1,num_patch
         begini=endi+1
         endi=(nx*patchtype)/num_patch
         do i=begini,endi
            if (type.eq.1) then
               zo(i,j)=zo1/z_i
               T_s(i,j)=theta_s1/t_scale
              ! q_s(i,j)=q_s1
               patch(i,j)=type
            end if
            if(type.eq.2) then
               zo(i,j)=zo2/z_i
               T_s(i,j)=theta_s2/t_scale
              ! q_s(i,j)=q_s2
               patch(i,j)=type
            end if
            patchnum(type)=patchnum(type)+1
          end do
          if (type.eq.1) then
             type=2
          elseif (type.eq.2) then
             type=1
          end if
       end do
      end do
end if

!JJF - for testing purposes -  output arrays of T_s and zo_s after patch formation
! write temp
open(unit=14, file="T_s_init.txt",action="write", status="replace")
do j=1,ny
  write(14,*)(T_s(i,j),i=1,nx)
end do
! write zo
open(unit=15, file="zo_init.txt",action="write", status="replace")
do j=1,ny
  write(15,*)(zo(i,j),i=1,nx)
end do
! write q_s
!open(unit=16, file="q_s_init.txt",action="write", status="replace")
!do j=1,ny
!  write(16,*)(q_s(i,j),i=1,nx)
!end do

end subroutine patches

!!!xxxxxxxxxxxx--------VIJ-------------XXXXXXXXXXXXXXXXXXX
subroutine remote_to_patch(T_s_in,T_or_z)
implicit none
integer:: ii,jj,kk,T_or_z
real(kind=rprec),dimension(nx,ny):: T_s_in,dummy
real(kind=rprec):: sigma_multiplier,crap1,crap2,crap3
real(kind=rprec):: patch_cold, patch_hot

sigma_multiplier=1._rprec
dummy=0._rprec
if (T_or_z == 1) then
    crap1=sum(T_s_in)/float(nx*ny)
elseif (T_or_z == 2) then
    crap1=exp(sum(dlog(T_s_in))/float(nx*ny))
else
    print *,'Wrong choice of T_or_z in remote_to_patch().. STOPPING'
    stop
end if

crap2=sqrt(sum((T_s_in-crap1)**2)/float(nx*ny))

patch_hot=crap1+sigma_multiplier*crap2;
patch_cold=crap1-sigma_multiplier*crap2;

print *,'mean, std, patch_hot, patch_cold = ',crap1,crap2,patch_hot,patch_cold
!First do the patch business for temperature
if ((patch_hot .lt. 0._rprec) .OR. (patch_cold .lt. 0._rprec)) then
print *,'Hot & Cold patch calculation yields negative T_s'
print *,'Trying sigma_multiplier = 0.75'
     if (patch_cold < 0._rprec) then
        sigma_multiplier=0.75_rprec
        patch_cold=crap1-sigma_multiplier*crap2;
        patch_hot=crap1+sigma_multiplier*crap2;
        crap3=0.5_rprec*(patch_cold+patch_hot)
        print *,'NEW:mean, patch_hot, patch_cold = ',crap3,patch_hot,patch_cold
        if (patch_cold < 0._rprec) then
           print *,'sigma = 0.75 FAILED, Trying sigma_= 0.5'
           sigma_multiplier=0.5_rprec
           patch_cold=crap1-sigma_multiplier*crap2;
           patch_hot=crap1+sigma_multiplier*crap2;
           crap3=0.5_rprec*(patch_cold+patch_hot)
           print *,'NEW:mean, patch_hot, patch_cold = ',crap3,patch_hot,patch_cold
           if (patch_cold < 0._rprec) then
              print *,'sigma = 0.5 FAILED, STOPPING NOW...'
              print *,'This message is from the subroutine patch_to_remote in scalars_module2.f90.'
           end if
        end if
     end if
end if
! Assign to patches
if (T_or_z .eq. 1) then
  dummy(1:nx/2,:)=patch_hot
  dummy(nx/2+1:nx,:)=patch_cold
   print *,'2 patch T_s field generated: patch_hot,patch_cold = ',patch_hot,patch_cold
elseif (T_or_z .eq. 2) then
! Assign cold patch to the first patch and hot patch to the second one
! Note that this is exactly opposite to what we do for temperature as
! in REALITY, hot surfaces have low roughness and vice-versa
  dummy(1:nx/2,:)=patch_cold
  dummy(nx/2+1:nx,:)=patch_hot
   print *,'2 patch roughnesss field generated: patch_smooth,patch_rough = ',patch_cold,patch_hot
end if

T_s_in=dummy

! testing if remote_to_patch is being used
open(unit=16, file="remote_to_patch_usage.txt",action="write", status="replace")
close(16)


end subroutine remote_to_patch
!!!xxxxxxxxxxxx--------VIJ-------------XXXXXXXXXXXXXXXXXXX

subroutine avgpatch(u_avg)
! computes the averaged value of a variable (at the wall) over a patch
! and assigns it to an nx X ny array

! sc:
! note: this is inefficient: should calculate directly to u_avg --no maybe its good
use sim_param,only:u
implicit none
real(kind=rprec),dimension(nx,ny),intent(inout)::u_avg
integer::i,j,k
real(kind=rprec),dimension(ptypes)::temp

temp=0._rprec
do j=1,ny
do i=1,nx
do k=1,ptypes
   if (patch(i,j).eq.k) then
      temp(patch(i,j))=temp(patch(i,j))+u(i,j,1)/real(patchnum(patch(i,j)))
   end if
end do
end do
end do

do j=1,ny
do i=1,nx
do k=1,ptypes
   if (patch(i,j).eq.k) then
      u_avg(i,j)=real(temp(k))
   end if
end do
end do
end do
end subroutine avgpatch

end module bottombc
