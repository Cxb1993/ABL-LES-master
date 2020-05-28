module scalars_module2
use types,only:rprec
use param !, jt_global => jt  !--rename to avoid name clashes
! could also modify all routines to access jt from param module, not argument list
use bottombc !makes obukhov functions available
use sim_param,only:u,v,w,theta,q,path, &
                   L11t,L22t,L33t,Q11t,Q22t,Q33t, &
                   dudz,dvdz,txx,txz,tyy,tyz,tzz,p
use sgsmodule,only: Nu_t
use scalars_module,only: L,wstar,dTdz,dqdz,sgs_t1,sgs_t2,sgs_t3,sgs_q3,dTdx,dTdy
!use main
 
implicit none

!! --------------------------------------------------------------------
integer,parameter :: average_dim_select_flag=1-(average_dim_num/2) 
! The variable average_dim_select_flag generates the following values based
! on the value of average_dim_num in param.f90 :- 
! a) average_dim_num = 2 : 
!    average_dim_select_flag = 0 ==> average the 3d array over x and y and output the z profile
! b) average_dim_num = 1 : 
!    average_dim_select_flag = 1 ==> average the 3d array over y and output the (x,z) profile
integer, parameter :: dim1_size=average_dim_select_flag*(nx-nz+1)+nz-1
integer, parameter :: dim2_size=average_dim_select_flag*(nz-2)+1
integer, parameter :: dim1_global=average_dim_select_flag*(nx-nz_tot+1)+nz_tot-1
integer, parameter :: dim2_global=average_dim_select_flag*(nz_tot-2)+1
real(kind=rprec):: T_s_min, T_s_max,z_o_min,z_o_max
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
! Part II of the scalar files .. also look at scalars_module.f90
!contains subroutines:
! ic_scal --- initialize velocity fields and scalars
! patch_or_remote - initialize surface boundary conditions
! scalar_in - read surface data from an external data file (remote-sensed)
! append_zeros_string - append zeros in front of a string - called by scalar_in (NOT USED !!)
! scalar_slice - outputs time-averaged x-z slices of scalar variables
! controlled by c_count & p_count from param; Uses file unit numbers from (36-47)
! obukhov_slice - outputs the obukhov variables (phi,psi, L,wstar);
! called from scalar_slice and toggled by parameter obukhov_output (specified in scalars_module.f)
! DIURNAL_FORCING - sets the surface temperature equal to the mean temperature from a data file via a 
! diurnal forcing over time changing the surface temperature at each time step
! 
! Authored by Vijayant Kumar
! Last modified - June 16, 2005
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX

contains

!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
subroutine ic_scal()
!subroutine ic_scal(u,v,w,theta,q,sign_wt_s,theta_mean)
!c...Log profile that is modified to flatten at z=z_i
!c .. Modified by Vijayant to put scalars back
! Last modified April 14, 2004

implicit none
real(kind=rprec),dimension(nz)::ubar,vbar,wbar
real(kind=rprec)::rms, noise, arg, arg2,theta_mean
real(kind=rprec)::z,w_star,T_star,q_star,ran3
real(kind=rprec)::z_turb_limit,perturb_height_factor,z_inv
integer::jx,jy,jz,seed,jz_abs
! SKS
real(kind=rprec),dimension(nz_tot-1)::utotinit,vtotinit,wtotinit,Ttotinit
real(kind=rprec),dimension(nz-1)::uinit,vinit,winit,Tinit
real(kind=rprec)::k_global

$if ($MPI)
  integer::recvcounts(nproc)
  integer::displs(nproc)
$endif
! SKS


! theta_mean=T_s_min-2._rprec !T_s_min is dimensional while T_s is non-dimensional
       if (lbc .eq. 0) then
! we set the mean temp to be T_init (defined in param.f90) for lbc=0 or lbc=1
           theta_mean=T_init
           print *,'theta_mean = ',theta_mean
       else
           theta_mean=T_init
           print *,'theta_mean = ',theta_mean
       end if

       if (wt_s .lt. 0._rprec) then
         ! SKS
         perturb_height_factor=0.50_rprec
         z_inv=0.8_rprec*z_i
         ! perturb_height_factor=0.10_rprec
         ! z_inv=0.10_rprec*z_i
         ! z_inv=0.162_rprec*z_i
         ! SKS
       else
         ! SKS
         perturb_height_factor=0.60_rprec
         ! perturb_height_factor=0.3_rprec
         z_inv=0.8_rprec*z_i
         ! perturb_height_factor=0.3_rprec
         ! SKS
       end if

       z_turb_limit=perturb_height_factor*z_i
      
      if (wt_s .eq. 0.0_rprec) then
! Compute the values of w_star etc. using the default value of
! wt_s = 0.06
      w_star=(g/theta_mean*0.06_rprec*z_i)**(1._rprec/3._rprec)
! w_star is of O(1) with z_i=500 and wt_s=0.06
      T_star=0.06_rprec/w_star
      q_star=T_star
      else
      w_star=sign((g/theta_mean*abs(wt_s)*z_i)**(1._rprec/3._rprec),wt_s)
      T_star=wt_s/w_star
      q_star=T_star
      end if

         $if ($MPI)
            print *,'Modified Log Profile for IC for coord = ',coord
         $else
            print *,'Modified Log Profile for IC'
         $endif
       do jz=1,nz
         $if ($MPI)
            z = (coord*(nz-1) + jz - 0.5_rprec) * dz
         $else
            z = (real(jz)-0.5_rprec)*dz
         $endif
!c IC in equilibrium with rough surface (rough dominates in effective zo)
        arg2=z/(sum(zo)/float(nx*ny))
        arg=(1._rprec/vonk)*log(arg2)!-1./(2.*vonk*z_i*z_i)*z*z

        if (coriolis_forcing) then
!MM ---
         k_global = jz + coord*(nz-1)
       if (k_global > (0.67*nproc*(nz-1))) then
                ubar(jz)=ug_dim/u_star
                vbar(jz)=vg_dim/u_star
       else
               ! ubar(jz)=(k_global/(0.6*nproc*(nz-1)))*ug_dim*0.35/u_star+0.65*ug_dim/u_star    
               ubar(jz)=ug_dim/u_star
               vbar(jz)=vg_dim/u_star       

        end if
       ! vbar(jz)=-(k_global/(nproc*(nz-1)))*5*1.2/u_star
        wbar(jz)=0._rprec
!--- MM
! Note that ug and vg have already been non-dimensionalized in param.f90
        else
        ubar(jz)=arg
        vbar(jz)=0._rprec
        wbar(jz)=0._rprec
        end if
!C sc: I changed around some parenthesis here

        if (z.gt.(z_turb_limit)) then
           ubar(jz)=ubar(jz-1)
        end if
       end do

  rms = 3._rprec
  do jz=1,nz
  $if ($MPI)
    jz_abs = coord * (nz-1) + jz
    z = (coord * (nz-1) + jz - 0.5_rprec) * dz * z_i
  $else
    jz_abs = jz
    z = (jz-.5_rprec) * dz * z_i
  $endif
  seed = -80 - jz_abs  !--trying to make consistent init for MPI
    do jy=1,ny
      do jx=1,nx
!c...Ran3 returns uniform RV between 0 and 1. (sigma_rv=0.289)
!c...Taking std dev of vel as 1 at all heights

!cVK Note that if you put wt_s = 0 symbolizing neutal conditions
!c u should also put L_z=z_i i.e. the inversion layer height should
!c be equal to the height of the domain and in that case the second
!c part of the subsequent if loop will never execute. This is
!c ensured by putting an OR statement in the if clause, which makes 
!c sure that only the first part of if block is executed and not the
!c block after else
       if (z .LE. z_turb_limit) then
         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
         u(jx,jy,jz)=noise*(1._rprec-z/z_i)*w_star/u_star+ubar(jz)
         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
         v(jx,jy,jz)=noise*(1._rprec-z/z_i)*w_star/u_star+vbar(jz) !noise
         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
         w(jx,jy,jz)=noise*(1._rprec-z/z_i)*w_star/u_star+wbar(jz)
         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
         theta(jx,jy,jz)=(theta_mean+10._rprec*noise*(1-z/z_i)*T_star)/T_scale
         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
         q(jx,jy,jz)=q_mix+50._rprec*noise*(1-z/z_i)*q_star
       else
         u(jx,jy,jz)=ubar(jz)
         v(jx,jy,jz)=vbar(jz)
         w(jx,jy,jz)=wbar(jz)
             if ((z .gt. z_turb_limit) .and. (z .le. z_inv)) then
                  theta(jx,jy,jz)=theta_mean/T_scale
             else
                  theta(jx,jy,jz)=(theta_mean+(z-z_inv)*inv_strength)/T_scale
             end if
         q(jx,jy,jz)=q_mix
        end if
       end do
     end do
  end do

  !...BC for W
  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    w(1:nx, 1:ny, 1) = 0._rprec
  end if
  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
    w(1:nx, 1:ny, nz) = 0._rprec
  endif

  !...BC for U, V & T
  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
    u(1:nx, 1:ny, nz) = u(1:nx, 1:ny, nz-1)
    v(1:nx, 1:ny, nz) = v(1:nx, 1:ny, nz-1)
    theta(1:nx, 1:ny, nz) = theta(1:nx, 1:ny, nz-1)+dTdz_top/T_scale*z_i*dz
  end if

!VK Display the mean vertical profiles of the initialized variables on the screen
open(unit=44,file=path//'output/init_profiles.dat',status="unknown",position="append")

do jz=1,nz
     $if ($MPI)
       z = (coord*(nz-1) + jz - 0.5_rprec) * dz
     $else
       z = (jz - 0.5_rprec) * dz
     $endif
! SKS
     write(6,7781) jz,z,(sum(u(:,:,jz))/float(nx*ny)),(sum(v(:,:,jz))/&
     float(nx*ny)),(sum(w(:,:,jz))/float(nx*ny)),&
     (sum(theta(:,:,jz))/float(nx*ny))*T_scale
!    write(6,7781) jz,z,(sum(u(:,:,jz))/float(nx*ny))*u_star,(sum(v(:,:,jz))/&
!    float(nx*ny))*u_star,(sum(w(:,:,jz))/float(nx*ny))*u_star,&
!    (sum(theta(:,:,jz))/float(nx*ny))*T_scale
! SKS
     write(44,7781) jz,z,(sum(u(:,:,jz))/float(nx*ny))*u_star,(sum(v(:,:,jz))/&
     float(nx*ny))*u_star,(sum(w(:,:,jz))/float(nx*ny))*u_star,&
     (sum(theta(:,:,jz))/float(nx*ny))*T_scale
end do
close(44)
7781 format('jz, z, ubar, vbar, wbar,T_bar:',(1x,I3,1x,F9.4,1x,F9.4,1x,F9.4,1x,F9.4,1x,F9.4))

! SKS
! This part written just to get the initial profiles 

open(unit=103,file=path//'output/initial_profiles.dat',status="unknown",position="append")

do jz=1,nz-1      
   uinit(jz) = (sum(u(:,:,jz))/float(nx*ny))*u_star
   vinit(jz) = (sum(v(:,:,jz))/float(nx*ny))*u_star
   winit(jz) = (sum(w(:,:,jz))/float(nx*ny))*u_star
   Tinit(jz) = (sum(theta(:,:,jz))/float(nx*ny))*T_scale
end do
      
$if ($MPI) 
  recvcounts = size(uinit)
  displs = coord_of_rank * recvcounts 
  call mpi_gatherv (uinit(1),size(uinit),MPI_RPREC, &
                    utotinit(1),recvcounts,displs,  &
                    MPI_RPREC,rank_of_coord(0),comm,ierr)
  recvcounts = size(vinit)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (vinit(1),size(vinit),MPI_RPREC, &
                    vtotinit(1),recvcounts,displs,  &
                    MPI_RPREC,rank_of_coord(0),comm,ierr)
  recvcounts = size(winit)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (winit(1),size(winit),MPI_RPREC, &
                    wtotinit(1),recvcounts,displs,  &
                    MPI_RPREC,rank_of_coord(0),comm,ierr)
  recvcounts = size(Tinit)
  displs = coord_of_rank * recvcounts 
  call mpi_gatherv (Tinit(1),size(Tinit),MPI_RPREC, &
                    Ttotinit(1),recvcounts,displs,  &
                    MPI_RPREC,rank_of_coord(0),comm,ierr)
$else
  utotinit = uinit
  vtotinit = vinit
  wtotinit = winit
  Ttotinit = Tinit
$endif

if((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  do jz=1,nz_tot-1
     write(103,8001) utotinit(jz),vtotinit(jz),wtotinit(jz),Ttotinit(jz)
  end do
  write(103,*)
end if
8001  format(1400(E14.5))
! SKS    

end subroutine ic_scal

!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
subroutine patch_or_remote()
implicit none
real(kind=rprec):: T_s_init_mean,z_o_init_mean
integer :: i,jy
! z_os already defined in scalars_module
! zo,T_s and q_s defined in bottombc
! April 20, 2004 - so far contains both patch_or_remote
! and scalar_in subroutine

if (patch_flag .eq. 1) then
print *, 'Assigning values to patches'
call patches()

z_os(:,:)=(1._rprec/10._rprec)*zo(:,:)

! Calculate minimum and maximum values of T_s and zo for use in
! ic_scal()
   if (lbc .eq. 0) then
   T_s_min=minval(T_s); T_s_max=maxval(T_s)   
   z_o_min=minval(zo); z_o_max=maxval(zo)
   print *,'Min and Max (T_s,zo)',T_s_min,T_s_max,z_o_min,z_o_max
   end if
!c sets temperature field and roughness field at the bottom , x-y plane
!c Added by Vijayant
else if (remote_flag .eq. 1) then
print *, 'Assigning remote-sensed values to the surface'
   call scalar_in() ! Updated T_s and zo loaded from bottombc
   
! Calculate minimum and maximum values of T_s and zo for use in ic_scal()
   if (lbc .eq. 0) then
   T_s_min=minval(T_s); T_s_max=maxval(T_s)   
   z_o_min=minval(zo); z_o_max=maxval(zo)
   print *,'Min and Max (T_s,zo)',T_s_min,T_s_max,z_o_min,z_o_max
   end if

   T_s_init_mean=sum(T_s)/float(nx*ny)
!  z_o_init_mean=sum(zo)/float(nx*ny)
!  perform logarithmic average for zo
   z_o_init_mean=exp(sum(dlog(zo(1:nx,1:ny)))/float(nx*ny))

if (remote_to_patch_flag) then
   call remote_to_patch(T_s,1)
   call remote_to_patch(zo,2)
end if

if (remote_homog_flag .eq. 1) then
print *,'Homogeinizing the remote b.c.s'
   T_s=T_s_init_mean
   zo=z_o_init_mean
end if

! Non-dimensionalize T
   T_s=(T_s)/T_scale
   zo=zo/z_i

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! z_os is the scalar roughness length. Divide momentum
! roughness length by 10 to get scalar roughness length
! data (look in the scalar_in routine above)
z_os(:,:)=(1._rprec/10._rprec)*zo(:,:)

      if (GABLS_test) then
        zo(:,:)=0.1_rprec/z_i
        z_os(:,:)=zo(:,:)
        print *,'setting zo and zo_s for GABLS case = 0.1 m'
      end if
end if

! Write surface bc input to a file
open(unit=57,file=path//'output/surface_bc_input.txt',status="unknown",position="append")
do jy=1,ny
write(57,5168) (T_s(i,jy),i=1,nx)
end do
do jy=1,ny
write(57,5168) (zo(i,jy),i=1,nx)
end do
close(57)

5168     format(1400(E14.5))

end subroutine patch_or_remote

!!!xxxxxxxx--------VIJ---------XXXXXXXXXXXX-------------
!!!xxxxxxxx--------VIJ---------XXXXXXXXXXXX-------------
subroutine scalar_in()
!c This reads in the scalar input from a interpolated scalar
!c file and assigns it to the x-y plane at z=0 i.e. at the ground
!c This for the moment is used to read the momentum roughness, z0
!c and temperature from the USDA remote-sensed data set. The data 
!c can be interpolated using bilinear/cubic/spline interpolation (MATLAB)
!c for the grid size(nx,ny)
!c Authored by Vijayant Kumar
!c Last modified on April 11th, 2004

use bottombc,only:T_s,zo !Load the variables from bottombc and update in here
implicit none
integer:: ii,jj
character(len=6):: suffix2,suffix3
      
       write(suffix2,'(i6.6)') nx ! create a string from nx   
       
      if (coarse_grain_flag) then
      write(suffix3,'(i6.6)') stencil_pts    
 
open(unit=77,file='../interp_data/coarse_grained/interp_temp_cg_'&
//suffix2(4:6)//'pts_'//suffix3(4:6)//'.out',status='unknown')
open(unit=78,file='../interp_data/coarse_grained/interp_z_m_cg_'&
//suffix2(4:6)//'pts_'//suffix3(4:6)//'.out',status='unknown')
     
print *,'interp_temp_cg_'//suffix2(4:6)//'pts_'&
//suffix3(4:6)//'.out loaded from scalar_in.f'
       
      do jj=1,ny
          read(77,5169) (T_s(ii,jj),ii=1,nx)
          read(78,5169) (zo(ii,jj),ii=1,nx)
      enddo
      T_s=T_s+273.15_rprec ! Convert T to Kelvin
      close(77)
      close(78)

       else

open(unit=77,file=path//'interp_data/interp_temp_'//suffix2(4:6)&
//'X'//suffix2(4:6)//'_cubic.out',status='unknown')
open(unit=78,file=path//'interp_data/interp_z_m_'//suffix2(4:6)&
//'X'//suffix2(4:6)//'_cubic.out',status='unknown')
print *,path//'interp_data/interp_temp_'//suffix2(4:6)//'X'//suffix2(4:6)//'_cubic.out'
       
      do jj=1,ny
          read(77,5169) (T_s(ii,jj),ii=1,nx)
          read(78,5169) (zo(ii,jj),ii=1,nx)
      enddo
      print *,'Mean T_s (K) = ',sum(T_s)/float(nx*ny)
      print *,'Mean zo (m) = ',sum(zo)/float(nx*ny)
      close(77)
      close(78)
          T_s=T_s+273.15_rprec ! Convert T to Kelvin
      end if
5169     format(1400(E16.10))
end subroutine scalar_in

!!!xxxxxxxxxxxx--------VIJ-------------XXXXXXXXXXXXXXXXXXX
!!!xxxxxxxxxx-------scalar output subroutine-----XXXXXXXXXXXXXXXXXXXXX

subroutine scalar_slice()
!c This is exactly the same like the subroutine avgslice with the
!c only difference being that it averages the scalar variables
!c to find the y-averaged instantaneous x-z slices of variables
!c t,q,sgs_t3,sgs_q3 and their variances such as t2, q2.
!c It also outputs the average covariance between wt and wq
!use scalars_module,only: dTdz,dqdz,sgs_t3,sgs_q3 
!use output_slice,only: collocate_MPI_averages
use sim_param,only:path,u,v,w,dudz,dvdz,txx,txz,tyy,tyz,tzz,p,txy,&
         dudx,dudy,dvdx,dvdy,dwdx,dwdy,dwdz
implicit none
integer:: i,j,k
real(kind=rprec),dimension(nx,nz-1),save:: atheta,t2,asgs_t1,asgs_t2,asgs_t3,awt,aut,avt,awwth,awuth,awvth
real(kind=rprec),dimension(nx,nz-1),save:: apthet,asigthet,apdtdx,apdtdy,apdtdz,asigdtdx,asigdtdy,asigdtdz
real(kind=rprec),dimension(nx,nz-1),save:: awpie3,aupie3,avpie3,apieux1,apieuy2,apieuz3,apiewx1,apiewy2,apiewz3
real(kind=rprec),dimension(nx,nz-1),save:: apievx1,apievy2,apievz3,at33t,at13t,at23t,at31tx,at32ty,at33tz
real(kind=rprec),dimension(nx,nz-1),save:: at11tx,at12ty,at13tz,at21tx,at22ty,at23tz
real(kind=rprec),dimension(nx,nz-1),save:: adTdz,adTdz2,adTdx,adTdy,anu_t,t3
real(kind=rprec):: ttheta1,tt2,tsgst1,tsgst2,tsgst3,twt,tut,tvt,twwth,twuth,twvth
real(kind=rprec)::tpthet,tsigthet,tpdtdx,tpdtdy,tpdtdz,tsigdtdx,tsigdtdy,tsigdtdz
real(kind=rprec)::twpie3,tupie3,tvpie3,tpieux1,tpieuy2,tpieuz3,tpiewx1,tpiewy2,tpiewz3
real(kind=rprec)::tpievx1,tpievy2,tpievz3,tt33t,tt13t,tt23t,tt31tx,tt32ty,tt33tz
real(kind=rprec)::tt11tx,tt12ty,tt13tz,tt21tx,tt22ty,tt23tz,sig1,sig1q
real(kind=rprec)::tdTdz,tdTdz2,tdTdx,tdTdy,tnu_t,tt3,arg1,arg2,arg3,arg3f,arg5,arg6
real(kind=rprec)::fr,ux,uy,uz,vx,vy,vz,wx,wy,wz,tx,ty,tz,arg4,arg7,arg8,arg9,arg10
real(kind=rprec),dimension(:,:),allocatable:: avg_scalar_out
real(rprec),parameter::delbar=2._rprec
real(rprec)::denomsig=0.5_rprec/((delbar**(2.0_rprec/3.0_rprec))-1._rprec)
real(rprec)::denomsigq=0.5_rprec/(((2._rprec*delbar)**(2.0_rprec/3.0_rprec))-1._rprec)

fr=(1._rprec/float(p_count))*float(c_count)

!if (jt .EQ. c_count) then
!  atheta=0._rprec;t2=0._rprec;asgs_t3=0._rprec;awt=0._rprec;adTdz=0._rprec
!  anu_t=0._rprec;t3=0._rprec;var_t=0._rprec
!end if

do k=1,nz-1
do i=1,nx
 ttheta1=0._rprec;tt2=0._rprec;tsgst1=0._rprec;tsgst2=0._rprec;tsgst3=0._rprec;
 twt=0._rprec;tut=0._rprec;tvt=0._rprec;twwth=0._rprec;twuth=0._rprec;twvth=0._rprec;
 tpthet=0._rprec;tsigthet=0._rprec;tpdtdx=0._rprec;tpdtdy=0._rprec;tpdtdz=0._rprec;
 tsigdtdx=0._rprec;tsigdtdy=0._rprec;tsigdtdz=0._rprec;
 twpie3=0._rprec;tupie3=0._rprec;tvpie3=0._rprec;tpieux1=0._rprec;tpieuy2=0._rprec;
 tpieuz3=0._rprec;tpiewx1=0._rprec;tpiewy2=0._rprec;tpiewz3=0._rprec;
 tpievx1=0._rprec;tpievy2=0._rprec;tpievz3=0._rprec;tt33t=0._rprec;tt13t=0._rprec;
 tt23t=0._rprec;tt31tx=0._rprec;tt32ty=0._rprec;tt33tz=0._rprec;
 tt11tx=0._rprec;tt12ty=0._rprec;tt13tz=0._rprec;tt21tx=0._rprec;tt22ty=0._rprec;
 tt23tz=0._rprec;
 tdTdz=0._rprec;tdTdx=0._rprec;tdTdy=0._rprec;tnu_t=0._rprec;tt3=0._rprec;

do j=1,ny

    if((k .eq. 1) .AND. ((.not. USE_MPI) .or. (USE_MPI .and. coord .eq. 0))) then

         arg1=u(i,j,k)
         arg2=v(i,j,k)
         arg4=0.25_rprec*w(i,j,k+1)
         arg3=p(i,j,k)-0.5_rprec*(arg1*arg1+arg2*arg2+arg4*arg4)
! KMG: Lii and Qii are already on cs nodes 
         arg3f= p(i,j,k)- 0.5_rprec*(arg1*arg1+arg2*arg2+arg4*arg4)-(1._rprec/3._rprec)*(denomsig)*(L11t(i,j,k)+L22t(i,j,k)+L33t(i,j,k))
!         arg5=txx(i,j,k)
!         arg6=txy(i,j,k)
!         arg7=tyy(i,j,k)
!         arg8=tzz(i,j,k)
         arg9=(txz(i,j,k)+txz(i,j,k+1))/2.0_rprec
         arg10=(tyz(i,j,k)+tyz(i,j,k+1))/2.0_rprec
        
         ux=dudx(i,j,k)
         uy=dudy(i,j,k)
         vx=dvdx(i,j,k)
         vy=dvdy(i,j,k)
         wz=dwdz(i,j,k)
         wx=(dwdx(i,j,k)+dwdx(i,j,k+1))/2.0_rprec
         wy=(dwdy(i,j,k)+dwdy(i,j,k+1))/2.0_rprec
         tx=dTdx(i,j,k)
         ty=dTdy(i,j,k)
         else

         arg1=(u(i,j,k)+u(i,j,k-1))/2.0_rprec
         arg2=(v(i,j,k)+v(i,j,k-1))/2.0_rprec
         arg4=w(i,j,k)
         arg3=(p(i,j,k)+p(i,j,k-1))/2.0_rprec - 0.5_rprec*(arg1*arg1+arg2*arg2+arg4*arg4)
         arg3f=(p(i,j,k)+p(i,j,k-1))/2.0_rprec - 0.5_rprec*(arg1*arg1+arg2*arg2+arg4*arg4)-(1._rprec/3._rprec)*(denomsig)*(L11t(i,j,k)+L22t(i,j,k)+L33t(i,j,k))
!         arg5=(txx(i,j,k)+txx(i,j,k-1))/2.0_rprec
!         arg6=(txy(i,j,k)+txy(i,j,k-1))/2.0_rprec
!         arg7=(tyy(i,j,k)+tyy(i,j,k-1))/2.0_rprec
!         arg8=(tzz(i,j,k)+tzz(i,j,k-1))/2.0_rprec
         arg9=txz(i,j,k)
         arg10=tyz(i,j,k)
 
         ux=(dudx(i,j,k)+dudx(i,j,k-1))/2.0_rprec
         uy=(dudy(i,j,k)+dudy(i,j,k-1))/2.0_rprec
         vx=(dvdx(i,j,k)+dvdx(i,j,k-1))/2.0_rprec
         vy=(dvdy(i,j,k)+dvdy(i,j,k-1))/2.0_rprec
         wz=(dwdz(i,j,k)+dwdz(i,j,k-1))/2.0_rprec
         wx=dwdx(i,j,k)
         wy=dwdy(i,j,k)
         tx=(dTdx(i,j,k)+dTdx(i,j,k-1))/2.0_rprec
         ty=(dTdy(i,j,k)+dTdy(i,j,k-1))/2.0_rprec

        end if
      uz=dudz(i,j,k)
      vz=dvdz(i,j,k)
      tz=dTdz(i,j,k)
      arg5=txx(i,j,k)
      arg6=txy(i,j,k)
      arg7=tyy(i,j,k)
      arg8=tzz(i,j,k)
      sig1=denomsig*(L11t(i,j,k)+L22t(i,j,k)+L33t(i,j,k))
      sig1q=denomsigq*(Q11t(i,j,k)+Q22t(i,j,k)+Q33t(i,j,k))
 
    ttheta1=ttheta1+theta(i,j,k)
    tt2=tt2+theta(i,j,k)*theta(i,j,k)
    tsgst1=tsgst1+sgs_t1(i,j,k)
    tsgst2=tsgst2+sgs_t2(i,j,k)
    tsgst3=tsgst3+sgs_t3(i,j,k)
    tdTdz=tdTdz+dTdz(i,j,k)
  !  tdTdz2=tdTdz2+tz
    tdTdx=tdTdx+dTdx(i,j,k)
    tdTdy=tdTdy+dTdy(i,j,k)
    tnu_t=tnu_t+Nu_t(i,j,k)
    twt=twt+arg4*theta(i,j,k)
    tut=tut+arg1*theta(i,j,k)
    tvt=tvt+arg2*theta(i,j,k)
    twwth=twwth+arg4*arg4*theta(i,j,k)
    twuth=twuth+arg4*arg1*theta(i,j,k)
    twvth=twvth+arg4*arg2*theta(i,j,k) 
    tt3=tt3+theta(i,j,k)*theta(i,j,k)*theta(i,j,k)
    tpthet=tpthet+theta(i,j,k)*arg3f
    tpdtdx=tpdtdx+tx*arg3f
    tpdtdy=tpdtdy+ty*arg3f
    tpdtdz=tpdtdz+tz*arg3f
    tsigdtdx=tsigdtdx+tx*sig1
    tsigdtdy=tsigdtdy+ty*sig1
    tsigdtdz=tsigdtdz+tz*sig1
    twpie3=twpie3+arg4*sgs_t3(i,j,k)
    tupie3=tupie3+arg1*sgs_t3(i,j,k)
    tvpie3=tvpie3+arg2*sgs_t3(i,j,k)
    tpieux1=tpieux1+ux*sgs_t1(i,j,k)
    tpieuy2=tpieuy2+uy*sgs_t2(i,j,k)
    tpieuz3=tpieuz3+uz*sgs_t3(i,j,k)
    tpiewx1=tpiewx1+wx*sgs_t1(i,j,k)
    tpiewy2=tpiewy2+wy*sgs_t2(i,j,k)
    tpiewz3=tpiewz3+wz*sgs_t3(i,j,k)
    tpievx1=tpievx1+vx*sgs_t1(i,j,k)
    tpievy2=tpievy2+vy*sgs_t2(i,j,k)
    tpievz3=tpievz3+vz*sgs_t3(i,j,k)
    tt33t=tt33t+theta(i,j,k)*arg8
    tt13t=tt13t+theta(i,j,k)*arg9
    tt23t=tt23t+theta(i,j,k)*arg10
    tt31tx=tt31tx+tx*arg9
    tt33tz=tt33tz+tz*arg8
    tt32ty=tt32ty+ty*arg10
    tt11tx=tt11tx+tx*arg5
    tt12ty=tt12ty+ty*arg6
    tt13tz=tt13tz+tz*arg9
    tt21tx=tt21tx+tx*arg6
    tt22ty=tt22ty+ty*arg7
    tt23tz=tt23tz+tz*arg10

end do

   atheta(i,k)=atheta(i,k)+(fr)*ttheta1/Ny
   t2(i,k)=t2(i,k)+(fr)*tt2/Ny
   asgs_t1(i,k)=asgs_t1(i,k)+(fr)*tsgst1/Ny
   asgs_t2(i,k)=asgs_t2(i,k)+(fr)*tsgst2/Ny
   asgs_t3(i,k)=asgs_t3(i,k)+(fr)*tsgst3/Ny
   awt(i,k)=awt(i,k)+(fr)*twt/Ny
   aut(i,k)=aut(i,k)+(fr)*tut/Ny
   avt(i,k)=avt(i,k)+(fr)*tvt/Ny
   awwth(i,k)=awwth(i,k)+(fr)*twwth/Ny
   awuth(i,k)=awuth(i,k)+(fr)*twuth/Ny
   awvth(i,k)=awvth(i,k)+(fr)*twvth/Ny
   apthet(i,k)=apthet(i,k)+(fr)*tpthet/Ny
   asigthet(i,k)=asigthet(i,k)+(fr)*tsigthet/Ny
   apdtdx(i,k)=apdtdx(i,k)+(fr)*tpdtdx/Ny
   apdtdy(i,k)=apdtdy(i,k)+(fr)*tpdtdy/Ny
   apdtdz(i,k)=apdtdz(i,k)+(fr)*tpdtdz/Ny

   asigdtdx(i,k)=asigdtdx(i,k)+(fr)*tsigdtdx/Ny
   asigdtdy(i,k)=asigdtdy(i,k)+(fr)*tsigdtdy/Ny
   asigdtdz(i,k)=asigdtdz(i,k)+(fr)*tsigdtdz/Ny

   awpie3(i,k)=awpie3(i,k)+(fr)*twpie3/Ny
   aupie3(i,k)=aupie3(i,k)+(fr)*tupie3/Ny 
   avpie3(i,k)=avpie3(i,k)+(fr)*tvpie3/Ny
 
   apieux1(i,k)=apieux1(i,k)+(fr)*tpieux1/Ny
   apieuy2(i,k)=apieuy2(i,k)+(fr)*tpieuy2/Ny
   apieuz3(i,k)=apieuz3(i,k)+(fr)*tpieuz3/Ny
   apiewx1(i,k)=apiewx1(i,k)+(fr)*tpiewx1/Ny
   apiewy2(i,k)=apiewy2(i,k)+(fr)*tpiewy2/Ny
   apiewz3(i,k)=apiewz3(i,k)+(fr)*tpiewz3/Ny
   apievx1(i,k)=apievx1(i,k)+(fr)*tpievx1/Ny
   apievy2(i,k)=apievy2(i,k)+(fr)*tpievy2/Ny
   apievz3(i,k)=apievz3(i,k)+(fr)*tpievz3/Ny
   at33t(i,k)=at33t(i,k)+(fr)*tt33t/Ny
   at13t(i,k)=at13t(i,k)+(fr)*tt13t/Ny
   at23t(i,k)=at23t(i,k)+(fr)*tt23t/Ny
   at31tx(i,k)=at31tx(i,k)+(fr)*tt31tx/Ny
   at32ty(i,k)=at32ty(i,k)+(fr)*tt32ty/Ny
   at33tz(i,k)=at33tz(i,k)+(fr)*tt33tz/Ny
   at11tx(i,k)=at11tx(i,k)+(fr)*tt11tx/Ny
   at12ty(i,k)=at12ty(i,k)+(fr)*tt12ty/Ny
   at13tz(i,k)=at13tz(i,k)+(fr)*tt13tz/Ny
   at21tx(i,k)=at21tx(i,k)+(fr)*tt21tx/Ny
   at22ty(i,k)=at22ty(i,k)+(fr)*tt22ty/Ny
   at23tz(i,k)=at23tz(i,k)+(fr)*tt23tz/Ny
   adTdz(i,k)=adTdz(i,k)+(fr)*tdTdz/Ny
  ! adTdz2(i,k)=adTdz2(i,k)+(fr)*tdTdz2/Ny
   adTdx(i,k)=adTdx(i,k)+(fr)*tdTdx/Ny
   adTdy(i,k)=adTdy(i,k)+(fr)*tdTdy/Ny
   anu_t(i,k)=anu_t(i,k)+(fr)*tnu_t/Ny
   t3(i,k)=t3(i,k)+(fr)*tt3/Ny

end do
end do
      
if (mod(jt,p_count)==0) then
  allocate(avg_scalar_out(1:nx,1:nz_tot-1));
  call collocate_MPI_averages_N(atheta,avg_scalar_out,35,'theta')
  call collocate_MPI_averages_N(t2,avg_scalar_out,36,'t2')
  call collocate_MPI_averages_N(asgs_t1,avg_scalar_out,37,'sgs_t1')
  call collocate_MPI_averages_N(asgs_t2,avg_scalar_out,38,'sgs_t2')
  call collocate_MPI_averages_N(asgs_t3,avg_scalar_out,39,'sgs_t3')
  call collocate_MPI_averages_N(awt,avg_scalar_out,40,'wt')
  call collocate_MPI_averages_N(aut,avg_scalar_out,41,'ut')
  call collocate_MPI_averages_N(avt,avg_scalar_out,42,'vt')
  call collocate_MPI_averages_N(awwth,avg_scalar_out,43,'wwt')
  call collocate_MPI_averages_N(awuth,avg_scalar_out,44,'wut')
  call collocate_MPI_averages_N(awvth,avg_scalar_out,45,'wvt')
  call collocate_MPI_averages_N(apthet,avg_scalar_out,46,'ptheta')
  call collocate_MPI_averages_N(asigthet,avg_scalar_out,47,'sigtheta')
  call collocate_MPI_averages_N(apdtdx,avg_scalar_out,48,'pdtdx')
  call collocate_MPI_averages_N(apdtdy,avg_scalar_out,49,'pdtdy')
  call collocate_MPI_averages_N(apdtdz,avg_scalar_out,50,'pdtdz')
  call collocate_MPI_averages_N(asigdtdx,avg_scalar_out,51,'sigdtdx')
  call collocate_MPI_averages_N(asigdtdy,avg_scalar_out,52,'sigdtdy')
  call collocate_MPI_averages_N(asigdtdz,avg_scalar_out,53,'sigdtdz')
  call collocate_MPI_averages_N(awpie3,avg_scalar_out,54,'wpie3')
  call collocate_MPI_averages_N(aupie3,avg_scalar_out,55,'upie3')
  call collocate_MPI_averages_N(avpie3,avg_scalar_out,56,'vpie3')
  call collocate_MPI_averages_N(apieux1,avg_scalar_out,57,'pie1ux')
  call collocate_MPI_averages_N(apieuy2,avg_scalar_out,58,'pie2uy')
  call collocate_MPI_averages_N(apieuz3,avg_scalar_out,59,'pie3uz')
  call collocate_MPI_averages_N(apiewx1,avg_scalar_out,60,'pie1wx')
  call collocate_MPI_averages_N(apiewy2,avg_scalar_out,61,'pie2wy')
  call collocate_MPI_averages_N(apiewz3,avg_scalar_out,62,'pie3wz')

  call collocate_MPI_averages_N(apievx1,avg_scalar_out,63,'pie1vx')
  call collocate_MPI_averages_N(apievy2,avg_scalar_out,64,'pie2vy')
  call collocate_MPI_averages_N(apievz3,avg_scalar_out,65,'pie3vz')

  call collocate_MPI_averages_N(at33t,avg_scalar_out,66,'t33t')
  call collocate_MPI_averages_N(at13t,avg_scalar_out,67,'t13t')
  call collocate_MPI_averages_N(at23t,avg_scalar_out,68,'t23t')


  call collocate_MPI_averages_N(at31tx,avg_scalar_out,69,'t31tx')
  call collocate_MPI_averages_N(at32ty,avg_scalar_out,70,'t32ty')
  call collocate_MPI_averages_N(at33tz,avg_scalar_out,71,'t33tz')

  call collocate_MPI_averages_N(at11tx,avg_scalar_out,72,'t11tx')
  call collocate_MPI_averages_N(at12ty,avg_scalar_out,73,'t12ty')
  call collocate_MPI_averages_N(at13tz,avg_scalar_out,74,'t13tz')

  call collocate_MPI_averages_N(at21tx,avg_scalar_out,301,'t21tx')
  call collocate_MPI_averages_N(at22ty,avg_scalar_out,302,'t22ty')
  call collocate_MPI_averages_N(at23tz,avg_scalar_out,304,'t23tz')

  call collocate_MPI_averages_N(adTdz,avg_scalar_out,305,'dTdz')
 ! call collocate_MPI_averages_N(adTdz2,avg_scalar_out,306,'dTdz2')
  call collocate_MPI_averages_N(adTdx,avg_scalar_out,307,'dTdx')
  call collocate_MPI_averages_N(adTdy,avg_scalar_out,308,'dTdy') 
  call collocate_MPI_averages_N(anu_t,avg_scalar_out,309,'Nu_t')
  call collocate_MPI_averages_N(t3,avg_scalar_out,310,'t3')
  deallocate(avg_scalar_out)

 atheta=0._rprec;t2=0._rprec;asgs_t1=0._rprec;asgs_t2=0._rprec;asgs_t3=0._rprec;
 awt=0._rprec;aut=0._rprec;avt=0._rprec;awwth=0._rprec;awuth=0._rprec;awvth=0._rprec;
 apthet=0._rprec;asigthet=0._rprec;apdtdx=0._rprec;apdtdy=0._rprec;apdtdz=0._rprec;
 asigdtdx=0._rprec;asigdtdy=0._rprec;asigdtdz=0._rprec;
 awpie3=0._rprec;aupie3=0._rprec;avpie3=0._rprec;apieux1=0._rprec;apieuy2=0._rprec;
 apieuz3=0._rprec;apiewx1=0._rprec;apiewy2=0._rprec;apiewz3=0._rprec;
 apievx1=0._rprec;apievy2=0._rprec;apievz3=0._rprec;at33t=0._rprec;at13t=0._rprec;
 at23t=0._rprec;at31tx=0._rprec;at32ty=0._rprec;at33tz=0._rprec;
 at11tx=0._rprec;at12ty=0._rprec;at13tz=0._rprec;at21tx=0._rprec;at22ty=0._rprec;
 at23tz=0._rprec;
 adTdz=0._rprec;adTdz2=0._rprec;adTdx=0._rprec;adTdy=0._rprec;anu_t=0._rprec;t3=0._rprec;


end if

5168      format(1400(E14.5))

end subroutine scalar_slice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--The following subroutine does the collocation of the MPI arrays for
! averaging in avgslice and scalar_slice (in scalars_module2.f90)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine collocate_MPI_averages_N(avg_var_proc,avg_var_tot_domain,file_ind,filename_str)
use param
$if ($MPI)
  integer :: recvcounts(nproc)
  integer :: displs(nproc)
$endif
integer :: ind1,ind2,file_ind

character (*),intent(in) :: filename_str
character (len=256) :: local_filename
real(kind=rprec),dimension(dim1_size,dim2_size)::avg_var_proc
real(kind=rprec),dimension(dim1_global,dim2_global)::avg_var_tot_domain

local_filename=path//'output/aver_'//trim(filename_str)//'.out'

  avg_var_tot_domain=0._rprec
$if ($MPI)
  recvcounts = size (avg_var_proc)
  displs = coord_of_rank * recvcounts 
  call mpi_gatherv (avg_var_proc(1,1), size (avg_var_proc), MPI_RPREC,                &
                    avg_var_tot_domain(1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)
$else
  avg_var_tot_domain=avg_var_proc
$endif

  if((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then  
  open(file_ind,file=trim(local_filename),status="unknown",position="append")
      if (average_dim_num .eq. 1) then
         do ind2=1,nz_tot-1
          write(file_ind,5168) jt*dt,(avg_var_tot_domain(ind1,ind2),ind1=1,nx)
         end do
      else if (average_dim_num .eq. 2) then
         write(file_ind,5168) jt*dt,(avg_var_tot_domain(ind1,1),ind1=1,nz_tot-1)
      end if
  close(file_ind)
  end if

 5168     format(1400(E14.5))
end subroutine collocate_MPI_averages_N
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--The following subroutine does the collocation of the MPI arrays for
! averaging in avgslice and scalar_slice (in scalars_module2.f90)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine collocate_MPI_averages(avg_var_proc,avg_var_tot_domain,file_ind)
use param
$if ($MPI)
  integer :: recvcounts(nproc)
  integer :: displs(nproc)
$endif
integer :: ind1,ind2,file_ind
real(kind=rprec),dimension(dim1_size,dim2_size)::avg_var_proc
real(kind=rprec),dimension(dim1_global,dim2_global)::avg_var_tot_domain

  avg_var_tot_domain=0._rprec
  $if ($MPI)
  recvcounts = size (avg_var_proc)
  displs = coord_of_rank * recvcounts 
  call mpi_gatherv (avg_var_proc(1,1), size (avg_var_proc), MPI_RPREC,                &
                    avg_var_tot_domain(1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)
$else
  avg_var_tot_domain=avg_var_proc
$endif

  if((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then  
      if (average_dim_num .eq. 1) then
         do ind2=1,nz_tot-1
          write(file_ind,5168) jt*dt,(avg_var_tot_domain(ind1,ind2),ind1=1,nx)
         end do
      else if (average_dim_num .eq. 2) then
         write(file_ind,5168) jt*dt,(avg_var_tot_domain(ind1,1),ind1=1,nz_tot-1)
      end if
         close(file_ind)
  end if

5168     format(1400(E14.5))
end subroutine collocate_MPI_averages
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
subroutine timestep_conditions(CFL,visc_stab)
! This subroutine computes CFl and viscous stability and is called every wbase timesteps from main.f90
implicit none
$if ($MPI)
  integer :: recvcounts(nproc)
  integer :: displs(nproc)
! SKS
  real(kind=rprec),dimension(4,nproc)::max_var_tot_domain
! real(kind=rprec),dimension(nproc,4)::max_var_tot_domain
! SKS
$endif

real(kind=rprec) :: delta, u_res_max, nu_max
real(kind=rprec),dimension(1,4) :: max_vels
real(kind=rprec),intent(out) :: CFL, visc_stab
 
$if ($MPI)
  recvcounts = size (max_vels)
  displs = coord_of_rank * recvcounts 
  max_vels(1,1)=maxval(u(1:nx,1:ny,1:nz-1))
  max_vels(1,2)=maxval(v(1:nx,1:ny,1:nz-1))
  max_vels(1,3)=maxval(w(1:nx,1:ny,1:nz-1))
  max_vels(1,4)=maxval(Nu_t(1:nx,1:ny,1:nz-1))
  call mpi_gatherv (max_vels(1,1), size(max_vels), MPI_RPREC,                &
                    max_var_tot_domain(1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)

  ! SKS
  u_res_max=sqrt(maxval(max_var_tot_domain(1,:))**2+maxval(max_var_tot_domain(2,:))**2+&
  maxval(max_var_tot_domain(3,:))**2)
  ! u_res_max=sqrt(maxval(max_var_tot_domain(:,1))**2+maxval(max_var_tot_domain(:,2))**2+&
  ! maxval(max_var_tot_domain(:,3))**2)
  ! SKS
  nu_max=maxval(max_var_tot_domain(:,4))
$else
  u_res_max = sqrt(maxval(u(1:nx,1:ny,1:nz-1)**2+v(1:nx,1:ny,1:nz-1)**2+&
  w(1:nx,1:ny,1:nz-1)**2)) ! don't bother with interp here
  nu_max=maxval(Nu_t(1:nx,1:ny,1:nz-1))
$endif  

if ((.not. USE_MPI) .or. (USE_MPI .and. rank == 0)) then
    delta = min(dx, dy, dz)
    CFL = u_res_max*dt/delta
    visc_stab = dt*nu_max/delta**2
end if

end subroutine timestep_conditions

!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
subroutine ic_scal_GABLS()
!subroutine ic_scal(u,v,w,theta,q,sign_wt_s,theta_mean)
!c...Log profile that is modified to flatten at z=z_i
!c .. Modified by Vijayant to put scalars back
! Last modified April 14, 2004

implicit none
real(kind=rprec),dimension(nz)::ubar,vbar,wbar
real(kind=rprec)::rms, noise, arg, arg2,theta_mean
real(kind=rprec)::z,w_star,T_star,q_star,ran3,perturb_height_factor
integer::jx,jy,jz,seed,jz_abs

       theta_mean=T_init
       perturb_height_factor=50._rprec/z_i
      
      if (wt_s .eq. 0.0_rprec) then
! Compute the values of w_star etc. using the default value of
! wt_s = 0.06
      w_star=(g/theta_mean*0.06_rprec*z_i)**(1._rprec/3._rprec)
! w_star is of O(1) with z_i=500 and wt_s=0.06
      T_star=0.06_rprec/w_star
      q_star=T_star
      else
      w_star=sign((g/theta_mean*abs(wt_s)*z_i)**(1._rprec/3._rprec),wt_s)
      T_star=wt_s/w_star
      q_star=T_star
      end if
        
      $if ($MPI)
         print *,'Modified Log Profile for IC for coord = ',coord
      $else
         print *,'Modified Log Profile for IC'
      $endif

      do jz=1,nz
         $if ($MPI)
            z = (coord*(nz-1) + jz - 0.5_rprec) * dz
         $else
            z=(real(jz)-0.5_rprec)*dz
         $endif
!c IC in equilibrium with rough surface (rough dominates in effective zo)
         arg2=z/(sum(zo)/float(nx*ny))
         arg=(1._rprec/vonk)*log(arg2)!-1./(2.*vonk*z_i*z_i)*z*z
         if (coriolis_forcing) then
         ubar(jz)=ug_dim/u_star
         vbar(jz)=vg_dim/u_star
         wbar(jz)=0._rprec
! Note that ug and vg have already been non-dimensionalized in param.f90
         else
         ubar(jz)=arg
         vbar(jz)=0._rprec
         wbar(jz)=0._rprec
         end if
!C sc: I changed around some parenthesis here
        if (z.gt.(perturb_height_factor*z_i)) then
        ubar(jz)=ubar(jz-1)
        end if
       end do

  rms = 3._rprec
  do jz=1,nz
  $if ($MPI)
    jz_abs = coord * (nz-1) + jz
    z = (coord * (nz-1) + jz - 0.5_rprec) * dz * z_i
  $else
    jz_abs = jz
    z = (jz-.5_rprec) * dz * z_i
  $endif
  seed = -80 - jz_abs  !--trying to make consistent init for MPI
    do jy=1,ny
      do jx=1,nx
!c...Ran3 returns uniform RV between 0 and 1. (sigma_rv=0.289)
!c...Taking std dev of vel as 1 at all heights

!cVK Note that if you put wt_s = 0 symbolizing neutal conditions
!c u should also put L_z=z_i i.e. the inversion layer height should
!c be equal to the height of the domain and in that case the second
!c part of the subsequent if loop will never execute. This is
!c ensured by putting an OR statement in the if clause, which makes 
!c sure that only the first part of if block is executed and not the
!c block after else

       if (z.le.perturb_height_factor*z_i) then
         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
         u(jx,jy,jz)=noise*(1._rprec-z/z_i)*w_star/u_star+ubar(jz) !noise
         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
         v(jx,jy,jz)=noise*(1._rprec-z/z_i)*w_star/u_star+vbar(jz) !noise
         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
         w(jx,jy,jz)=noise*(1._rprec-z/z_i)*w_star/u_star+wbar(jz)
         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
         theta(jx,jy,jz)=(theta_mean+10._rprec*noise*(1-z/z_i)*T_star)/T_scale
         q(jx,jy,jz)=q_mix
       else
         u(jx,jy,jz)=ubar(jz)
         v(jx,jy,jz)=vbar(jz)
         w(jx,jy,jz)=wbar(jz)
         theta(jx,jy,jz)=(theta_mean+(z-perturb_height_factor*z_i)*inv_strength)/T_scale
         q(jx,jy,jz)=q_mix
        end if
       end do
     end do
  end do

  !...BC for W
  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    w(1:nx, 1:ny, 1) = 0._rprec
  end if
  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
    w(1:nx, 1:ny, nz) = 0._rprec
  endif

  !...BC for U, V & T
  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
    u(1:nx, 1:ny, nz) = u(1:nx, 1:ny, nz-1)
    v(1:nx, 1:ny, nz) = v(1:nx, 1:ny, nz-1)
    theta(1:nx, 1:ny, nz) = theta(1:nx, 1:ny, nz-1)+dTdz_top/T_scale*z_i*dz
  end if

!VK Display the mean vertical profiles of the initialized variables on the screen
open(unit=44,file=path//'output/init_profiles.dat',status="unknown",position="append")
do jz=1,nz
     $if ($MPI)
       z = (coord*(nz-1) + jz - 0.5_rprec) * dz
     $else
       z = (jz - 0.5_rprec) * dz
     $endif
     write(6,7781) jz,z,(sum(u(:,:,jz))/float(nx*ny))*u_star,(sum(v(:,:,jz))/&
     float(nx*ny))*u_star,(sum(w(:,:,jz))/float(nx*ny))*u_star,&
     (sum(theta(:,:,jz))/float(nx*ny))*T_scale
     write(44,7781) jz,z,(sum(u(:,:,jz))/float(nx*ny))*u_star,(sum(v(:,:,jz))/&
     float(nx*ny))*u_star,(sum(w(:,:,jz))/float(nx*ny))*u_star,&
     (sum(theta(:,:,jz))/float(nx*ny))*T_scale
end do
close(44)
7781 format('jz, z, ubar, vbar, wbar,T_bar:',(1x,I3,1x,F9.4,1x,F9.4,1x,F9.4,1x,F9.4,1x,F9.4))


end subroutine ic_scal_GABLS

!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine budget_TKE_scalar
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! To compute the budget of TKE, we need the following terms:
! <u'w'>*d/dz(<U>), <v'w'>*d/dz(<V>), (g/<T>)*<w'T'>, [0.5*d/dz(<w'e>) where
! e=<u'^2>+<v'^2>+<w'^2>], d/dz(<w'p'>), dissipation (Nu_t*S^2),d/dz(<u'\tau{13}>)+
! d/dz(v'\tau{23}), <\tau{13}>d/dz(<U>)+<\tau{23}>d/dz(<V>)
! Of the eight terms above, we can directly compute terms 1,2,3,8 from the variables
! calculated/outputted in avgslice.f90 and scalar_slice.f90
! So, the rest 4 terms will be computed/outputted here
! Similarly for temperature variance budget, we need
! <w'T'>*d/dz(<T>), d/dz(<w*T^2>) and we already have term 1 from 
! scalar_slice.f90.. so we just need to compute term 2 here 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use sim_param,only:path,u,v,w,theta,p,txz,tyz
use param,only:dz,p_count,c_count,jt
use sgsmodule,only:dissip
implicit none
integer::i,j,k
real(kind=rprec),dimension(nx,nz-1),save::awu2,awv2,awT2,auv,awp
real(kind=rprec),dimension(nx,nz-1),save::autau13,avtau23,adissip
real(kind=rprec)::twu2,twv2,twt2,tuv,twp,arg1,arg2,arg3,arg4,arg5,arg6,fr
real(kind=rprec)::arg7
real(kind=rprec)::tutau13,tvtau23,tdissip
real(kind=rprec),dimension(:,:),allocatable::avg_budget_out
real(kind=rprec),dimension(nz-1)::ubar_profile,vbar_profile,Tbar_profile

fr=(1._rprec/real(p_count,kind=rprec))*real(c_count,kind=rprec)

do k=1,nz-1
ubar_profile(k)=sum(u(1:nx,1:ny,k))/(nx*ny)
vbar_profile(k)=sum(v(1:nx,1:ny,k))/(nx*ny)
Tbar_profile(k)=sum(theta(1:nx,1:ny,k))/(nx*ny)
end do

do k=1,Nz-1
do i=1,Nx
   twu2=0._rprec;twv2=0._rprec;twT2=0._rprec;tuv=0._rprec;twp=0._rprec;
   tutau13=0._rprec;tvtau23=0._rprec;tdissip=0._rprec

   do j=1,Ny
      if((k .eq. 1) .AND. ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0))) then  
         arg1=0._rprec;arg2=0._rprec;arg3=0._rprec;arg4=0._rprec;arg5=0._rprec;arg6=0._rprec
      else  
         arg4=(p(i,j,k)+p(i,j,k-1))/2.
         arg5=(u(i,j,k)-ubar_profile(k)+u(i,j,k-1)-ubar_profile(k-1))/2.
         arg6=(v(i,j,k)-vbar_profile(k)+v(i,j,k-1)-vbar_profile(k-1))/2.
         arg7=(theta(i,j,k)-Tbar_profile(k)+theta(i,j,k-1)-Tbar_profile(k-1))/2.
      end if
      twu2=twu2+w(i,j,k)*arg5*arg5 !directly computes <w(u')^2>
      twv2=twv2+w(i,j,k)*arg6*arg6 !directly computes <w(v')^2>
! Also note that since <w> = 0, there is no need to calculate it as
! <w^3> = <w'^3> and we are outputting <w^3> in avgslice
      twT2=twT2+w(i,j,k)*arg7
      twp=twp+w(i,j,k)*arg4
! <u'v'> is not as simple as <u'w'> since <u> .ne. 0 whereas <w>=0
! therefore, we directly calculate <u'v'> here
      tuv=tuv+(u(i,j,k)-ubar_profile(k))*(v(i,j,k)-vbar_profile(k))
      tutau13=tutau13+arg5*txz(i,j,k) !computes SGS transport of TKE i.e. <u'\tau_{13}>
      tvtau23=tvtau23+arg6*tyz(i,j,k) !computes SGS transport of TKE i.e. <v'\tau_{13}>
      tdissip=tdissip+dissip(i,j,k) ! outputs dissip calculated in sgs stag..
   end do
   awu2(i,k)=awu2(i,k)+(fr)*twu2/Ny
   awv2(i,k)=awv2(i,k)+(fr)*twv2/Ny
   awT2(i,k)=awT2(i,k)+(fr)*twT2/Ny
   awp(i,k)=awp(i,k)+(fr)*twp/Ny
   auv(i,k)=auv(i,k)+(fr)*tuv/Ny
   autau13(i,k)=autau13(i,k)+(fr)*tutau13/Ny
   avtau23(i,k)=avtau23(i,k)+(fr)*tvtau23/Ny
   adissip(i,k)=adissip(i,k)+(fr)*tdissip/Ny
end do
end do

if (mod(jt,p_count)==0) then
  allocate(avg_budget_out(1:nx,1:(nz_tot-1)));
 ! call collocate_MPI_averages_N(awu2,avg_budget_out,61,'awu2')
 ! call collocate_MPI_averages_N(awv2,avg_budget_out,62,'awv2')
  call collocate_MPI_averages_N(awT2,avg_budget_out,63,'awT2')
  call collocate_MPI_averages_N(awp,avg_budget_out,64,'awp_sc')
 ! call collocate_MPI_averages_N(auv,avg_budget_out,65,'auv')
 ! call collocate_MPI_averages_N(autau13,avg_budget_out,66,'autau13')
 ! call collocate_MPI_averages_N(avtau23,avg_budget_out,67,'avtau23')
  call collocate_MPI_averages_N(adissip,avg_budget_out,68,'adissip');deallocate(avg_budget_out)
!VK Zero out the outputted averages !!
  awu2=0._rprec;awv2=0._rprec;awT2=0._rprec;awp=0._rprec;auv=0._rprec
  autau13=0._rprec;avtau23=0._rprec;adissip=0._rprec
end if

end subroutine budget_TKE_scalar

!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
subroutine ic_scal_GABLS_diurnal()
!subroutine ic_scal(u,v,w,theta,q,sign_wt_s,theta_mean)
!c...Log profile that is modified to flatten at z=z_i
!c .. Modified by Vijayant to put scalars back
! Last modified April 14, 2004

implicit none
real(kind=rprec),dimension(nz)::ubar,vbar,wbar,T_bar,tke_bar,tke_sum
real(kind=rprec)::rms, noise, arg, arg2,theta_mean
real(kind=rprec)::z,w_star,T_star,q_star,ran3,perturb_height_factor
real(kind=rprec),dimension(:,:,:),allocatable::wtemp
integer::jx,jy,jz,seed,jz_abs,ii
character (64) :: fname, temp

tke_bar=0._rprec !initialize mean TKE profile as 0
tke_sum=0._rprec !initialize mean TKE profile as 0

      theta_mean=T_init
      perturb_height_factor=800._rprec/z_i
      
      if (wt_s .eq. 0.0_rprec) then
      w_star=(g/theta_mean*0.06_rprec*z_i)**(1._rprec/3._rprec)
! w_star is of O(1) with z_i=500 and wt_s=0.06
      T_star=0.06_rprec/w_star
      q_star=T_star
      else
      w_star=sign((g/theta_mean*abs(wt_s)*z_i)**(1._rprec/3._rprec),wt_s)
      T_star=wt_s/w_star
      q_star=T_star
      end if

         $if ($MPI)
            print *,'Modified Log Profile for IC for coord = ',coord
         $else
            print *,'Modified Log Profile for IC'
         $endif
       do jz=1,nz
         $if ($MPI)
            z = (coord*(nz-1) + jz - 0.5_rprec) * dz
         $else
            z=(real(jz)-0.5_rprec)*dz
         $endif
!c IC in equilibrium with rough surface (rough dominates in effective zo)
        arg2=z/(sum(zo)/float(nx*ny))
        arg=(1._rprec/vonk)*log(arg2)!-1./(2.*vonk*z_i*z_i)*z*z
        if (coriolis_forcing) then
        ubar(jz)=ug_dim/u_star
        vbar(jz)=vg_dim/u_star
        wbar(jz)=0._rprec
! Note that ug and vg have already been non-dimensionalized in param.f90
        else
        ubar(jz)=arg
        vbar(jz)=0._rprec
        wbar(jz)=0._rprec
        end if
!C sc: I changed around some parenthesis here
        if (z.gt.(perturb_height_factor*z_i)) then
          ubar(jz)=ubar(jz-1)
        end if
       end do
       
5169     format(1400(E16.10))

  do jz=1,nz
  $if ($MPI)
    jz_abs = coord * (nz-1) + jz
    z = (coord * (nz-1) + jz - 0.5_rprec) * dz * z_i
  $else
    jz_abs = jz
    z = (jz-.5_rprec) * dz * z_i
  $endif
  seed = -80 - jz_abs  !--trying to make consistent init for MPI
  if (z .lt. 800._rprec) tke_bar(jz)=0.5_rprec*(1-z/800._rprec)

  rms=(2._rprec*tke_bar(jz)/3._rprec)**0.5_rprec !2/3(1/2(u^2))=1/3(u^2) for each (u,v,w)
    do jy=1,ny
      do jx=1,nx
!c...Ran3 returns uniform RV between 0 and 1. (sigma_rv=0.289)
!c...Taking std dev of vel as 1 at all heights

!cVK Note that if you put wt_s = 0 symbolizing neutal conditions
!c u should also put L_z=z_i i.e. the inversion layer height should
!c be equal to the height of the domain and in that case the second
!c part of the subsequent if loop will never execute. This is
!c ensured by putting an OR statement in the if clause, which makes 
!c sure that only the first part of if block is executed and not the
!c block after else

       if (z.le.perturb_height_factor*z_i) then
         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
         u(jx,jy,jz)=noise/u_star+ubar(jz) !noise
         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
         v(jx,jy,jz)=noise/u_star+vbar(jz) !noise
         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
         w(jx,jy,jz)=noise/u_star+wbar(jz) !noise
         
         call T_pos_gabls(theta_mean,z)
         theta(jx,jy,jz)=theta_mean/T_scale
         q(jx,jy,jz)=q_mix
       else
         u(jx,jy,jz)=ubar(jz)
         v(jx,jy,jz)=vbar(jz)
         w(jx,jy,jz)=wbar(jz)
         call T_pos_gabls(theta_mean,z)
         theta(jx,jy,jz)=theta_mean/T_scale
         q(jx,jy,jz)=q_mix
        end if
       end do
     end do
  end do


  !...BC for W
  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    w(1:nx, 1:ny, 1) = 0._rprec
  end if
  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
    w(1:nx, 1:ny, nz) = 0._rprec
  endif

  !...BC for U, V & T
  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
    u(1:nx, 1:ny, nz) = u(1:nx, 1:ny, nz-1)
    v(1:nx, 1:ny, nz) = v(1:nx, 1:ny, nz-1)
    theta(1:nx, 1:ny, nz) = theta(1:nx, 1:ny, nz-1)+dTdz_top/T_scale*z_i*dz
  end if

!calculate initial TKE profile
allocate(wtemp(ld,ny,nz)) !wtemp allocated
wtemp=0._rprec
wtemp(:,:,1:nz-1)=0.5_rprec*(w(:,:,1:nz-1)+w(:,:,2:nz));
  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
   wtemp(:,:,nz)=w(:,:,nz);
  end if
do jz=1,nz
   tke_sum(jz)=tke_sum(jz)+0.5_rprec*sum((u(1:nx,1:ny,jz)-sum(u(1:nx,1:ny,jz))/(nx*ny))**2)/(nx*ny)
   tke_sum(jz)=tke_sum(jz)+0.5_rprec*sum((v(1:nx,1:ny,jz)-sum(v(1:nx,1:ny,jz))/(nx*ny))**2)/(nx*ny)
   tke_sum(jz)=tke_sum(jz)+0.5_rprec*sum((wtemp(1:nx,1:ny,jz)-sum(wtemp(1:nx,1:ny,jz))/(nx*ny))**2)/(nx*ny)
end do
deallocate(wtemp)

!VK Display the mean vertical profiles of the initialized variables on the
!screen
    write(fname,'(A,i6.6,A)')path//'output/init_profiles.dat'
    $if ($MPI)
      write (temp, '(".c",i0)') coord
      fname = trim (fname) // temp
    $endif
open(unit=44,file=fname,status="unknown",position="append")
do jz=1,nz
     $if ($MPI)
       z = (coord*(nz-1) + jz - 0.5_rprec) * dz
     $else
       z = (jz - 0.5_rprec) * dz
     $endif
     write(6,7781) jz,z,(sum(u(:,:,jz))/float(nx*ny))*u_star,(sum(v(:,:,jz))/&
     float(nx*ny))*u_star,(sum(w(:,:,jz))/float(nx*ny))*u_star,&
     (sum(theta(:,:,jz))/float(nx*ny))*T_scale,tke_sum(jz)*u_star*u_star
     write(44,7781) jz,z,(sum(u(:,:,jz))/float(nx*ny))*u_star,(sum(v(:,:,jz))/&
     float(nx*ny))*u_star,(sum(w(:,:,jz))/float(nx*ny))*u_star,&
     (sum(theta(:,:,jz))/float(nx*ny))*T_scale,real(tke_sum(jz))*u_star*u_star
end do
close(44)

    write(fname,'(A,i6.6,A)')path//'output/init_profiles_3d.bin'
    $if ($MPI)
      write (temp, '(".c",i0)') coord
      fname = trim (fname) // temp
    $endif
open(unit=44,file=fname,form="unformatted")
write(44) real(u),real(v),real(w),real(theta);close(44)

7781 format('jz, z, ubar, vbar, wbar,T_bar,TKE:',(1x,I3,1x,F9.4,1x,F9.4,1x,F9.4,1x,F9.4,1x,F9.4,1x,F16.10))

end subroutine ic_scal_GABLS_diurnal

!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
subroutine T_pos_gabls(T1,z_loc)
implicit none
real(kind=rprec):: T1,z_loc
if (z_loc<=200._rprec) then
    T1=288._rprec-z_loc*(288._rprec-286._rprec)/200._rprec;
else if (z_loc>200._rprec .AND. z_loc<=850._rprec) then
    T1=286._rprec
else if (z_loc>850._rprec .AND. z_loc<=1000._rprec) then
    T1=286._rprec+(z_loc-850._rprec)*(292._rprec-286._rprec)/(1000._rprec-850._rprec);
else if (z_loc>1000._rprec .AND. z_loc<=2000._rprec) then
    T1=292._rprec+(z_loc-1000._rprec)*(300._rprec-292._rprec)/(2000._rprec-1000._rprec);
else if (z_loc>2000._rprec .AND. z_loc<=3500._rprec) then
    T1=300._rprec+(z_loc-2000._rprec)*(310._rprec-300._rprec)/(3500._rprec-2000._rprec);
else if (z_loc>3500._rprec .AND. z_loc<=4000._rprec) then
    T1=310._rprec+(z_loc-3500._rprec)*(312._rprec-310._rprec)/(4000._rprec-3500._rprec);
else
    print *,'z not contained in the if block range !!!'
    stop
end if
end subroutine T_pos_gabls
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX

end module scalars_module2
