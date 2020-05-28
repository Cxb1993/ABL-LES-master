module scalars_module
! HUMIDITY subroutines in place but not yet turned on !!
use types,only:rprec
use param 
use sim_param,only:u,v,w,theta,q,path,u_old
use bottombc ! Includes patches subroutine
use sgsmodule,only:Nu_t,Pr_t
implicit none
integer, parameter:: tag_counter = 200
logical, parameter:: SCALAR_DEBUG=.FALSE.
!!!!!!--------------------------------------------
! Part I of the scalar files - contains the basic subroutines
! Also look at scalars_module2.f90 for other subroutines !! 
! CONTAINS subroutines :
! theta_all_in_one - Performs the various steps for theta calculation
! humidity_all_in_one - Performs the various steps for humidity
! scalar_RHS_calc - computes the RHS side of the scalar evolution equation
! calcbeta - computes the buoyancy term for temperature
! step_scalar - time evolves the scalar
! obukhov - computes the obukhov similarity terms for use in scalars,wallstress and derivwall
! Authored by Vijayant Kumar
! Last modified - April 24, 2004
!!!!!!--------------------------------------------
 
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif
!! KMG: added sgs_t1 and sgs_t2 for horizontal sgs heat fluxes Pi_1 and Pi_2
real(kind=rprec),dimension(ld,ny,$lbz:nz):: beta_scal,Pr_
real(kind=rprec),dimension(ld,ny,$lbz:nz):: dTdz,dqdz,dTdx,dTdy
! Only ones needed for output..Might need to add x and y derivatives here in case they need to be outputted
! Right now they are in the "scalar"_all_in_one routines below !!
real(kind=rprec),dimension(ld,ny,$lbz:nz):: RHS_Tf,RHS_T,RHS_qf,RHS_q
real(kind=rprec), dimension(ld,ny,$lbz:nz):: sgs_t3,sgs_q3,sgs_t1,sgs_t2 ! defines the surface sgs flux

real(kind=rprec),dimension(nx,ny)::L,wstar ! defines obukhov length and convective vel scale, w_star
real(kind=rprec),dimension(ld,ny)::T_s_filtered ! filtered T_s for calc of wT_s

! Now define local u_star in bottombc.f90 and calculate in wallstress.f90 and use that value everywhere else
integer, parameter:: obukhov_output=0 !Controls whether obukhov variables are outputted by scalar_slice
integer, parameter:: wt_s_vector_dim1=no_days*86400/300+1

real(kind=rprec),dimension(wt_s_vector_dim1,1) :: wt_s_vector
! Variables for heterogeneity analysis
! hetero_array_freqz = number of time steps equivalent to 20 seconds
integer,parameter:: hetero_array_freqz=int(20/dt_dim),hetero_count_out=p_count
integer,save::time_ind

!--MM 2015: Added Artificially Moving(Heating) Inversion:
integer, parameter :: nz_global = (nz-1) * nproc
integer ::k_global
real ::dt_avg
real(kind=rprec),dimension(ld,ny,$lbz:nz)::theta_old

contains
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
subroutine theta_all_in_one
use topbc,only:sponge
use test_filtermodule
implicit none
real::wt_s_current,dummy_t
integer::jz,ios,counter
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif
real(kind=rprec),dimension(ld,ny,$lbz:nz)::dTdx,dTdy
real(kind=rprec),dimension($lbz:nz)::sponge_theta

 if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
     if (jan_diurnal_run) then !perform time variable heat flux
         if (jt .eq. SCAL_init) then
            counter=1
            open(unit=56,file=path//'heatflux_diurnal_run.txt')
            do
               read(56,*,iostat=ios) wt_s_vector(counter,1)
               if (ios/=0) exit
               counter=counter+1
            end do
            print *,'length of wt_s_vector = ',size(wt_s_vector)
         end if
! note that jt_total is always jt-1 as it is updated only at the end of each step
! and hence the formula below reflects that !!
         if (GABLS_diurnal_test) then !read a continuous WT tseries for GABLS diurnal
         wt_s_current=wt_s_vector(jt_total+1,1)
         else ! for the actual jan diurnal run 
         wt_s_current=wt_s_vector(floor(((jt_total)*dt*z_i/u_star)/300._rprec)+1,1)
         end if

         print '(A,F12.6,A,I3.3)','wt_s = ',wt_s_current&
         ,' at time (hour) = ',floor(((jt_total+1)*dt*z_i/u_star)/3600)
     else    
        wt_s_current=wt_s
        ! SKS
        ! wt_s = 0.20 ( w*T* or w* theta* )
        ! SKS
     end if
 end if

! SKS
  if (model == 1 .or. model == 2 .or. model == 3 .or. model == 4 .or. model == 5) then
     Pr_ = Pr
  else
     do jz = $lbz,nz ! ubc_jz
        Pr_(:,:,jz) = Pr_t(:,:,jz) ! Pr
     end do
  end if
! SKS

! Right now set the Prandtl num matrix equal to a constant Prandtl
! number as specified in param. could use the dynamic model ideal to compute Pr as well !!
! The plan-averaged dynamic Prandtl number model is already coded. just need to put it in !!

call filt_da(theta,dTdx,dTdy)
call ddz_uv (dTdz,theta)

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
   dTdz(:,:,Nz)=dTdz_top/T_scale*z_i ! Valid for temperature
end if

$if ($MPI)
! print *,'synchronizing in theta_all_in for coord = ',coord
! Need to synchronize w and dTdz across the processors for uvp node
! based computation in scalar_RHS_calc (calculation of RHS_m)
  call mpi_sendrecv (w(1, 1, 1), ld*ny, MPI_RPREC, down, tag_counter+1,  &
                     w(1, 1, nz), ld*ny, MPI_RPREC, up, tag_counter+1,   &
                     comm, status, ierr)   
  call mpi_sendrecv (dTdz(1, 1, 1), ld*ny, MPI_RPREC, down, tag_counter+3,  &
                     dTdz(1, 1, nz), ld*ny, MPI_RPREC, up, tag_counter+3,   &
                     comm, status, ierr)   

! Also need to synchronize Nu_t across the processors for 
! computation in scalar_RHS_calc (calculation of RHS_m)
  call mpi_sendrecv (Nu_t(1, 1, 1), ld*ny, MPI_RPREC, down, tag_counter+3,  &
                     Nu_t(1, 1, nz), ld*ny, MPI_RPREC, up, tag_counter+3,   &
                     comm, status, ierr)   
  call mpi_sendrecv (Nu_t(1, 1, nz-1), ld*ny, MPI_RPREC, up, tag_counter+4,  &
                     Nu_t(1, 1, 0), ld*ny, MPI_RPREC, down, tag_counter+4,   &
                     comm, status, ierr)
$endif

if (S_FLAG) then
   RHS_Tf=RHS_T

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then !Do only for the first processor
if (GABLS_test) then
    if (jt .ge. max(SCAL_init,DYN_init)) then
       T_s=265._rprec-0.25_rprec*(jt_total*dt*z_i/u_star/3600) !continuous change in T_s
    else
       T_s=265._rprec
    end if
! SKS
if(mod(jt,1) == 100) then        ! 100 because wbase = 100
! SKS
print*,'(A,F12.6,A,I3.3)','T_s = ',sum(T_s)/real(nx*ny)&
,' at time (hour) = ',floor((jt_total*dt*z_i/u_star)/3600)
! SKS
end if
! SKS
T_s=T_s/T_scale
elseif (GABLS_diurnal_test) then ! T_s conditions for the new GABLS Diurnal run
  dummy_t=(jt_total*dt*z_i/u_star)/3600
  print *,'dummy_t,jt_total',dummy_t,jt_total
  if (dummy_t .le. 17.4_rprec) then
     T_s=-10._rprec-25._rprec*cos(dummy_t*0.22_rprec+0.2_rprec) !DAY 1 
  elseif ((dummy_t .gt. 17.4_rprec) .and. (dummy_t .le. 30._rprec)) then 
     T_s=-0.54_rprec*dummy_t+15.2_rprec !NIGHT 1
  elseif ((dummy_t .gt. 30._rprec) .and. (dummy_t .le. 41.9_rprec)) then
     T_s=-7._rprec-25._rprec*cos(dummy_t*0.21_rprec+1.8_rprec)!DAY 2 
  elseif ((dummy_t .gt. 41.9_rprec) .and. (dummy_t .le. 53.3_rprec)) then
     T_s=-0.37_rprec*dummy_t+18._rprec !NIGHT 2
  elseif ((dummy_t .gt. 53.3_rprec) .and. (dummy_t .le. 65.6_rprec)) then 
     T_s=-4._rprec-25._rprec*cos(dummy_t*0.22_rprec+2.5_rprec)!DAY 3
  elseif (dummy_t .gt. 65.6_rprec) then
     T_s=4.4_rprec !NIGHT 3
  end if
  print '(A,F12.6,A,I3.3)','T_s = ',sum(T_s)/real(nx*ny)&
  ,' at time (hour) = ',floor((jt_total*dt*z_i/u_star)/3600)
  T_s=(T_s+273.15_rprec)/T_scale
end if
end if !end loop for surface boundary conditions only to be performed by process 0

! Perform test filtering of T_s for calculation of surf fluxes
     if ((jt_total .eq. SCAL_init) .and. (lbc .eq. 0) .and. (remote_homog_flag .eq. 0)) then
       print *,'T_s b4 filtering',sqrt(sum((T_s-sum(T_s)/float(nx*ny))**2))/float(nx*ny)
       T_s_filtered(1:nx,1:ny)=T_s
       call test_filter(T_s_filtered,G_test)
       T_s=T_s_filtered(1:nx,1:ny)
       print *,'T_s after filtering',sqrt(sum((T_s-sum(T_s)/float(nx*ny))**2))/float(nx*ny)
     end if
     if ((jt_total .eq. nsteps-1) .and. (lbc .eq. 0) .and. (remote_homog_flag .eq. 0)) then
        print *,'T_s after filtering',sqrt(sum((T_s-sum(T_s)/float(nx*ny))**2))/float(nx*ny)
     end if

! call scalar_RHS_calc(T_s,z_os,RHS_T,sgs_t3,jt,psi_h,phi_h,Pr_,wt_s_current)
theta_old(:,:,:)=theta(:,:,:)

call scalar_RHS_calc(theta,dTdx,dTdy,dTdz,T_s,z_os,RHS_T,sgs_t1,sgs_t2,sgs_t3,wt_s_current)

! subroutine looks like this later
! subroutine scalar_RHS_calc(scalar,dsdx,dsdy,dsdz,S_Surf,z_os,RHS,sgs_horx,sgs_hory,sgs_vert,surf_flux_current)


if (ubc==1 .and. damping_method==2) then ! add the damping term to the scalar equation
   do jz=1,nz-1
   RHS_T(1:nx,1:ny,jz)=RHS_T(1:nx,1:ny,jz)-0.5_rprec*(sponge(jz)+sponge(jz+1))*&
                 (theta(1:nx,1:ny,jz)-sum(theta(1:nx,1:ny,jz))/(nx*ny))
   end do
end if

call calcbeta(theta) ! Calculates the buoyancy term which gets added to the vertical momentum equation
! call calcbeta(theta,beta_scal)

if (jt .eq. SCAL_INIT .and. (.not. initsc)) then
   RHS_Tf=RHS_T
end if

call step_scalar(theta,RHS_T,RHS_Tf)
end if

!--- MM: Artificial heating of the inversion layer:

!$if ($MPI)
!if (coord==2) then
!        dt_avg=sum(theta(1:nx,1:ny,1:nz-1))/(nx*ny*(nz-1))
!        dt_avg=dt_avg-sum(theta_old(1:nx,1:ny,1:nz-1))/(nx*ny*(nz-1))
! do jz=nproc/2-1,nproc-1
!        call mpi_send (dt_avg, 1, MPI_REAL,jz , 8101, comm, ierr)
! end do
!end if


!if (coord>nproc/2-2) then
!  !--wait for u(n), gam(n) from "up"
!  call mpi_recv (dt_avg, 1, MPI_REAL, 2, 8101, comm, status, ierr)
!end if
!$endif

!if (jt<nsteps/60) then
!        u_old=0
!end if

!do jz=1,nz-1
!    k_global = jz + coord *(nz-1)
!    if (k_global > nz_global*0.67) then              
!               theta(:,:,jz)=theta(:,:,jz)+dt_avg+sum(theta_old(1:nx,1:ny,jz)-theta(1:nx,1:ny,jz))/(nx*ny)
!    end if
!end do

!u_old=-(sum(u(1:nx,1:ny,1:nz-1))/(nx*ny*(nz-1)))*((-0.005)*coriol*T_scale)/g
!MM -----


$if ($MPI)
  call mpi_sendrecv (theta(1, 1, 1), ld*ny, MPI_RPREC, down, tag_counter+7,  &
                     theta(1, 1, nz), ld*ny, MPI_RPREC, up, tag_counter+7,   &
                     comm, status, ierr)   
  call mpi_sendrecv (theta(1, 1, nz-1), ld*ny, MPI_RPREC, up, tag_counter+8,  &
                     theta(1, 1, 0), ld*ny, MPI_RPREC, down, tag_counter+8,   &
                     comm, status, ierr)
$endif

end subroutine theta_all_in_one

!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
subroutine calcbeta (scalar)
! subroutine calcbeta (scalar, beta_scal)
! This calculates the buoyancy term (beta_scal) to be added to the vertical momentum equation for temperature
! Authored by Vijayant Kumar
! Last updated April 14, 2004
implicit none
integer::i, j, k, jz_min
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif
real(kind=rprec),dimension(ld,ny,$lbz:nz),intent(in)::scalar
real(kind=rprec),dimension(nz)::scalar_bar
real(kind=rprec)::g_hat,above, below
!..Non-dimensionalize gravity
g_hat=g*(z_i/(u_star**2))
beta_scal=0._rprec

! Note Beta is stored on W nodes, but Theta is on UVP nodes
! We do not time-advance the ground nodes, so start at k=2
! VK: Inserted the averaging code inside this file itself rather than doing it in prof
do k=$lbz,nz
scalar_bar(k)=0.0    
   do j=1,ny
      do i=1,nx
        scalar_bar(k)=scalar_bar(k)+scalar(i,j,k)
      end do
   end do
scalar_bar(k)=scalar_bar(k)/(nx*ny)
end do
! We do not time-advance the ground nodes, so start at k=2
! For the MPI case, this means that we start from jz=2 for
! coord=0 and jz=1 otherwise... enable by an if statment

 if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    jz_min = 2
 else
    jz_min = 1
 end if

do k=jz_min,Nz-1
      do j=1,Ny
             do i=1,nx
                above=(scalar(i,j,k)-scalar_bar(k))/scalar_bar(k)
                below=(scalar(i,j,k-1)-scalar_bar(k-1))/scalar_bar(k-1)
                beta_scal(i,j,k)=g_hat*(above + below)/2._rprec
             end do
      end do
end do
return
end subroutine calcbeta
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------

subroutine scalar_RHS_calc(scalar,dsdx,dsdy,dsdz,S_Surf,z_os,RHS,sgs_horx,sgs_hory,sgs_vert,surf_flux_current)
use test_filtermodule
implicit none
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif
integer::i,j,k,jz
integer::jz_min,ubc_jz
real:: surf_flux_current,crap2
real(kind=rprec),dimension(ld,ny,$lbz:nz):: dsdx,dsdy,dsdz
real(kind=rprec),dimension(ld,ny,$lbz:nz):: RHS,temp
real(kind=rprec),dimension(ld,ny,$lbz:nz):: scalar
real(kind=rprec),dimension(ld,ny,$lbz:nz):: dtemp,sgs_vert,sgs_horx,sgs_hory
real(kind=rprec),dimension(ld_big,ny2,$lbz:nz):: u_m,v_m,w_m,dsdx_m,dsdy_m,dsdz_m
 real(kind=rprec),dimension(ld_big,ny2,$lbz:nz):: RHS_m
real(kind=rprec),dimension(nx,ny):: ustar_local,S_Surf,surf_flux,z_os
real(kind=rprec),dimension(ld,ny):: scalar_node_1 ! test filtered and used for computing surface flux
real(kind=rprec),dimension (ptypes):: ustar,wt_s2
character (64) :: fname_hetero

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then

   ustar=0.
if (patch_flag==1) then !if block 101
   if (spatial_flux_flag) then !SPATIAL_FLUX IF BLOCK! if block 102
       if (jt .eq. SCAL_init) then
        print *,'reading spatial flux data from file !!'
        open(unit=77,file=path//'spatial_flux.dat',status='unknown')
          do j=1,ny
              read(77,5168) (surf_flux(i,j),i=1,nx)
          end do
          surf_flux=surf_flux/u_star/T_scale !The readin surf_flux is dimensional - so non-dimensionalize
       else if(jt .GT. SCAL_init) then
          surf_flux=sgs_t3(1:nx,1:ny,1)
       end if !end for if loop for jt .eq. SCAL_init
       ustar_local=ustar_avg !set ustar as value computed in obukhov

        do j=1,ny
        do i=1,nx
          dsdz(i,j,1) =-phi_h(i,j)*surf_flux(i,j)/(ustar_local(i,j)*vonk*DZ/2._rprec) !set the gradient at the first point
        end do
        end do

!!VIJ!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! Need to work on this later to remove the above do loop
! dsdz(:,:,1) =-phi_h(:,:)*surf_flux(:,:)/(ustar_local(:,:)*vonk*DZ/2._rprec) !set the gradient at the first point
!!!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

5169     format(1400(E16.10))
   else    ! SPATIAL FLUX IF BLOCK CONTINUES !block 102 conts.
    do k=1,ptypes
      wt_s2(k)=(-1.)**(k+1)*surf_flux_current 
!VK Creates patches of type ptypes with alternating signs of heat flux
    end do
!c This computes the average value of the scalar
!c S_surf refers to the value assigned to the surface according to
!c routine patches.f based on number of patches, and parameters in
!c dimen.h
! Also calculated is a local averaged value of u_star from the subgrid
! stress at the wall
    do j=1,ny
    do i=1,nx
    do k=1,ptypes
    if (patch(i,j)==k) then
     ustar(patch(i,j))=ustar(patch(i,j))+ustar_avg(i,j)/patchnum(patch(i,j))
    end if
    end do
    end do
    end do
!c  distribute ustar value to patches
    do j=1,ny
    do i=1,nx
    do k=1,ptypes
        if (patch(i,j)==k) then
          ustar_local(i,j)=ustar(k)
        end if
    end do
    end do
    end do

!c..Compute surface flux and dsdz at z=DZ/2
    do j=1,Ny
    do i=1,Nx
! lbc=1 is used for prescribing the surface flux while
! lbc=0 has been used to prescribe the temperature
 if (lbc==1.and.scalar(1,1,1)<2) then
   do k=1,ptypes
     if (patch(i,j)==k) then
     surf_flux(i,j)=wt_s2(k)/T_scale/u_star
! The u_star above is coming from dimen.h = Ug for coriolis and is not
! the local u_star computed from stress at the surface.
     end if
  end do
  
 else if (lbc==0.and.scalar(1,1,1)<2) then
     ustar_local=ustar_avg
     surf_flux(i,j)=(S_Surf(i,j)-scalar(i,j,1))*vonk*ustar_local(i,j)&
    /(dlog(dz/(2._rprec*z_os(i,j)))-psi_h(i,j))
!c equ. 4.28. in brutsaert
! Impose the gradient value at the first node.
! Note that ustar2 is a variable coming from sgs stresses at the wall,
! ustar2 has a single value for all the nodes in a single patch.
! phi_h is obtained from the routine obukhov.f, surf_flux is as computed above
! and vonk and dz are constants. Everything is in dimensionless form
 end if

!c....Now we have the lowest dsdz on the UVP nodes all others on w nodes
dsdz(i,j,1) =-phi_h(i,j)*surf_flux(i,j)/(ustar_local(i,j)*vonk*DZ/2._rprec)

    end do
    end do

 end if !end spatial heat flux flag !end if block 102
  elseif (remote_flag==1) then !if block 101 conts.
     if (lbc .eq. 0) then

     ustar_local=ustar_avg

! ustar_local_mean=SUM(ustar_local(1:16,1:16))/float(nx*ny)
! use variable scalar_node_1 to store the scalar field at node 1
! test filter this variable to bring the surface flux calculated later in a local
! formulation more in line with the average formulation as prescribed by
! Monin-Obukhov theory

     scalar_node_1=scalar(:,:,1)

  if (remote_flux_homog_flag .eq. 1) then
     if (mod(jt,100)==0) print *,'apply MO in avg sense: remote_flux homog_flag ON'
     crap2=sum(scalar_node_1(1:nx,1:ny))/real(nx*ny);
     scalar_node_1(1:nx,1:ny)=crap2
     ustar_local=sum(sqrt(u(1:nx,1:ny,1)**2+v(1:nx,1:ny,1)**2))/float(nx*ny)*vonK/&
     (dlog(0.5_rprec*dz/exp(sum(dlog(zo(1:nx,1:ny)))/float(nx*ny)))-sum(psi_m(1:nx,1:ny))/float(nx*ny))
     if (mod(jt,1000)==0) then
        print *,'Testing remote flux homog calc:u*,theta_1,theta_sfc,psi_m,psi_h,phi_h,zo:(5:6,5)',&
        ustar_local(5:6,5),scalar_node_1(5:6,5),S_surf(5:6,5),psi_m(5:6,5),psi_h(5:6,5),phi_h(5:6,5),&
        zo(5:6,5)
     end if
  else
! Filtering only makes sense if the flux calculation is not based on mean values (remote_flux_homog_flga .ne. 1)

     call test_filter(scalar_node_1,G_test)
  end if

  do j=1,ny
  do i=1,nx
     surf_flux(i,j)=(S_Surf(i,j)-scalar_node_1(i,j))*vonk*ustar_local(i,j)&
     /(dlog(dz/(2._rprec*z_os(i,j)))-psi_h(i,j))
     dsdz(i,j,1) =-phi_h(i,j)*surf_flux(i,j)/(ustar_local(i,j)*vonk*DZ/2._rprec)
  end do
  end do

  if (jt==SCAL_init+1) then
    open(unit=56,file=path//'output/surf_flux_step1.txt',status="unknown",position="append")
    do j=1,ny
      write(56,5168) (surf_flux(i,j)*u_star*T_scale,i=1,nx)
    end do
    close(56)
  end if

5168     format(1400(E14.5))
  else if ((lbc .eq. 1) .and. (spatial_flux_flag) ) then
       if (jt .eq. SCAL_init) then
        print *,'reading spatial flux data from file !!'
        open(unit=77,file='./spatial_flux.dat',status='unknown')
          do j=1,ny
              read(77,5168) (surf_flux(i,j),i=1,nx)
          end do
          if (remote_to_patch_flag) then
             print *,'mean surf flux BEFORE remote to patch: ',sum(surf_flux)/float(nx*ny) 
             print *,'Creating 2 patch spatial heat flux field from spatial flux data ...'
             call remote_to_patch(surf_flux,1) 
             print *,'mean surf flux AFTER remote to patch: ',sum(surf_flux)/float(nx*ny)
          else if (remote_homog_flag == 1) then
             print *,'Homogenizing the spatial flux field as remote_homog_flag == 1'
             surf_flux=sum(surf_flux)/float(nx*ny)
          end if
          print *,'surf flux read',surf_flux(2:4,2:4) 
          surf_flux=surf_flux/u_star/T_scale !The readin surf_flux is dimensional - so non-dimensionalize
       else if(jt .GT. SCAL_init) then
          surf_flux=sgs_t3(1:nx,1:ny,1)
      end if !end for if loop for jt .eq. SCAL_init
       ustar_local=ustar_avg !set ustar as value computed in obukhov

        do j=1,ny
        do i=1,nx
          dsdz(i,j,1) =-phi_h(i,j)*surf_flux(i,j)/(ustar_local(i,j)*vonk*DZ/2._rprec) !set the gradient at the first point
        end do
        end do

     end if
  end if !end if block 101

end if !end for if (USE_MPI ...) block

 call dealias1(u,u_m)
 call dealias1(v,v_m)
 call dealias1(w,w_m)
 call dealias1(dsdx,dsdx_m)
 call dealias1(dsdy,dsdy_m)
 call dealias1(dsdz,dsdz_m)

 if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    jz_min = 2
 else
    jz_min = 1
 end if
 
    ubc_jz = nz-1

! Now compute the RHS term of the filtered scalar equation. 
! Note that this is the advection term with the scalar as 
! the diffusion term has been thrown away. This is done step 
! by step for each of the expanded arrays from dealias1 separately
! for the first node & last node AND the rest of the nodes.

! xxxxx ------Comments valid for MPI case only ---------XXXX
! For MPI case, all the nodes have fluid nodes (1:nz-1) except for
! near-the wall processes (2:nz-1 for w nodes and 1:nz-1 for uvp nodes)
! and the top nodes (1:nz)
! The following loop executes from 2:nz-1 for the near-wall process
! and 1:nz-1 for the other processes. Also, the subsequent 2 loops
! take care of the first node for the near-wall process (coord = 0)
! and the topmost node for the top process (coord = nproc-1).
! Note also that we need to make ghost data available for dTdz and
! w for the topmost node (jz=n) within each process and therefore,
! this synchronization (MPI_SENDRECV()) has been done in the subroutine
! theta_all_in_one ()
! xxxxx --------- MPI Comment block ends ------------------XXXX

do k=jz_min,nz-1
  do j=1,Ny2
    do i=1,Nx2
    RHS_m(i,j,k)=u_m(i,j,k)*dsdx_m(i,j,k)+v_m(i,j,k)*dsdy_m(i,j,k)&
    +(w_m(i,j,k)*dsdz_m(i,j,k)+w_m(i,j,k+1)*dsdz_m(i,j,k+1))/2._rprec
    end do
  end do
end do


 if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
   do j=1,Ny2
     do i=1,Nx2
      RHS_m(i,j,1)=u_m(i,j,1)*dsdx_m(i,j,1)+v_m(i,j,1)*dsdy_m(i,j,1)&
      +(0.5_rprec*w_m(i,j,2))*dsdz_m(i,j,2)
     end do
   end do
 end if
 
 if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
   do j=1,Ny2
     do i=1,Nx2
      RHS_m(i,j,Nz)=u_m(i,j,Nz)*dsdx_m(i,j,Nz)+v_m(i,j,Nz)*dsdy_m(i,j,Nz)
     end do
   end do
 end if

! call dealias2(RHS,RHS_m)

!MM--Adding extra terms for baroclinicity after dealiasing_May2016

do k=1,nz-1
  do j=1,Ny2
    do i=1,Nx2
 !       RHS_m(i,j,k)=RHS_m(i,j,k)+u_m(i,j,k)*((-0.005)*coriol*u_star)/g-v_m(i,j,k)*(0.0*coriol*u_star)/g
    end do
  end do
end do

 call dealias2(RHS,RHS_m)

!MM -----------

!c...Now building the SGS part of the RHS.
! Here the sgs_term for scalars is built up using Nu_t from sgs_stag_W.f
! and dividing it by the turbulent Prandtl # specified in dimen.h
!c....Note: Since we bring the Convective term to RHS its sign changes.
!c....Below "Temp" is used for SGS flux; its divergence is added to RHS
!VK.. Nu_t is on w nodes everywhere except at z=dz/2.. while
!VK dsdx is on uvp nodes.. so, interpolate Nu_t as we want temp to
!VK be on uvp nodes
! All this valid only till Pr_ is a constant..
! This will need a makeover once Pr_ becomes dynamic as well...

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
   do j=1,Ny
     do i=1,Nx
       temp(i,j,1)=(1./Pr_(i,j,1))*Nu_t(i,j,1)*dsdx(i,j,1)
        sgs_horx(i,j,1)=-1.*temp(i,j,1)
     end do
   end do
end if

  do k=jz_min,ubc_jz
    do j=1,Ny
      do i=1,Nx
        temp(i,j,k)=(1./(0.5_rprec*(Pr_(i,j,k)+Pr_(i,j,k+1))))* &
                0.5_rprec*(Nu_t(i,j,k)+Nu_t(i,j,k+1))*dsdx(i,j,k)
         sgs_horx(i,j,k)=-1.*temp(i,j,k)
        ! temp(i,j,k)=(1./Pr_(i,j,k))*0.5_rprec*(Nu_t(i,j,k)+Nu_t(i,j,k+1))*dsdx(i,j,k)
      end do
    end do
  end do


  call DDX (dtemp, temp)  

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
   do j=1,Ny
     do i=1,Nx
       RHS(i,j,1) = (-1.*RHS(i,j,1) + dtemp(i,j,1))
       temp(i,j,1)=(1./Pr_(i,j,1))*Nu_t(i,j,1)*dsdy(i,j,1)
       sgs_hory(i,j,1)=-1.*temp(i,j,1)
     end do
   end do
end if

 do k=jz_min,ubc_jz
! Nu_t is on w nodes and dsdy is on uvp nodes for jz=2 to nz
   do j=1,Ny
     do i=1,Nx
       RHS(i,j,k) = (-1.*RHS(i,j,k) + dtemp(i,j,k))
       temp(i,j,k)=(1./(0.5_rprec*(Pr_(i,j,k)+Pr_(i,j,k+1))))* &
               0.5_rprec*(Nu_t(i,j,k)+Nu_t(i,j,k+1))*dsdy(i,j,k)
       sgs_hory(i,j,k)=-1.*temp(i,j,k)
     end do
   end do
 end do


 call DDY (dtemp, temp)   
  
!c...Use MO flux at wall for the scalar sgs term !
!c Note that the total contribution to the scalar sgs term at
!c the first node comes from the surface flux computed above from
!c the specified heat flux, wt_s

 if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
   do j=1,Ny
     do i=1,Nx
      RHS(i,j,1) = RHS(i,j,1) + dtemp(i,j,1)
      temp(i,j,1) = -1.*surf_flux(i,j)
      sgs_vert(i,j,1) = surf_flux(i,j)
     end do
   end do
 end if
! Note sgs_vert is -1*temp because sgs_vert is modeled as -Nu_t*dsdz/Pr
! while temp is the same term but w/o the minus sign due to the additional
! minus outside the scalar fluctuation flux term in RHS
! need to run this loop nz due to the nature of the differenetiation in ddz_w

 do k=jz_min,nz
   do j=1,Ny
     do i=1,Nx
       RHS(i,j,k) = RHS(i,j,k) + dtemp(i,j,k)
       temp(i,j,k)=(1./Pr_(i,j,k))*Nu_t(i,j,k)*dsdz(i,j,k)
       sgs_vert(i,j,k)=-1.*temp(i,j,k)
     end do
   end do
 end do
!c...The SGS_z flux is on the W nodes, but DDZ_W will put it back on UVP nodes! 
!c Also note that sgs_vert(i,j,k) influences the computations in 
! OBUKHOV.f and is not involved in any computations in this routine.
! sgs_t3(i,j,1) (<w'theta'> is used for computing wt at the surface in OBUKHOV)

  call DDZ_w (dtemp, temp)
  do k=1,ubc_jz
    Do j=1,Ny
    do i=1,Nx
    RHS(i,j,k) = RHS(i,j,k) + dtemp(i,j,k)
    end do
    end do
  end do

end subroutine scalar_RHS_calc
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
subroutine step_scalar(scalar,RHS_pre,RHS_post)
implicit none
integer:: i,j,k
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif
real(kind=rprec),dimension(ld,ny,$lbz:nz)::scalar, RHS_pre, RHS_post
!real(kind=rprec)::wt_s_current
!cVK - This routine moves the scalar field (scalar in this case)
!cVK - forward in time using the scalar from previous time step
!cVK - and the RHS terms from the previous two time steps 
!cVK - using second order Adams-Bashforth scheme

! Note that in the last staments within this file, we set the value
! of scalar at the topmost node based on prescribed bc (inv_strength)
! and so, we could save some computation by only performing
! the scalar computation till Nz-1 global node...

do k=1,nz-1
      do j=1,ny
             do i=1,nx
                 scalar(i,j,k)= scalar(i,j,k)+dt*(1.5_rprec*RHS_pre(i,j,k)-0.5_rprec*RHS_post(i,j,k))
             end do
      end do
end do     

!VK Note that this subroutine was designed to be working with a set of scalars (incl.
!VK temperature and humidity and so, these boundary conditions as given below should
!VK be interpreted in the right context and not just for temperature
!VK For example, the first if block refers to a condition with humidity while the
!VK second and third statements are focussed to temperature

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
! if MPI - then clicks and is used for the process dealing wih the top nodes
! else in general is used for temp bc
! Note that L_z*nproc is the total vertical extent of the domain for the MPI and
! non-MPI cases ..(with MPI in place, we can not use L_z>z_i anymore)
     if ((L_z*nproc)>z_i) then ! for temperature and non-neutral case
         scalar(:,:,Nz)=scalar(:,:,Nz-1)+dTdz_top/T_scale*z_i*dz !ubc 
! inv_strength - refers to the slope of the inversion ~ 0.003K/Km for temperature
     else ! for everything else - neutral and passive scalars (may be modified depending on need)
         scalar(:,:,Nz)=scalar(:,:,Nz-1)
     end if
end if

end subroutine step_scalar

!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
subroutine obukhov
use types,only:rprec
use sim_param,only:u,v,theta,path
use bottombc !Includes patches subroutine, phi_m,psi_m,phi_h,psi_h,ustar_avg
use test_filtermodule
implicit none
integer:: jx,jy
real(kind=rprec), dimension(ld,ny):: wt_avg,wq_avg,theta_avg,u1,v1
real(kind=rprec), dimension(nx,ny):: x,zeta ! wstar, L already defined above
real(kind=rprec) g_,wt_,wq_,ustar_,theta_,L_,zo_,u_avg(nx,ny),fr,wstar_avg
real(kind=rprec),save:: obuk_L,obuk_ustar,obuk_phi_m,obuk_phi_h,obuk_psi_m   
real(kind=rprec),save:: obuk_wt_sfc,obuk_psi_h,obuk_zo,obuk_wstar   

if (jt .LT. SCAL_init .OR. (.NOT. S_FLAG)) then  
    if (.not. initsc) psi_m=0._rprec ! psi_m is being read in initial.f90
    phi_m=1._rprec
    psi_h=0._rprec
    phi_h=1._rprec
    ustar_avg=sum(sqrt(u(1:nx,1:ny,1)**2+v(1:nx,1:ny,1)**2))/float(nx*ny)*vonK/&
    (dlog(0.5_rprec*dz/zo)-sum(psi_m(1:nx,1:ny))/float(nx*ny))
    L=0._rprec
    wstar=0._rprec
    return
end if  

! Do the following only for jt .eq. SCAL_init
! This becomes curcial for the supercomputing runs as we need to break the
! simulation into smaller chunks and thereby need an accurate way to continue
! the simulation from the vel_sc.out file..therefore, added sgs_t3(:,:,1) i.e. 
! the surface flux to the list of output variables in io.f90
   if (jt .EQ. SCAL_init) then
    obuk_L=0._rprec
    obuk_ustar=0._rprec
    obuk_wstar=0._rprec
    obuk_phi_m=0._rprec
    obuk_phi_h=0._rprec
    obuk_psi_m=0._rprec
    obuk_psi_h=0._rprec
    obuk_zo=0._rprec
        if (.not. initsc) then
        sgs_t3(:,:,1)=wt_s/u_star/T_scale
        end if
   end if

!  nondimensionalize g
   g_=g/(u_star**2/z_i)  ! end if

theta_avg=theta(:,:,1) 
wt_avg=sgs_t3(:,:,1) ! We need only the surface flux - defined by sgs
zo_=exp(sum(dlog(zo(1:nx,1:ny)))/float(nx*ny))

! averages over x-y plane @ z = 1
wt_=sum(sgs_t3(1:nx,1:ny,1))/float(nx*ny)
ustar_=sum(sqrt(u(1:nx,1:ny,1)**2+v(1:nx,1:ny,1)**2))/float(nx*ny)*vonK/&
(dlog(0.5_rprec*dz/zo_)-sum(psi_m(1:nx,1:ny))/float(nx*ny))
theta_=sum(theta_avg(1:nx,1:ny))/float(nx*ny)

! SKS
! print*,jt,psi_m ! Turns out that psi_m hasnt been initialized and still its 0 at
! jt = 1
! SKS

if ((patch_flag==1 .and. num_patch==1) .OR. (OB_homog_flag)) then  
  do jx=1,nx
    do jy=1,ny
      wt_avg(jx,jy)=wt_
      ustar_avg(jx,jy)=ustar_
      theta_avg(jx,jy)=theta_
    end do
  end do
else
  u1=u(:,:,1)
  v1=v(:,:,1)
  call test_filter(u1,G_test)
  call test_filter(v1,G_test)
  do jx=1,nx
    do jy=1,ny
     u_avg(jx,jy)=sqrt(u1(jx,jy)**2+v1(jx,jy)**2)
    end do
  end do
  ustar_avg(1:nx,:)=u_avg(:,:)*vonK/(dlog(0.5_rprec*dz/zo(:,:))-psi_m(:,:))
end if

! Compute Obukhov Length
    
   do jx=1,nx
   do jy=1,ny
    L(jx,jy)=-ustar_avg(jx,jy)**3/(vonk*g_/theta_avg(jx,jy)*wt_avg(jx,jy))
! SKS
! Note that ustar_ that we calculate is non-dimensionalized by u_g = ustar = 8
! and theta_avg is non-dimensionalized by T_scale.
! If you work it out..i.e. if you substitute the formulae for each of the 
! quantities above, we get obukhov length normalized by z_i
! Note: g_=g/(u_star**2/z_i)
! SKS

! w_star is defined as [(g/<T_0>)*<w'T'>*z_i]**(1/3) (refer to 
! Nieuwstadt et al., Turbulent Shear flows, 1991)
! Therefore, for our case, where we are computing the non-dimensional w_star,
! the formula transforms to [(g_nd/<T_0_nd>)*<w'T'>_nd]**(1/3) where the suffix
! _nd refers to being non-dimensionalized using Ug (coriolis velocity), T_scale (300K)
! and z_i
    wstar(jx,jy)=sign((g_/theta_avg(jx,jy)*abs(wt_avg(jx,jy)))**(1./3.),wt_avg(jx,jy))

! wstar(jx,jy)=sign((g_/theta_avg(jx,jy)*abs(wt_avg(jx,jy))*z_i)**(1./3.),wt_avg(jx,jy))
! The above is the earlier wrong formula where there has been this additional z_i which for
! the post-processing means division by 10 as the usual value of z_i=1000 which with cube root on
! w_star just becomes a multiplication by 10. So, in case you feel that the value of w_star is quite
! high and the data is dated before Dec. 7th, 2004, please make sure to divide by 10

! for unstable conditions
      if ((L(jx,jy)<0._rprec) .and. (wt_avg(jx,jy) .ne. 0._rprec)) then
             x(jx,jy)=(1._rprec-16._rprec*dz/2._rprec/L(jx,jy))**.25_rprec
             psi_m(jx,jy)=2._rprec*dlog((1.+x(jx,jy))/2._rprec)+&
             dlog((1._rprec+x(jx,jy)**2)/2._rprec)-2._rprec*datan(x(jx,jy))+pi/2._rprec
             psi_h(jx,jy)=2._rprec*dlog((1._rprec+x(jx,jy)**2)/2._rprec)
             phi_m(jx,jy)=x(jx,jy)**(-1)
             phi_h(jx,jy)=x(jx,jy)**(-2)
      else if ((L(jx,jy)>0._rprec).and.(wt_avg(jx,jy).ne. 0._rprec)) then
! Implementing new formulations for phi and psi for stable case
! using Cheng & Brutsaert (2004): source - Brutsaert's book from
! Marc's Hydrology course..the new relations are from the GABLS study
             zeta(jx,jy)=0.5_rprec*dz/L(jx,jy)
! %%%%%%%%%%%%%%%%%%%%%%%% Formulations used in GABLS (Good for z/L < 1)%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%% Also apeearing in Hogstorm, BLM, 1988)%%%%%%%%%%%%%%%%
        if ((jan_diurnal_run) .OR. (GABLS_diurnal_test)) then !USE MO functions for z/L < 1
! %%%%%%%%%%%%%%%%%%%%%%%% Formulations used in Cheng & Brutsaert (BLM,2005) %%%%%%%%%%%%
! X=z/L where z in present case is dz/2 and L is the Obukhov length
! The equations: phi_m(X) = 1+a*[(X+(X^b)*{(1+X^b)^(-1+1/b)})/(X+{(1+X^b)^(1/b)})];a=6.1;b=2.5
! The equations: phi_h(X) = 1+c*[(X+(X^d)*{(1+X^d)^(-1+1/d)})/(X+{(1+X^d)^(1/d)})];c=5.3;d=1.1
! The equations: psi_m(X) = -a*ln[X+{(1+X^b)^(1/b)}];a=6.1;b=2.5
! The equctions: psi_h(X) = -c*ln[X+{(1+X^d)^(1/d)}];c=5.3;d=1.1
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        phi_m(jx,jy)=1._rprec+6.1_rprec*(zeta(jx,jy)+zeta(jx,jy)**2.5_rprec*&
        ((1._rprec+zeta(jx,jy)**2.5_rprec)**(-1._rprec+1._rprec/2.5_rprec)))/(zeta(jx,jy)+&
        ((1._rprec+zeta(jx,jy)**2.5_rprec)**(1._rprec/2.5_rprec)))
!             
        psi_m(jx,jy)=-1._rprec*6.1_rprec*dlog(zeta(jx,jy)+((1._rprec+&
        zeta(jx,jy)**2.5_rprec)**(1._rprec/2.5_rprec)))
!            
        phi_h(jx,jy)=1._rprec+5.3_rprec*(zeta(jx,jy)+zeta(jx,jy)**1.1_rprec*&
        ((1._rprec+zeta(jx,jy)**1.1_rprec)**(-1._rprec+1._rprec/1.1_rprec)))/(zeta(jx,jy)+&
        ((1._rprec+zeta(jx,jy)**1.1_rprec)**(1._rprec/1.1_rprec)))
!
        psi_h(jx,jy)=-1._rprec*5.3_rprec*dlog(zeta(jx,jy)+((1._rprec+&
        zeta(jx,jy)**1.1_rprec)**(1._rprec/1.1_rprec)))
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else

! SKS
! Flux profile relationships taken from paper by A.J.Dyer - "A Review of flux-profile
! relationships" - Boundary Layer Meteorology 7 (1974) 363-372
          phi_m(jx,jy)=1._rprec+5.0_rprec*zeta(jx,jy)
          phi_h(jx,jy)=1._rprec+5.0_rprec*zeta(jx,jy)
          psi_m(jx,jy)=-1._rprec*5.0_rprec*zeta(jx,jy)
          psi_h(jx,jy)=-1._rprec*5.0_rprec*zeta(jx,jy)

! The ones written below are from the original code
!         phi_m(jx,jy)=1._rprec+4.8_rprec*zeta(jx,jy)
!         phi_h(jx,jy)=1._rprec+7.8_rprec*zeta(jx,jy)
!         psi_m(jx,jy)=-1._rprec*4.8_rprec*zeta(jx,jy)
!         psi_h(jx,jy)=-1._rprec*7.8_rprec*zeta(jx,jy)
! SKS

! %%%%%%%%%%%%%%%%%%%%%%%% Formulations used in Brutsaert (1982) %%%%%%%%%%%%%%%%%%%%
!         phi_m(jx,jy)=1._rprec+5.0_rprec*zeta(jx,jy)
!         phi_h(jx,jy)=1._rprec+5.0_rprec*zeta(jx,jy)
!         psi_m(jx,jy)=-1._rprec*5.0_rprec*zeta(jx,jy)
!         psi_h(jx,jy)=-1._rprec*5.0_rprec*zeta(jx,jy)
        end if ! end for the jan_diurnal_run if block     
      else
             psi_m(jx,jy)=0._rprec
             psi_h(jx,jy)=0._rprec
             phi_m(jx,jy)=1._rprec
             phi_h(jx,jy)=1._rprec
      end if ! (Loop5 ends)
  end do
  end do

  L_=-(ustar_**3)/(vonk*(g_/theta_)*wt_)
  wstar_avg=sign((g_/theta_*abs(wt_))**(1./3.),wt_)

! SKS
  if (mod(jt,100) == 0) then     ! 100 because wbase = 100
! SKS
  write (6,7780) L_*z_i,ustar_*u_star,theta_*T_scale,(sum(T_s)/float(nx*ny))*T_scale,dz/2/L_,&
  wt_*u_star*T_scale
! SKS
  end if
! SKS

7780 format ('L(m),ustar(m/s),theta_1(K),T_s(K),z/L,wt_s(Km/s):',&
(5(1x,E12.6),1x,E12.6))

!-------------------- OUTPUT ------------------------------
! Output the heat flux time series to a file to be used later
    open (unit=47,file=path//'output/WT_sfc_tseries.out',status="unknown",position="append")
    write(47,5168) (jt_total+1)*dt,wt_
    close(47)
!----------------------------------------------------------

if (((initsc) .AND. ((jt_total+1) .ge. SCAL_init)) .OR. &
   ((.not. initsc) .AND. ((jt_total+1) .ge. SCAL_init+1))) then
  fr=1._rprec/c_count
  obuk_L=obuk_L+fr*L_
  obuk_ustar=obuk_ustar+fr*ustar_
  obuk_wstar=obuk_wstar+fr*wstar_avg
  obuk_phi_m=obuk_phi_m+fr*sum(phi_m(:,:))/float(nx*ny)
  obuk_psi_m=obuk_psi_m+fr*sum(psi_m(:,:))/float(nx*ny)
  obuk_phi_h=obuk_phi_h+fr*sum(phi_h(:,:))/float(nx*ny)
  obuk_psi_h=obuk_psi_h+fr*sum(psi_h(:,:))/float(nx*ny)
  obuk_zo=obuk_zo+fr*zo_
  obuk_wt_sfc=obuk_wt_sfc+fr*wt_

  if (mod(jt,c_count)==0) then
    open (unit=47,file=path//'output/mo.out',status="unknown",position="append")
    write(47,5168) (jt_total+1)*dt,obuk_L,obuk_ustar,obuk_wstar,obuk_phi_m,obuk_psi_m,&
                 obuk_phi_h,obuk_psi_h,obuk_zo,obuk_wt_sfc
    close(47)

    obuk_L=0._rprec;obuk_ustar=0._rprec;obuk_wstar=0._rprec;obuk_phi_m=0._rprec
    obuk_phi_h=0._rprec;obuk_psi_m=0._rprec;obuk_psi_h=0._rprec;obuk_zo=0._rprec
    obuk_wt_sfc=0._rprec;
  end if
end if

5168     format(E14.5,9(1x,E14.5))
  return
end subroutine obukhov 

end module scalars_module
