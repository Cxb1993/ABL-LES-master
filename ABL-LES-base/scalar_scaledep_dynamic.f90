subroutine scalar_scaledep_dynamic(Cs2_Pr,S11,S12,S13,S22,S23,S33)
! SKS
! Subroutine written by Stimit on 15th June 2010
! Returns Cs^2/Pr_sgs
! SKS
use types,only:rprec
use param
use sim_param,only:path,u,v,w,theta
use sgsmodule,only:Cs_opt2,Pr_t
use test_filtermodule
use bottombc

implicit none
integer::jz,tagPr=300
real(kind=rprec),dimension(ld,ny,nz)::S11,S12,S13,S22,S23,S33
real(kind=rprec),dimension(nz-1),intent(out)::Cs2_Pr
real(kind=rprec),dimension(nz-1)::Cs2_Pr4,Cs2_Pr2
! SKS
! 4 and 2 in the above variables (after Pr i.e. Pr4 and Pr2)
! correspond to the test and grid filters respectively
! SKS
real(kind=rprec),dimension(ld,ny)::L1,L2,L3
real(kind=rprec),dimension(ld,ny)::Q1,Q2,Q3
real(kind=rprec),dimension(ld,ny)::M1,M2,M3
real(kind=rprec),dimension(ld,ny)::N1,N2,N3
real(kind=rprec),dimension(ld,ny)::S_bar,S11_bar,S12_bar, &
     S13_bar,S22_bar,S23_bar,S33_bar
real(kind=rprec),dimension(ld,ny)::S_hat,S11_hat,S12_hat, &
     S13_hat,S22_hat,S23_hat,S33_hat
real(kind=rprec),dimension(ld,ny)::S_dTdx_bar,S_dTdy_bar,S_dTdz_bar
real(kind=rprec),dimension(ld,ny)::S_dTdx_hat,S_dTdy_hat,S_dTdz_hat
real(kind=rprec),dimension(ld,ny)::dTdx_bar,dTdy_bar,dTdz_bar, &
                                   dTdx_hat,dTdy_hat,dTdz_hat

real(kind=rprec),dimension(ld,ny)::u_bar,v_bar,w_bar,theta_bar
real(kind=rprec),dimension(ld,ny)::u_hat,v_hat,w_hat,theta_hat
real(kind=rprec),dimension(ld,ny)::S
real(kind=rprec),dimension(nz-1)::beta,avgCs2,avgCs,avgPr_t
real(kind=rprec),dimension(nz_tot-1)::Cs2_Pr_tot,Cs2_Pr4_tot, &
                                      Cs2_Pr2_tot,avgCs_tot, &
                                      avgPr_t_tot

real(kind=rprec)::delta,const
real(kind=rprec)::betaclip
character(len=24)::fname
real(kind=rprec)::z,surface_flux,phi_h_avg,zo_avg,ustar_

$if ($MPI)
  $define $lbz 0
  integer::recvcounts(nproc)
  integer::displs(nproc)
$else
  $define $lbz 1
$endif
! Following definition is below all others because it involves $lbz
real(kind=rprec),dimension(ld,ny,$lbz:nz)::dTdx,dTdy,dTdz

delta = filter_size*(dx*dy*dz)**(1._rprec/3._rprec)

call filt_da(theta,dTdx,dTdy)
call ddz_uv(dTdz,theta)
if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
  dTdz(:,:,nz)=dTdz_top/T_scale*z_i ! Valid for temperature
end if

! ** IMP : The model only works for prescribed surface flux i.e. lbc = 1
! ** and also only for homogeneous boundary conditions (at the surface)
if (lbc .ne. 1) then
   print*,'Please check the bottom boundary condition'
   stop
end if

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  zo_avg       = exp(sum(dlog(zo(1:nx,1:ny)))/float(nx*ny))
  ustar_       = sum(sqrt(u(1:nx,1:ny,1)**2+v(1:nx,1:ny,1)**2))/float(nx*ny)*vonk/ &
                 (dlog(0.5_rprec*dz/zo_avg)-sum(psi_m(1:nx,1:ny))/float(nx*ny))
  phi_h_avg    = sum(phi_h(1:nx,1:ny))/float(nx*ny)
  surface_flux = wt_s/T_scale/u_star
  dTdz(:,:,1)  = phi_h_avg*surface_flux/(ustar_*vonk*dz/2._rprec)
end if

do jz=1,nz-1
!  using L_j as temp storage here
   if ( ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) .and. &
        (jz == 1) ) then
     ! watch the 0.25's: w = c*z*z near wall, so get 0.25
     L1(:,:) = u(:,:,1)*theta(:,:,1)             ! uvptheta-node
     L2(:,:) = v(:,:,1)*theta(:,:,1)             ! uvptheta-node
     L3(:,:) = 0.25_rprec*w(:,:,2)*theta(:,:,1)  ! uvptheta-node
     u_bar(:,:) = u(:,:,1)
     v_bar(:,:) = v(:,:,1)
     w_bar(:,:) = 0.25_rprec*w(:,:,2)
     theta_bar(:,:) = theta(:,:,1)
   else  ! w-nodes
     L1(:,:) = 0.5_rprec*(u(:,:,jz) + u(:,:,jz-1))* &
               0.5_rprec*(theta(:,:,jz) + theta(:,:,jz-1))
     L2(:,:) = 0.5_rprec*(v(:,:,jz) + v(:,:,jz-1))* &
               0.5_rprec*(theta(:,:,jz) + theta(:,:,jz-1))
     L3(:,:) = w(:,:,jz)* &
               0.5_rprec*(theta(:,:,jz) + theta(:,:,jz-1))
     u_bar(:,:) = 0.5_rprec*(u(:,:,jz) + u(:,:,jz-1))
     v_bar(:,:) = 0.5_rprec*(v(:,:,jz) + v(:,:,jz-1))
     w_bar(:,:) = w(:,:,jz)
     theta_bar(:,:) = 0.5_rprec*(theta(:,:,jz) + theta(:,:,jz-1))
   end if
   u_hat     = u_bar
   v_hat     = v_bar
   w_hat     = w_bar
   theta_hat = theta_bar
   
   call test_filter(u_bar,G_test)
   call test_filter(v_bar,G_test)
   call test_filter(w_bar,G_test)
   call test_filter(theta_bar,G_test)

   call test_filter(L1,G_test) ! in-place filtering
   L1 = L1 - u_bar*theta_bar
   call test_filter(L2,G_test)
   L2 = L2 - v_bar*theta_bar
   call test_filter(L3,G_test)
   L3 = L3 - w_bar*theta_bar
   
   Q1 = u_bar*theta_bar
   Q2 = v_bar*theta_bar
   Q3 = w_bar*theta_bar

   call test_filter(u_hat,G_test_test)
   call test_filter(v_hat,G_test_test)
   call test_filter(w_hat,G_test_test)
   call test_filter(theta_hat,G_test_test)

   call test_filter(Q1,G_test_test)
   Q1 = Q1 - u_hat*theta_hat
   call test_filter(Q2,G_test_test)
   Q2 = Q2 - v_hat*theta_hat
   call test_filter(Q3,G_test_test)
   Q3 = Q3 - w_hat*theta_hat

!  calculate |S|
   S(:,:) = sqrt(2._rprec*(S11(:,:,jz)**2 + S22(:,:,jz)**2 +& 
                 S33(:,:,jz)**2 + 2._rprec*(S12(:,:,jz)**2 +&
                 S13(:,:,jz)**2 + S23(:,:,jz)**2)))

!  already on w-nodes
   S11_bar(:,:) = S11(:,:,jz)  
   S12_bar(:,:) = S12(:,:,jz)  
   S13_bar(:,:) = S13(:,:,jz)  
   S22_bar(:,:) = S22(:,:,jz)  
   S23_bar(:,:) = S23(:,:,jz)  
   S33_bar(:,:) = S33(:,:,jz)  
   
   S11_hat = S11_bar
   S12_hat = S12_bar
   S13_hat = S13_bar
   S22_hat = S22_bar
   S23_hat = S23_bar
   S33_hat = S33_bar
   
   call test_filter(S11_bar,G_test)
   call test_filter(S12_bar,G_test)
   call test_filter(S13_bar,G_test)
   call test_filter(S22_bar,G_test)
   call test_filter(S23_bar,G_test)
   call test_filter(S33_bar,G_test)

   call test_filter(S11_hat,G_test_test)
   call test_filter(S12_hat,G_test_test)
   call test_filter(S13_hat,G_test_test)
   call test_filter(S22_hat,G_test_test)
   call test_filter(S23_hat,G_test_test)
   call test_filter(S33_hat,G_test_test)

   S_bar = sqrt(2._rprec*(S11_bar**2 + S22_bar**2 + S33_bar**2 + &
                2._rprec*(S12_bar**2 + S13_bar**2 + S23_bar**2)))

   S_hat = sqrt(2._rprec*(S11_hat**2 + S22_hat**2 + S33_hat**2 + &
                2._rprec*(S12_hat**2 + S13_hat**2 + S23_hat**2)))

   S_dTdx_bar = S*dTdx(:,:,jz)
   S_dTdy_bar = S*dTdy(:,:,jz)
   S_dTdz_bar = S*dTdz(:,:,jz)

   S_dTdx_hat = S_dTdx_bar
   S_dTdy_hat = S_dTdy_bar
   S_dTdz_hat = S_dTdz_bar

   call test_filter(S_dTdx_bar,G_test)
   call test_filter(S_dTdy_bar,G_test)
   call test_filter(S_dTdz_bar,G_test)

   call test_filter(S_dTdx_hat,G_test_test)
   call test_filter(S_dTdy_hat,G_test_test)
   call test_filter(S_dTdz_hat,G_test_test)

   dTdx_bar   = dTdx(:,:,jz)
   dTdy_bar   = dTdy(:,:,jz)
   dTdz_bar   = dTdz(:,:,jz)

   dTdx_hat   = dTdx_bar
   dTdy_hat   = dTdy_bar
   dTdz_hat   = dTdz_bar

   call test_filter(dTdx_bar,G_test)
   call test_filter(dTdy_bar,G_test)
   call test_filter(dTdz_bar,G_test)

   call test_filter(dTdx_hat,G_test_test)
   call test_filter(dTdy_hat,G_test_test)
   call test_filter(dTdz_hat,G_test_test)

   const = delta**2
   M1 = const*(S_dTdx_bar - 4._rprec*S_bar*dTdx_bar)
   M2 = const*(S_dTdy_bar - 4._rprec*S_bar*dTdy_bar)
   M3 = const*(S_dTdz_bar - 4._rprec*S_bar*dTdz_bar)

   N1 = const*(S_dTdx_hat - 16._rprec*S_hat*dTdx_hat)
   N2 = const*(S_dTdy_hat - 16._rprec*S_hat*dTdy_hat)
   N3 = const*(S_dTdz_hat - 16._rprec*S_hat*dTdz_hat)
 
   ! M1 = max(real(M1),real(1E-24))
   ! M2 = max(real(M2),real(1E-24))
   ! M3 = max(real(M3),real(1E-24))

   ! N1 = max(real(N1),real(1E-24))
   ! N2 = max(real(N2),real(1E-24))
   ! N3 = max(real(N3),real(1E-24))

   avgCs(jz) = sum(sqrt(Cs_opt2(1:nx,1:ny,jz)))/(nx*ny)

   Cs2_Pr2(jz) = sum(L1*M1 + L2*M2 + L3*M3)/sum(M1**2 + M2**2 + M3**2)
   Cs2_Pr4(jz) = sum(Q1*N1 + Q2*N2 + Q3*N3)/sum(N1**2 + N2**2 + N3**2)

   Cs2_Pr2(jz) = max(real(1E-24),real(Cs2_Pr2(jz)))
   Cs2_Pr4(jz) = max(real(1E-24),real(Cs2_Pr4(jz)))

   beta(jz)    = Cs2_Pr4(jz)/Cs2_Pr2(jz)
   betaclip    = max(real(beta(jz)),real(1._rprec/8._rprec))

   Cs2_Pr(jz)  = Cs2_Pr2(jz)/betaclip
   Cs2_Pr(jz)  = max(0._rprec,real(Cs2_Pr(jz),kind=rprec))

   Pr_t(:,:,jz) = Cs_opt2(:,:,jz)/Cs2_Pr(jz)
   avgPr_t(jz)  = sum(Pr_t(1:nx,1:ny,jz))/(nx*ny)
end do

$if ($MPI)
  call mpi_sendrecv (Pr_t(1,1,1), ld*ny,MPI_RPREC,down,tagPr,  &
                     Pr_t(1,1,nz),ld*ny,MPI_RPREC,up,  tagPr,  &
                     comm,status,ierr)
  call mpi_sendrecv (Pr_t(1,1,nz-1),ld*ny,MPI_RPREC,up,  tagPr+1, &
                     Pr_t(1,1,0),   ld*ny,MPI_RPREC,down,tagPr+1, &
                     comm,status,ierr)
$endif


$if ($MPI)
   recvcounts = size(avgPr_t)
   displs = coord_of_rank * recvcounts
   call mpi_gatherv (avgPr_t(1),size(avgPr_t),MPI_RPREC,avgPr_t_tot(1), &
                     recvcounts,displs,MPI_RPREC,rank_of_coord(0),comm,ierr)

   recvcounts = size(avgCs)
   displs = coord_of_rank * recvcounts
   call mpi_gatherv (avgCs(1),size(avgCs),MPI_RPREC,avgCs_tot(1), &
                     recvcounts,displs,MPI_RPREC,rank_of_coord(0),comm,ierr)

   recvcounts = size(Cs2_Pr)
   displs = coord_of_rank * recvcounts
   call mpi_gatherv (Cs2_Pr(1),size(Cs2_Pr),MPI_RPREC,Cs2_Pr_tot(1), &
                     recvcounts,displs,MPI_RPREC,rank_of_coord(0),comm,ierr)

   recvcounts = size(Cs2_Pr2)
   displs = coord_of_rank * recvcounts
   call mpi_gatherv (Cs2_Pr2(1),size(Cs2_Pr2),MPI_RPREC,Cs2_Pr2_tot(1), &
                     recvcounts,displs,MPI_RPREC,rank_of_coord(0),comm,ierr)

   recvcounts = size(Cs2_Pr4)
   displs = coord_of_rank * recvcounts
   call mpi_gatherv (Cs2_Pr4(1),size(Cs2_Pr4),MPI_RPREC,Cs2_Pr4_tot(1), &
                     recvcounts,displs,MPI_RPREC,rank_of_coord(0),comm,ierr)
$else
   avgCs_tot   = avgCs
   avgPr_t_tot = avgPr_t
   Cs2_Pr_tot  = Cs2_Pr
   Cs2_Pr2_tot = Cs2_Pr2
   Cs2_Pr4_tot = Cs2_Pr4
$endif

if (mod(jt,p_count) .eq. 0) then
   if((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
      do jz=1,nz_tot-1
         if (jz .eq. 1) then
            z = 0.5_rprec*dz
         else
            z = (jz-1)*dz
         end if
         write(14,7780)(jt_total-1)*dt,z,Cs2_Pr4_tot(jz)/Cs2_Pr2_tot(jz), &
            Cs2_Pr2_tot(jz),Cs2_Pr4_tot(jz),Cs2_Pr_tot(jz),avgCs_tot(jz), &
            avgPr_t_tot(jz)
      end do
      write(14,*)
   end if
end if

7780 format(1400(E14.5))

end subroutine scalar_scaledep_dynamic
