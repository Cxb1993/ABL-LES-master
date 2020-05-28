subroutine lagrange_scalar_Sdep(Cs2_Pr,S11,S12,S13,S22,S23,S33)
! Subroutine written by Stimit Shah on 27th Sept 2010
! Calculates Cs2byPr at delta with averaging along fluid pathlines.
use types,only:rprec
use param
use sim_param,only:u,v,w,theta
use sgsmodule,only:G_LM,G_MM,G_QN,G_NN,Cs_opt2,Pr_t,opftime, &
                   Beta,Beta_avg,Betaclip_avg,TnLMMM,TnQNNN
use test_filtermodule
use bottombc

implicit none

integer::jx,jy,jz
integer::counter
integer::istart,iend

integer::tagPr=320
real(kind=rprec),dimension(ld,ny,nz)::S11,S12,S13,S22,S23,S33
real(kind=rprec),dimension(ld,ny,nz),intent(out)::Cs2_Pr
real(kind=rprec),dimension(ld,ny,nz)::Cs2_Pr4,Cs2_Pr2

real(kind=rprec)::tf1,tf2,tf1_2,tf2_2
real(kind=rprec)::fractus
real(kind=rprec)::Betaclip

real(kind=rprec),dimension(ld,ny)::S,tempos
real(kind=rprec),dimension(ld,ny)::L1,L2,L3,M1,M2,M3
real(kind=rprec),dimension(ld,ny)::Q1,Q2,Q3,N1,N2,N3
real(kind=rprec),dimension(ld,ny)::LM,MM,QN,NN,Tn,epsi,dumfac

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

real(kind=rprec),dimension(nz)::avgPr_t,avgCs

real(kind=rprec),dimension(nz-1)::dumavgPr_t,dumCs2_Pr,dumCs2_Pr2, &
                                  dumCs2_Pr4,dumavgCs
real(kind=rprec),dimension(nz_tot-1)::Cs2_Pr_tot,Cs2_Pr4_tot, &
                                      Cs2_Pr2_tot,avgCs_tot,  &
                                      avgPr_t_tot
real(kind=rprec)::delta,const
real(kind=rprec)::lagran_dt,opftdelta,powcoeff
character(len=24)::fname
real(kind=rprec)::z,zo_avg,ustar_
real(kind=rprec),dimension(nx,ny)::surface_flux

$if ($MPI)
  $define $lbz 0
  integer::recvcounts(nproc)
  integer::displs(nproc)
$else
  $define $lbz 1
$endif

! Following definitions are below all others because it involves $lbz
real(kind=rprec),dimension(ld,ny,$lbz:nz)::dTdx,dTdy,dTdz

real(rprec),dimension(ld,ny,$lbz:nz)::u_temp,v_temp,theta_temp
logical,save::G_LM_MM_init=.false.
logical,save::G_QN_NN_init=.false.

delta = filter_size*(dx*dy*dz)**(1._rprec/3._rprec)

tf1=2._rprec
tf2=4._rprec
tf1_2=tf1**2
tf2_2=tf2**2

opftdelta = opftime*delta
powcoeff  = -1._rprec/8._rprec

u_temp     = u
v_temp     = v
theta_temp = theta

call interpolag_scalar_Sdep()

! this is the time step used in the lagrangian computations
lagran_dt = dt*real(cs_count,kind=rprec)
fractus   = 1._rprec/real(ny*nx,kind=rprec)
const     = delta**2

call filt_da(theta,dTdx,dTdy)
call ddz_uv(dTdz,theta)
if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
  dTdz(:,:,nz)=dTdz_top/T_scale*z_i ! Valid for temperature
end if

!==============================================================================
! ** IMP : The model works only for homogeneous BCs (at the surface)
! ** IMP : Should be checked wehn running with heterogeneities
!==============================================================================
! Compute surface flux and dTdz at z = dz/2
! lbc = 1 is used for prescribe the surface flux while
! lbc = 0 has been used to prescribe the temperature
if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
     if (lbc == 1) then
        zo_avg = exp(sum(dlog(zo(1:nx,1:ny)))/float(nx*ny))
        ustar_ = sum(sqrt(u(1:nx,1:ny,1)**2+v(1:nx,1:ny,1)**2))/float(nx*ny)*vonk/ &
                 (dlog(0.5_rprec*dz/zo_avg)-sum(psi_m(1:nx,1:ny))/float(nx*ny))
        do jy=1,ny
        do jx=1,nx
           surface_flux(jx,jy) = wt_s/T_scale/u_star
        end do
        end do
        ! The u_star above is coming from param.f90 = ug for coriolis and isnt
        ! the local u_star computed from stress at the surface.
     else if (lbc == 0) then
        zo_avg = exp(sum(dlog(zo(1:nx,1:ny)))/float(nx*ny))
        ustar_ = sum(sqrt(u(1:nx,1:ny,1)**2+v(1:nx,1:ny,1)**2))/float(nx*ny)*vonk/ &
                 (dlog(0.5_rprec*dz/zo_avg)-sum(psi_m(1:nx,1:ny))/float(nx*ny))
        do jy=1,ny
        do jx=1,nx
           surface_flux(jx,jy) = (T_s(jx,jy) - theta(jx,jy,1))*vonk*ustar_ &
           /(dlog(dz/(2._rprec*z_os(jx,jy)))-psi_h(jx,jy))
        end do
        end do
     end if
     ! Now we have the lowest dTdz on the uvp nodes all others on w nodes
   do jy=1,ny
   do jx=1,nx
     dTdz(jx,jy,1) = -phi_h(jx,jy)*surface_flux(jx,jy)/(ustar_*vonk*dz/2._rprec)
   end do
   end do
end if
!==============================================================================

do jz=1,nz

  if ( ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) .and. (jz == 1) ) then
    ! watch the 0.25's : recall w = c*z^2 close to wall, so get 0.25
    u_bar(:,:)     = u_temp(:,:,1)       ! first uvtheat-node
    v_bar(:,:)     = v_temp(:,:,1)
    w_bar(:,:)     = 0.25_rprec*w(:,:,2)
    theta_bar(:,:) = theta_temp(:,:,1)
  else ! w-nodes
    u_bar(:,:)     = 0.5_rprec*(u_temp(:,:,jz) + u_temp(:,:,jz-1))
    v_bar(:,:)     = 0.5_rprec*(v_temp(:,:,jz) + v_temp(:,:,jz-1))
    w_bar(:,:)     = w(:,:,jz)
    theta_bar(:,:) = 0.5_rprec*(theta_temp(:,:,jz) + theta_temp(:,:,jz-1))
  end if
  u_hat     = u_bar
  v_hat     = v_bar
  w_hat     = w_bar
  theta_hat = theta_bar

  L1 = u_bar*theta_bar
  L2 = v_bar*theta_bar
  L3 = w_bar*theta_bar

  Q1 = u_bar*theta_bar
  Q2 = v_bar*theta_bar
  Q3 = w_bar*theta_bar

  call test_filter(u_bar,G_test)
  call test_filter(v_bar,G_test)
  call test_filter(w_bar,G_test)
  call test_filter(theta_bar,G_test)

  call test_filter(L1,G_test)  ! in-place filtering
  L1 = L1 - u_bar*theta_bar
  call test_filter(L2,G_test)
  L2 = L2 - v_bar*theta_bar
  call test_filter(L3,G_test)
  L3 = L3 - w_bar*theta_bar

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

  ! calculate |S|
  S(:,:) = sqrt(2._rprec*(S11(:,:,jz)**2 + S22(:,:,jz)**2 + S33(:,:,jz)**2 +&
                2._rprec*(S12(:,:,jz)**2 + S13(:,:,jz)**2 + S23(:,:,jz)**2)))

  ! S_ij already on w-nodes
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
  
  S_bar = sqrt(2._rprec*(S11_bar**2 + S22_bar**2 + S33_bar**2 +&
               2._rprec*(S12_bar**2 + S13_bar**2 + S23_bar**2)))
  
  S_hat = sqrt(2._rprec*(S11_hat**2 + S22_hat**2 + S33_hat**2 +&
               2._rprec*(S12_hat**2 + S13_hat**2 + S23_hat**2)))

  if (((.not. USE_MPI) .or. (USE_MPI.and.coord==0)) .and. (jz==1)) then
     S_dTdx_bar = S*dTdx(:,:,jz)
     S_dTdy_bar = S*dTdy(:,:,jz)
     S_dTdz_bar = S*(dTdz(:,:,jz) + dTdz(:,:,jz+1))/2._rprec
  else
     S_dTdx_bar = S(:,:)*(dTdx(:,:,jz) + dTdx(:,:,jz-1))/2._rprec
     S_dTdy_bar = S(:,:)*(dTdy(:,:,jz) + dTdy(:,:,jz-1))/2._rprec
     S_dTdz_bar = S(:,:)*dTdz(:,:,jz)
  end if

  S_dTdx_hat = S_dTdx_bar
  S_dTdy_hat = S_dTdy_bar
  S_dTdz_hat = S_dTdz_bar

  call test_filter(S_dTdx_bar,G_test)
  call test_filter(S_dTdy_bar,G_test)
  call test_filter(S_dTdz_bar,G_test)

  call test_filter(S_dTdx_hat,G_test_test)
  call test_filter(S_dTdy_hat,G_test_test)
  call test_filter(S_dTdz_hat,G_test_test)
  
  if (((.not. USE_MPI) .or. (USE_MPI.and.coord==0)) .and. (jz==1)) then
     dTdx_bar   = dTdx(:,:,jz)
     dTdy_bar   = dTdy(:,:,jz)
     dTdz_bar   = (dTdz(:,:,jz) + dTdz(:,:,jz+1))/2._rprec
  else
     dTdx_bar = (dTdx(:,:,jz) + dTdx(:,:,jz-1))/2._rprec
     dTdy_bar = (dTdy(:,:,jz) + dTdy(:,:,jz-1))/2._rprec
     dTdz_bar = dTdz(:,:,jz)
  end if

  dTdx_hat = dTdx_bar
  dTdy_hat = dTdy_bar
  dTdz_hat = dTdz_bar

  call test_filter(dTdx_bar,G_test)
  call test_filter(dTdy_bar,G_test)
  call test_filter(dTdz_bar,G_test)

  call test_filter(dTdx_hat,G_test_test)
  call test_filter(dTdy_hat,G_test_test)
  call test_filter(dTdz_hat,G_test_test)
  
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

  LM = L1*M1 + L2*M2 + L3*M3
  MM = M1**2 + M2**2 + M3**2
  QN = Q1*N1 + Q2*N2 + Q3*N3
  NN = N1**2 + N2**2 + N3**2

  if (inilag) then
    if ((.not. G_LM_MM_init) .and. (jt == cs_count .or. jt == DYN_init)) then
      print *,'G_MM and G_LM initialized' 
      G_MM (:,:,jz) = MM
      G_LM (:,:,jz) = 0.03_rprec*MM
      G_MM(ld-1:ld,:,jz)=1._rprec
      G_LM(ld-1:ld,:,jz)=1._rprec

      if (jz == nz) G_LM_MM_init = .true.
    end if
  end if

  if(inflow)then
    iend   = floor(buff_end*nx + 1._rprec)
    iend   = modulo(iend-1,nx) + 1
    istart = floor((buff_end-buff_len)*nx + 1._rprec)
    istart = modulo(istart-1,nx) + 1
      
    Tn = merge(.1_rprec*const*S**2,MM,MM.le..1_rprec*const*S**2)
    MM = Tn
    LM(istart+1:iend,1:ny)      = 0._rprec
    G_LM(istart+1:iend,1:ny,jz) = 0._rprec
    Tn = merge(0.1_rprec*const*S**2,NN,NN.le.0.1_rprec*const*S**2)
    NN = Tn
    QN(istart+1:iend,1:ny)      = 0._rprec
    G_QN(istart+1:iend,1:ny,jz) = 0._rprec
  endif

  ! Tn = max(real(G_LM(:,:,jz)*G_MM(:,:,jz)),real(1E-24))
  Tn(:,:) = TnLMMM(:,:,jz)
  Tn = opftdelta*(Tn**powcoeff)
  Tn(:,:) = max(real(1E-24),real(Tn(:,:)))
  dumfac = lagran_dt/Tn 
  epsi = dumfac/(1._rprec + dumfac)
  
  G_LM(:,:,jz)=(epsi*LM + (1._rprec-epsi)*G_LM(:,:,jz))
  G_MM(:,:,jz)=(epsi*MM + (1._rprec-epsi)*G_MM(:,:,jz))
  G_LM(:,:,jz)= max(real(1E-24),real(G_LM(:,:,jz)))
  
  Cs2_Pr2(:,:,jz)    = G_LM(:,:,jz)/G_MM(:,:,jz)
  Cs2_Pr2(ld,:,jz)   = 1E-24
  Cs2_Pr2(ld-1,:,jz) = 1E-24
  Cs2_Pr2(:,:,jz)    = max(real(1E-24),real(Cs2_Pr2(:,:,jz)))

  if (inilag) then
    if ((.not. G_QN_NN_init) .and. (jt == cs_count .or. jt == DYN_init)) then
      print *,'G_NN and G_QN initialized'
      G_NN(:,:,jz) = NN
      G_QN(:,:,jz) = 0.03_rprec*NN
      G_NN(ld-1:ld,:,jz) = 1._rprec
      G_QN(ld-1:ld,:,jz) = 1._rprec
  
      if (jz == nz) G_QN_NN_init = .true.
    end if
  end if

  ! Tn = max(real(G_QN(:,:,jz)*G_NN(:,:,jz)),real(1E-24))
  Tn(:,:) = TnQNNN(:,:,jz)
  Tn = opftdelta*(Tn**powcoeff)
  Tn(:,:) = max(real(1E-24),real(Tn(:,:)))
  dumfac = lagran_dt/Tn
  epsi = dumfac/(1._rprec+dumfac)
  
  G_QN(:,:,jz) = (epsi*QN + (1._rprec-epsi)*G_QN(:,:,jz))
  G_NN(:,:,jz) = (epsi*NN + (1._rprec-epsi)*G_NN(:,:,jz))
  G_QN(:,:,jz) =  max(real(1E-24),real(G_QN(:,:,jz)))
  
  Cs2_Pr4(:,:,jz)    = G_QN(:,:,jz)/G_NN(:,:,jz)
  Cs2_Pr4(ld,:,jz)   = 1E-24
  Cs2_Pr4(ld-1,:,jz) = 1E-24
  Cs2_Pr4(:,:,jz)    = max(real(1E-24),real(Cs2_Pr4(:,:,jz)))
  
  Beta(:,:,jz) = (Cs2_Pr4(:,:,jz)/Cs2_Pr2(:,:,jz))**(log(tf1)/(log(tf2)-log(tf1)))

  counter = 0      
  do jy=1,Ny
    do jx=1,Nx
      if (Beta(jx,jy,jz).le. 1/(tf1*tf2)) counter = counter+1
    end do
  end do

  if ( ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) .and.  &
       (jz == nz) ) then
       Beta(:,:,jz)=1._rprec
  end if

  do jy = 1,ny
    do jx = 1,ld
      Betaclip=max(real(Beta(jx,jy,jz)),real(1.0_rprec/(tf1*tf2)))
      Cs2_Pr(jx,jy,jz)=Cs2_Pr2(jx,jy,jz)/Betaclip
    end do
  end do
  
  Cs2_Pr(ld,:,jz)   = 1E-24
  Cs2_Pr(ld-1,:,jz) = 1E-24
  Cs2_Pr(:,:,jz)    = max(real(1E-24),real(Cs2_Pr(:,:,jz)))

  Pr_t(:,:,jz) = Cs_opt2(:,:,jz)/Cs2_Pr(:,:,jz)
  avgPr_t(jz)  = sum(Pr_t(1:nx,1:ny,jz))/(nx*ny)
  avgCs(jz)    = sum(sqrt(Cs_opt2(1:nx,1:ny,jz)))/(nx*ny)
  
  ! Beta_avg(jz)     = sum(Cs2_Pr2(1:nx,1:ny,jz))/sum(Cs2_Pr(1:nx,1:ny,jz))
  ! Betaclip_avg(jz) = sum(Cs2_Pr4(1:nx,1:ny,jz))/sum(Cs2_Pr2(1:nx,1:ny,jz))

end do  ! this ends the main jz=1-nz loop          

$if ($MPI)
  call mpi_sendrecv (Pr_t(1,1,1), ld*ny,MPI_RPREC,down,tagPr,  &
                     Pr_t(1,1,nz),ld*ny,MPI_RPREC,up,  tagPr,  &
                     comm,status,ierr)
  call mpi_sendrecv (Pr_t(1,1,nz-1),ld*ny,MPI_RPREC,up,  tagPr+1, &
                     Pr_t(1,1,0),   ld*ny,MPI_RPREC,down,tagPr+1, &
                     comm,status,ierr)
$endif

do jz=1,nz-1
   dumavgPr_t(jz) = avgPr_t(jz)
   dumCs2_Pr4(jz) = sum(Cs2_Pr4(1:nx,1:ny,jz))/(nx*ny)
   dumCs2_Pr2(jz) = sum(Cs2_Pr2(1:nx,1:ny,jz))/(nx*ny)
   dumCs2_Pr(jz)  = sum(Cs2_Pr(1:nx,1:ny,jz))/(nx*ny)
   dumavgCs(jz)   = avgCs(jz)
end do

$if ($MPI)
   recvcounts = size(dumavgCs)
   displs = coord_of_rank * recvcounts
   call mpi_gatherv (dumavgCs(1),size(dumavgCs),MPI_RPREC,avgCs_tot(1), &
                     recvcounts,displs,MPI_RPREC,rank_of_coord(0),comm,ierr)

   recvcounts = size(dumavgPr_t)
   displs = coord_of_rank * recvcounts
   call mpi_gatherv (dumavgPr_t(1),size(dumavgPr_t),MPI_RPREC,avgPr_t_tot(1), &
                     recvcounts,displs,MPI_RPREC,rank_of_coord(0),comm,ierr)

   recvcounts = size(dumCs2_Pr)
   displs = coord_of_rank * recvcounts
   call mpi_gatherv (dumCs2_Pr(1),size(dumCs2_Pr),MPI_RPREC,Cs2_Pr_tot(1), &
                     recvcounts,displs,MPI_RPREC,rank_of_coord(0),comm,ierr)

   recvcounts = size(dumCs2_Pr2)
   displs = coord_of_rank * recvcounts
   call mpi_gatherv (dumCs2_Pr2(1),size(dumCs2_Pr2),MPI_RPREC,Cs2_Pr2_tot(1), &
                     recvcounts,displs,MPI_RPREC,rank_of_coord(0),comm,ierr)

   recvcounts = size(dumCs2_Pr4)
   displs = coord_of_rank * recvcounts
   call mpi_gatherv (dumCs2_Pr4(1),size(dumCs2_Pr4),MPI_RPREC,Cs2_Pr4_tot(1), &
                     recvcounts,displs,MPI_RPREC,rank_of_coord(0),comm,ierr)
$else
   avgCs_tot   = dumavgCs
   avgPr_t_tot = dumavgPr_t
   Cs2_Pr_tot  = dumCs2_Pr
   Cs2_Pr2_tot = dumCs2_Pr2
   Cs2_Pr4_tot = dumCs2_Pr4
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

end subroutine lagrange_scalar_Sdep
