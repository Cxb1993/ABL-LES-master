subroutine tecplot(plcount)
! Makes directory and writes data for post-processing in tecplot
! For movies and stuff

use types,only:rprec
use param
use sim_param,only:path,u,v,w,p,theta

implicit none
integer::jx,jy,jz
integer,intent(in)::plcount

$if ($MPI)
  integer::recvcounts(nproc)
  integer::displs(nproc)
$endif

real(kind=rprec),dimension(nx,ny,nz_tot-1)::u_tot,v_tot,w_tot,p_tot, &
                                        up_tot,vp_tot,wp_tot,pp_tot, &
                                        theta_tot,thetap_tot
real(kind=rprec),dimension(nx,ny,nz-1)::utemp,vtemp,wtemp,ptemp,thetatemp, &
                                   uptemp,vptemp,wptemp,pptemp,thetaptemp, &
                                   uprime,vprime,wprime,pprime,thetaprime
real(kind=rprec),dimension(nz-1)::uavg,vavg,wavg,pavg,thetaavg
real(kind=rprec)::frac
character(128)::folder,makedirectory,folderpath,filename1,filename2

frac = nx*ny ! To calculate average
do jz=1,nz-1
   uavg(jz) = 0._rprec;    vavg(jz) = 0._rprec
   wavg(jz) = 0._rprec;    pavg(jz) = 0._rprec
   thetaavg(jz) = 0._rprec
   do jx=1,nx
   do jy=1,ny
      uavg(jz)     = uavg(jz)     + u(jx,jy,jz)/frac
      vavg(jz)     = vavg(jz)     + v(jx,jy,jz)/frac
      wavg(jz)     = wavg(jz)     + w(jx,jy,jz)/frac
      pavg(jz)     = pavg(jz)     + p(jx,jy,jz)/frac
      thetaavg(jz) = thetaavg(jz) + theta(jx,jy,jz)/frac
   end do
   end do

   do jx=1,nx
      uprime(jx,:,jz)     = u(jx,:,jz)     - uavg(jz)
      vprime(jx,:,jz)     = v(jx,:,jz)     - vavg(jz)
      wprime(jx,:,jz)     = w(jx,:,jz)     - wavg(jz)
      pprime(jx,:,jz)     = p(jx,:,jz)     - pavg(jz)
      thetaprime(jx,:,jz) = theta(jx,:,jz) - thetaavg(jz)
   end do
end do

do jz=1,nz-1
   do jx=1,nx
      utemp(jx,:,jz)      = u(jx,:,jz)
      vtemp(jx,:,jz)      = v(jx,:,jz)
      wtemp(jx,:,jz)      = 0.5_rprec * ( w(jx,:,jz) + w(jx,:,jz+1) )
      ptemp(jx,:,jz)      = p(jx,:,jz)
      thetatemp(jx,:,jz)  = theta(jx,:,jz)

      uptemp(jx,:,jz)     = uprime(jx,:,jz)
      vptemp(jx,:,jz)     = vprime(jx,:,jz)
      wptemp(jx,:,jz)     = 0.5_rprec * ( wprime(jx,:,jz) + wprime(jx,:,jz+1) )
      pptemp(jx,:,jz)     = pprime(jx,:,jz)
      thetaptemp(jx,:,jz) = thetaprime(jx,:,jz)
   end do
end do

$if ($MPI)
   recvcounts = size(utemp)
   displs = coord_of_rank*recvcounts
   call mpi_gatherv (utemp(1,1,1),size(utemp),MPI_RPREC, &
                     u_tot(1,1,1),recvcounts,displs,     &
                     MPI_RPREC,rank_of_coord(0),comm,ierr)

   recvcounts = size(vtemp)
   displs = coord_of_rank*recvcounts
   call mpi_gatherv (vtemp(1,1,1),size(vtemp),MPI_RPREC, &
                     v_tot(1,1,1),recvcounts,displs,     &
                     MPI_RPREC,rank_of_coord(0),comm,ierr)

   recvcounts = size(wtemp)
   displs = coord_of_rank*recvcounts
   call mpi_gatherv (wtemp(1,1,1),size(wtemp),MPI_RPREC, &
                     w_tot(1,1,1),recvcounts,displs,     &
                     MPI_RPREC,rank_of_coord(0),comm,ierr)

   recvcounts = size(ptemp)
   displs = coord_of_rank*recvcounts
   call mpi_gatherv (ptemp(1,1,1),size(ptemp),MPI_RPREC, &
                     p_tot(1,1,1),recvcounts,displs,     &
                     MPI_RPREC,rank_of_coord(0),comm,ierr)

   recvcounts = size(thetatemp)
   displs = coord_of_rank*recvcounts
   call mpi_gatherv (thetatemp(1,1,1),size(thetatemp),MPI_RPREC, &
                     theta_tot(1,1,1),recvcounts,displs,         &
                     MPI_RPREC,rank_of_coord(0),comm,ierr)

   ! Fluctuating components
   recvcounts = size(uptemp)
   displs = coord_of_rank*recvcounts
   call mpi_gatherv (uptemp(1,1,1),size(uptemp),MPI_RPREC, &
                     up_tot(1,1,1),recvcounts,displs,      &
                     MPI_RPREC,rank_of_coord(0),comm,ierr)

   recvcounts = size(vptemp)
   displs = coord_of_rank*recvcounts
   call mpi_gatherv (vptemp(1,1,1),size(vptemp),MPI_RPREC, &
                     vp_tot(1,1,1),recvcounts,displs,      &
                     MPI_RPREC,rank_of_coord(0),comm,ierr)

   recvcounts = size(wptemp)
   displs = coord_of_rank*recvcounts
   call mpi_gatherv (wptemp(1,1,1),size(wptemp),MPI_RPREC, &
                     wp_tot(1,1,1),recvcounts,displs,      &
                     MPI_RPREC,rank_of_coord(0),comm,ierr)

   recvcounts = size(pptemp)
   displs = coord_of_rank*recvcounts
   call mpi_gatherv (pptemp(1,1,1),size(pptemp),MPI_RPREC, &
                     pp_tot(1,1,1),recvcounts,displs,      &
                     MPI_RPREC,rank_of_coord(0),comm,ierr)

   recvcounts = size(thetaptemp)
   displs = coord_of_rank*recvcounts
   call mpi_gatherv (thetaptemp(1,1,1),size(thetaptemp),MPI_RPREC, &
                     thetap_tot(1,1,1),recvcounts,displs,          &
                     MPI_RPREC,rank_of_coord(0),comm,ierr)
$else
   u_tot      = utemp
   v_tot      = vtemp
   w_tot      = wtemp
   p_tot      = ptemp
   theta_tot  = thetatemp

   up_tot     = uptemp
   vp_tot     = vptemp
   wp_tot     = wptemp
   pp_tot     = pptemp
   thetap_tot = thetaptemp
$endif

if((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
   folder = 'output/tecplot_data/'
   folderpath = path// trim(folder)
   makedirectory = 'mkdir ' // trim(folderpath)
   call system(makedirectory)

   write (filename1,'(a,i6.6,a)') trim(folderpath)//'stable_wt_01_inst',jt_total,'.dat'
   ! print*, trim(folderpath)
   ! print*, trim(makedirectory)
   ! print*, trim(filename1)
   open(unit=150,file=filename1,status="unknown",position="append")
   write(150,*)' VARIABLES = "X","Y","U","V","W","P","T"'
   write(150,*)' ZONE F=POINT, I= ',nx, ' J= ',0.8*nz_tot
   ! write(150,*)' ZONE F=POINT, I= ',nx, ' J= ',nz_tot-1+2
   ! Note that though VARIABLES are X & Y, data written is X & Z
   ! Similar is the case with ZONE F=POINT I & J
   ! +2 corresponds to the bottom-most and top-most layers

   ! ----------------------------------------------------------
   ! Top and bottom layers of the domain extrapolation
   ! ----------------------------------------------------------
   utemp(:,:,1)  = 0._rprec;  vtemp(:,:,1)  = 0._rprec
   wtemp(:,:,1)  = 0._rprec;  ptemp(:,:,1)  = 0._rprec
   thetatemp(:,:,1) = 1.5_rprec * theta_tot(:,:,1) - &
                      0.5_rprec * theta_tot(:,:,2)
   thetatemp(:,:,2) = 1.5_rprec * theta_tot(:,:,nz_tot-1) - &
                      0.5_rprec * theta_tot(:,:,nz_tot-2)

   uptemp(:,:,1) = 0._rprec;  vptemp(:,:,1) = 0._rprec
   wptemp(:,:,1) = 0._rprec;  pptemp(:,:,1) = 0._rprec
   thetaptemp(:,:,1) = 1.5_rprec * thetap_tot(:,:,1) - &
                       0.5_rprec * thetap_tot(:,:,2)
   thetaptemp(:,:,2) = 1.5_rprec * thetap_tot(:,:,nz_tot-1) - &
                       0.5_rprec * thetap_tot(:,:,nz_tot-2)
   ! ----------------------------------------------------------

   do jz = 0,0.8*nz_tot 
     ! 0.8*nz_tot corresponds to 80% of the domain height
     ! nz_tot-1+1 corresponds to start from 0th layer = wall
     ! do jy=1,ny
     jy = ny/2 ! Writing just the y = L_y/2 plane
     do jx=1,nx
        if (jz == 0) then
           write(150,100) (jx-1)*dx,0._rprec,            &
           utemp(jx,jy,1),vtemp(jx,jy,1),wtemp(jx,jy,1), &
           ptemp(jx,jy,1),thetatemp(jx,jy,1)
        else if (jz == nz_tot-1+1) then
           write(150,100) (jx-1)*dx,1._rprec,            &
           utemp(jx,jy,1),vtemp(jx,jy,1),wtemp(jx,jy,1), &
           ptemp(jx,jy,1),thetatemp(jx,jy,2)
           ! Note the difference in jz-indices of thetatemp
           ! at the top and bottom boundaries 
        else
           write(150,100) (jx-1)*dx,(2*jz-1)*dz/2._rprec,   &
           u_tot(jx,jy,jz),v_tot(jx,jy,jz),w_tot(jx,jy,jz), &
           p_tot(jx,jy,jz),theta_tot(jx,jy,jz)
        end if
     end do
     ! end do
   end do
   100  format (6e14.6)
   close(150)

   write (filename2,'(a,i6.6,a)') trim(folderpath)//'stable_wt_01_fluc',jt_total,'.dat'
   open(unit=151,file=filename2,status="unknown",position="append")
   write(151,*)' VARIABLES = "X","Y","Up","Vp","Wp","Pp","Tp"'
   write(151,*)' ZONE F=POINT, I= ',nx, ' J= ',0.8*nz_tot
   ! write(151,*)' ZONE F=POINT, I= ',nx, ' J= ',nz_tot-1+2
   ! Note that though VARIABLES are X & Y, data written is X & Z
   ! Similar is the case with ZONE F=POINT I & J

   do jz = 0,0.8*nz_tot 
     ! 0.8*nz_tot corresponds to 80% of the domain height
     ! nz_tot-1+1 corresponds to start from 0th layer = wall
     ! do jy=1,ny
     jy = ny/2 ! Writing just the y = L_y/2 plane
     do jx=1,nx
        if (jz == 0) then
           write(151,101) (jx-1)*dx,0._rprec,               &
           uptemp(jx,jy,1),vptemp(jx,jy,1),wptemp(jx,jy,1), &
           pptemp(jx,jy,1),thetaptemp(jx,jy,1)
        else if (jz == nz_tot-1+1) then
           write(151,101) (jx-1)*dx,1._rprec,               &
           uptemp(jx,jy,1),vptemp(jx,jy,1),wptemp(jx,jy,1), &
           pptemp(jx,jy,1),thetaptemp(jx,jy,2)
           ! Note the difference in jz-indices of thetaptemp 
           ! at the top and bottom boundaries 
        else
           write(151,101) (jx-1)*dx,(2*jz-1)*dz/2._rprec,      &
           up_tot(jx,jy,jz),vp_tot(jx,jy,jz),wp_tot(jx,jy,jz), &
           pp_tot(jx,jy,jz),thetap_tot(jx,jy,jz)
        end if
     end do
     ! end do
   end do
   101  format (7e14.5)
   close(151)

end if

end subroutine tecplot
