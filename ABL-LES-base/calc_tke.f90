subroutine calc_tke()
! Subroutine written by Stimit Shah on 18th Sept 2010
! Calculates total TKE in all x-y planes
! comments added by Joseph Fogarty

use types,only:rprec
use param
use sim_param,only:path,u,v,w,theta
implicit none

integer::jx,jy,jz
real(kind=rprec)::u2,v2,w2,t2,arg1,arg2,arg3,arg4
real(kind=rprec),dimension(nz-1)::up2,vp2,wp2,thetap2 ! theta prime square
real(kind=rprec),dimension(nz_tot-1)::tot_up2,tot_vp2,tot_wp2, &
                                      tot_thetap2

$if ($MPI)
  integer::recvcounts(nproc)
  integer::displs(nproc)
  $define $lbz 0
$else
  $define $lbz 1
$endif

real(kind=rprec),dimension($lbz:nz)::ubar,vbar,wbar,Tbar

! calculate velocity averages
do jz = $lbz,nz
   ubar(jz) = sum(u(1:nx,1:ny,jz))/(nx*ny)
   vbar(jz) = sum(v(1:nx,1:ny,jz))/(nx*ny)
   wbar(jz) = sum(w(1:nx,1:ny,jz))/(nx*ny)
   Tbar(jz) = sum(theta(1:nx,1:ny,jz))/(nx*ny)
end do

! use velocity averages and instantaneous to find perturbation
do jz = 1,nz-1
   u2 = 0._rprec
   v2 = 0._rprec
   w2 = 0._rprec
   t2 = 0._rprec

   do jx = 1,nx
   do jy = 1,ny
      arg1 = (u(jx,jy,jz) - ubar(jz))
      arg2 = (v(jx,jy,jz) - vbar(jz))
      arg3 = (w(jx,jy,jz)   - wbar(jz) + &
              w(jx,jy,jz+1) - wbar(jz+1))/2._rprec
      arg4 = (theta(jx,jy,jz) - Tbar(jz))

      u2 = u2 + arg1*arg1
      v2 = v2 + arg2*arg2
      w2 = w2 + arg3*arg3
      t2 = t2 + arg4*arg4
   end do
   end do

   up2(jz)     = u2/(nx*ny)
   vp2(jz)     = v2/(nx*ny)
   wp2(jz)     = w2/(nx*ny)
   thetap2(jz) = t2/(nx*ny)
end do

$if ($MPI)
  recvcounts = size(up2)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (up2(1),size(up2),MPI_RPREC,   &
                    tot_up2(1),recvcounts,displs, &
                    MPI_RPREC,rank_of_coord(0),comm,ierr)

  recvcounts = size(vp2)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (vp2(1),size(vp2),MPI_RPREC,   &
                    tot_vp2(1),recvcounts,displs, &
                    MPI_RPREC,rank_of_coord(0),comm,ierr)

  recvcounts = size(wp2)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (wp2(1),size(wp2),MPI_RPREC,   &
                    tot_wp2(1),recvcounts,displs, &
                    MPI_RPREC,rank_of_coord(0),comm,ierr)

  recvcounts = size(thetap2)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (thetap2(1),size(thetap2),MPI_RPREC, &
                    tot_thetap2(1),recvcounts,displs,   &
                    MPI_RPREC,rank_of_coord(0),comm,ierr)
$else
  tot_up2     = up2
  tot_vp2     = vp2
  tot_wp2     = wp2
  tot_thetap2 = thetap2
$endif

if((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  do jz=1,nz_tot-1
    write(111,8000) (2*jz-1)*dz/2._rprec,tot_up2(jz), &
                tot_vp2(jz),tot_wp2(jz),tot_thetap2(jz)
  end do
  write(111,*)
end if
8000  format(1400(E14.5))

end subroutine calc_tke
