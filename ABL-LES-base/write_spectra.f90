subroutine write_spectra()
! Writes spectra of u,v,w,and T calculated in x-direction

use sim_param,only:path,avg_spectra_uvwT
use fft
use param
use types,only:rprec
implicit none
   
real(kind=rprec),dimension(4,nx/2,nz_tot-1)::avg_spectra_uvwT_tot
integer::k,jx,jy,jz,z

$if ($MPI)
  $define $lbz 0
  integer::recvcounts(nproc)
  integer::displs(nproc)
$else    
  $define $lbz 1
$endif

avg_spectra_uvwT(:,:,:) = avg_spectra_uvwT(:,:,:)/(nsteps-spectraCALC)

$if ($MPI)
   recvcounts = size(avg_spectra_uvwT)
   displs = coord_of_rank*recvcounts
   call mpi_gatherv (avg_spectra_uvwT(1,1,1),size(avg_spectra_uvwT),MPI_RPREC, &
                     avg_spectra_uvwT_tot(1,1,1),recvcounts,displs,        &
                     MPI_RPREC,rank_of_coord(0),comm,ierr)
$else
   avg_spectra_uvwT_tot = avg_spectra_uvwT
$endif

if((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
   do jz=1,nz_tot-1
      do jx=1,nx/2
         write(84,8000) avg_spectra_uvwT_tot(1,jx,jz)
         write(85,8000) avg_spectra_uvwT_tot(2,jx,jz)
         write(86,8000) avg_spectra_uvwT_tot(3,jx,jz)
         write(87,8000) avg_spectra_uvwT_tot(4,jx,jz)
      end do
      write(84,*)
      write(85,*)
      write(86,*)
      write(87,*)
   end do
end if
8000  format(1400(E14.5))

if((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
   if (jt == spectraCALC) then
      do jz=1,nz_tot-1
         z = (jz-0.5_rprec)*dz*z_i
         write(88,8002) (real(kx(jx,1)/z_i*z),jx=1,nx/2)
      end do
   end if
end if
8002  format(1400(E14.5))

end subroutine write_spectra
