subroutine spectra()
! Calls function to calculate spectra and the writes it

use sim_param,only:path,u,v,w,theta,avg_spectra_uvwT
use param
use fft
use types,only:rprec
implicit none

real(kind=rprec),dimension(4,nx/2,nz-1)::spectra_uvwT
real(kind=rprec),dimension(4,nx/2,nz_tot-1)::spectra_uvwT_tot
integer::k,jx,jy,jz,z

$if ($MPI)
  $define $lbz 0
  integer::recvcounts(nproc)
  integer::displs(nproc)
$else
  $define $lbz 1
$endif

do jz=1,nz-1
   call calc_spectra(u(:,:,jz),spectra_uvwT(1,:,jz))
   call calc_spectra(v(:,:,jz),spectra_uvwT(2,:,jz))
   call calc_spectra(w(:,:,jz),spectra_uvwT(3,:,jz))
   call calc_spectra(theta(:,:,jz),spectra_uvwT(4,:,jz))
enddo

do jy=1,4
  do jx=1,nx/2
    do jz=1,nz-1
       avg_spectra_uvwT(jy,jx,jz) = avg_spectra_uvwT(jy,jx,jz) &
                                  + spectra_uvwT(jy,jx,jz)
    end do
 end do
end do

end subroutine spectra
