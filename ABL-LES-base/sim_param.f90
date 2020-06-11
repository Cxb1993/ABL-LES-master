module sim_param
use types,only:rprec
use param,only:lh,ld,nx,ny,nz,path
implicit none

save
public

$if ($MPI)
  !--this dimensioning adds a ghost layer for finite differences
  !--its simpler to have all arrays dimensioned the same, even though
  !--some components do not need ghost layer
  $define $lbz 0
$else
  $define $lbz 1
$endif

real(rprec),dimension(ld,ny,$lbz:nz)::u,v,w
real(rprec),dimension(ld,ny,$lbz:nz)::dudx,dudy,dudz, &
                                      dvdx,dvdy,dvdz, &
                                      dwdx,dwdy,dwdz, &
                                      RHSx,RHSy,RHSz, &
                                      RHSx_f,RHSy_f,RHSz_f

real(rprec),dimension(ld,ny,nz)::dpdx=0._rprec, &
                                 dpdy=0._rprec, &
                                 dpdz=0._rprec

real(rprec),dimension(ld,ny,$lbz:nz)::txx,txy,tyy
real(rprec),dimension(ld,ny,$lbz:nz)::txz,tyz,tzz
real(kind=rprec),dimension(ld,ny,0:nz)::p
real(rprec),dimension(ld,ny,$lbz:nz)::divtx,divty,divtz
real(kind=rprec),dimension(ld,ny,$lbz:nz)::theta,q
!MM
real ::u_old
!MM


!MM
real(kind=rprec),dimension(ld,ny,nz)::L11t,L22t,L33t,Q11t,Q22t,Q33t
! SKS
real(kind=rprec),dimension(4,nx/2,nz-1)::avg_spectra_uvwT=0._rprec
! SKS

end module sim_param
