module topbc
use types,only:rprec
use param,only:nz,nz_tot,damping_method
implicit none
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif
!real(kind=rprec),dimension(nz)::sponge
real(kind=rprec),dimension($lbz:nz)::sponge

contains
subroutine setsponge()
use param
implicit none
real(kind=rprec)::factor

real(kind=rprec) :: z_d,cfrdmp,z_local,z
!real(kind=rprec),dimension($lbz:nz) :: sponge
integer :: k


$if ($MPI)
  integer, parameter :: nz_global = nz * nproc
  integer :: k_global
  real (rprec) :: sponge_top
$endif

if (damping_method==2) then
         z_d = 0.75_rprec*(L_z*real(nproc))
         cfrdmp = 3.9_rprec
         sponge=0._rprec
         do k=1,nz-1
         $if ($MPI)
            z_local = (coord*(nz-1) + k - 0.5_rprec) * dz*z_i
         $else
            z_local=(real(k)-0.5_rprec)*dz*z_i
         $endif

           if (z_local .ge. z_d .and. z_local .le. L_z*real(nproc)) then
               sponge(k)=cfrdmp*0.5_rprec*(1._rprec-cos(pi*(z_local-z_d)/(L_z-z_d)))
           else
               sponge(k)=0._rprec
           end if
         end do            
          
          if ((.not. USE_MPI) .or. (coord == nproc-1)) then
                sponge(nz)=sponge(nz-1)
          end if  
elseif (damping_method==1) then
! sets relaxation term to vertical momentum equation in top quarter
! of domain, relaxation timescale on top is 50s with a factor of 5 if we
! had Nz=40 for the layers 40...31, Nieuwstadt et al. 1991, turbulent shear 
! flows

$if ($MPI)
  !--the non-MPI recursive form in inconvenient, so instead replace with
  !  analytically evaluated version (should be this way anyway...)
  sponge=0._rprec
!  factor=9._rprec/(nz_global - 3*nz_global/4 + 1)
  factor=9._rprec/(nz_tot - 0.75*nz_tot + 1)
  sponge_top = z_i / (50._rprec * u_star)
  do k = 1, nz
    k_global = k + coord * nz
    if (k_global > 0.75*nz_global + 1) then
      sponge(k) = sponge_top * 5._rprec**((k_global-nz_global) * factor)
    end if
  end do
$else
  sponge=0._rprec
  factor=9._rprec/(nz-0.75*nz+1)
  sponge(nz)=z_i/(50._rprec*u_star)
  do k=nz-1,0.75*nz+1,-1
     sponge(k)=sponge(k+1)/5._rprec**factor
  end do
$endif

end if

open(unit=6,file=path//'output/sponge.dat',status="unknown",position="append")

do k=1,nz-1
     $if ($MPI)
       z = (coord*(nz-1) + k - 0.5_rprec) * dz
     $else
       z = (k - 0.5_rprec) * dz
     $endif
     write(6,7781) k,z,sponge(k),z_d,L_z
end do
close(6)
7781 format('jz, z, sponge, z_d, L_x:',(1x,I3,1x,F11.4,1x,F11.4,1x,F11.4,1x,F11.4))


end subroutine setsponge
end module topbc
