subroutine energy (ke)
use types,only:rprec
use param
use sim_param,only:u,v,w
use messages
! SKS
! USE, INTRINSIC :: IEEE_ARITHMETIC, IEEE_EXCEPTIONS
! SKS
implicit none

character (*), parameter :: sub_name = 'energy'

logical, parameter :: DEBUG = .true.
logical, parameter :: flush = .true.

integer::jx,jy,jz

logical :: nan_flag

real(kind=rprec)::KE,denom,temp_w
$if ($MPI)
  real (rprec) :: ke_global
$endif

nan_flag = .false.

ke=0._rprec
denom=2._rprec*(nx*ny*(nz-1))

do jz=1,nz-1
do jy=1,ny
do jx=1,nx
   temp_w=.5_rprec*(w(jx,jy,jz)+w(jx,jy,jz+1))
   ke=ke+(u(jx,jy,jz)**2+v(jx,jy,jz)**2+temp_w**2)/denom

   $if ($IFORT || $IFC)
   $endif

end do
end do
end do

if (nan_flag) call error (sub_name,'NaN found')

$if ($MPI)
  call mpi_reduce (ke, ke_global, 1, MPI_RPREC, MPI_SUM, 0, comm, ierr)
  if (rank == 0) then  !--note its rank here, not coord
    ke = ke_global/nproc
    write (13, *) (jt_total-1) * dt, ke
  end if
  !if (rank == 0) ke = ke_global/nproc  !--its rank here, not coord
$else
  write (13,*) (jt_total-1)*dt, ke
$endif

if ( flush ) then
  if ((.not. USE_MPI) .or. (USE_MPI .and. rank == 0)) then
    close (13)
    open (13,file=path//'output/check_ke.out', position='append')
  end if
end if

end subroutine energy
