! Log profile that is modified to flatten at z=z_i
subroutine ic()
use types,only:rprec
use param
use sim_param,only:u,v,w
use bottombc
implicit none

logical, parameter :: DEBUG = .false.

integer::jx,jy,jz,seed
integer :: jz_abs

real(kind=rprec),dimension(nz)::ubar,vbar
real(kind=rprec)::rms, noise, arg, arg2
real(kind=rprec)::z,w_star,T_star,q_star,ran3
real(kind=rprec)::gam

if ((inflow) .and. (.not. read_inflow_file)) then  !--no turbulence
  u = face_avg 
  v = 0._rprec
  w = 0._rprec
else 
! MM I changed here from lower one to upper one 
   w_star=(9.81_rprec/T_init*wt_s*z_i)**(1._rprec/3._rprec)
  !  w_star=3._rprec
  !      T_star=wt_s/w_star
  !      q_star=T_star

  !$if ($MPI)
  !  seed = -80 - coord
  !$else
  !  seed=-80
  !$endif

  print *,'Modified Log Profile for IC'
  do jz=1,nz

    $if ($MPI)
      z = (coord*(nz-1) + jz - 0.5_rprec) * dz
    $else
      z = (jz-.5_rprec)*dz
    $endif
    ! IC in equilibrium with rough surface (rough dominates in effective zo)
    arg2=z/(sum(zo)/float(nx*ny))
    arg=(1._rprec/vonk)*log(arg2)!-1./(2.*vonk*z_i*z_i)*z*z

        if (coriolis_forcing) then
!       ubar(jz)=arg/30

! MM Initializgin with Ekman Steady Solution:
	gam=sqrt(u_star*coriol/(10*z_i))
	!MM for now lets ug=1
        !ug=1._rprec
	ubar(jz)=ug_dim*(1-exp(-gam*z*z_i)*cos(gam*z*z_i))/u_star
        vbar(jz)=ug_dim*exp(-gam*z*z_i)*sin(gam*z*z_i)/u_star
! MM

! Note that ug and vg have already been non-dimensionalized in param.f90
!       ubar(jz)=arg/30._rprec
        else
        ubar(jz)=arg
        vbar(jz)=0._rprec
        end if

     if ((coriolis_forcing).and.(z.gt.(.5_rprec*z_i))) ubar(jz)=ug_dim/u_star

  end do

  rms = 3._rprec
  seed = -80
  do jz=1,nz
  $if ($MPI)
    jz_abs = coord * (nz-1) + jz
    z = (coord * (nz-1) + jz - 0.5_rprec) * dz * z_i
    seed = -80 - jz_abs  !--trying to make consistent init for MPI
  $else
    jz_abs = jz
    z = (jz-.5_rprec) * dz * z_i
  $endif
  ! SKS
  ! seed = -80 - jz_abs  !--trying to make consistent init for MPI
  ! SKS
  do jy=1,ny
  do jx=1,nx
    !...Ran3 returns uniform RV between 0 and 1. (sigma_rv=0.289)
    !...Taking std dev of vel as 1 at all heights
    if (z.le.z_i) then
       noise=rms/.289_rprec*(ran3(seed)-.5_rprec)
       u(jx,jy,jz)=noise*(1._rprec-z/z_i)*w_star/u_star+ubar(jz)
       noise=rms/.289_rprec*(ran3(seed)-0.5_rprec)
       v(jx,jy,jz)=noise*(1._rprec-z/z_i)*w_star/u_star+vbar(jz) !noise
       noise=rms/.289_rprec*(ran3(seed)-.5_rprec)
       w(jx,jy,jz)=noise*(1._rprec-z/z_i)*w_star/u_star
    else
       noise=rms/.289_rprec*(ran3(seed)-.5_rprec)
       u(jx,jy,jz)=noise*w_star/u_star*.01_rprec+ubar(jz)
       noise=rms/.289_rprec*(ran3(seed)-0.5_rprec)
       v(jx,jy,jz)=noise*w_star/u_star*.01_rprec+vbar(jz)
       noise=rms/.289_rprec*(ran3(seed)-0.5_rprec)
       w(jx,jy,jz)=noise*w_star/u_star*.01_rprec
    end if
  end do
  end do
  end do

  !...BC for W
  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    w(1:nx, 1:ny, 1) = 0._rprec
  end if
  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
    w(1:nx, 1:ny, nz) = 0._rprec
  endif

  !...BC for U, V
  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
    u(1:nx, 1:ny, nz) = u(1:nx, 1:ny, nz-1)
    v(1:nx, 1:ny, nz) = v(1:nx, 1:ny, nz-1)
  end if

end if

!VK Display the mean vertical profiles of the initialized variables on the
!screen
do jz=1,nz
     $if ($MPI)
       z = (coord*(nz-1) + jz - 0.5_rprec) * dz
     $else
       z = (jz - 0.5_rprec) * dz
     $endif
     write(6,7780) jz,z,sum(u(1:nx,:,jz))/float(nx*ny),sum(v(1:nx,:,jz))/&
     float(nx*ny),sum(w(1:nx,:,jz))/float(nx*ny)
end do
7780 format('jz,z,ubar,vbar,wbar:',(1x,I3,4(1x,F9.4)))

end subroutine ic
