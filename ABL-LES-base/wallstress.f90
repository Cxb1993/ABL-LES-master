! For use with staggered grid LES
! JDA, 23 Jan 96
! zo is nondimensionalized, zo1 not!
!--provides txz, tyz, dudz, dvdz at jz=1
subroutine wallstress ()
use types,only:rprec
use param,only:jt,nsteps,dz,ld,lh,nx,ny,nz,vonk,lbc_mom,WALL_HOMOG_FLAG
use sim_param,only:u,v,dudz,dvdz,txz,tyz
use bottombc,only:zo,psi_m,phi_m,ustar_avg
use test_filtermodule
implicit none
integer::jx,jy
real(kind=rprec),dimension(nx,ny)::u_avg,u_aver,denom
real(kind=rprec),dimension(ld,ny)::u1,v1
real(kind=rprec)::const

select case (lbc_mom)

  case ('wall')

    u1=u(:,:,1)
    v1=v(:,:,1)
!   if ((patch_flag .eq. 1) .AND. (num_patch .eq. 1)) then
    if (WALL_HOMOG_FLAG) then
!     calculate u_star in the average sense !!
!     u1=sum(u1(1:nx,1:ny))/float(nx*ny);v1=sum(v1(1:nx,1:ny))/float(nx*ny)
      denom=log(0.5_rprec*dz/zo)-sum(psi_m)/float(nx*ny)
      u_avg=sum(sqrt(u1(1:nx,1:ny)**2+v1(1:nx,1:ny)**2))/float(nx*ny)
      if (jt .eq. nsteps) print *,'USED WALL HOMOG conds in wallstress'
    else
      call test_filter(u1,G_test)
      call test_filter(v1,G_test)
      denom=log(0.5_rprec*dz/zo)-psi_m
      u_avg=sqrt(u1(1:nx,1:ny)**2+v1(1:nx,1:ny)**2)
    end if
!   ustar=u_avg*vonk/denom
    ustar_avg=u_avg*vonk/denom

    do jy=1,ny
    do jx=1,nx
       ! const=-(ustar(jx,jy)**2) /u_avg(jx,jy)
       const=-(ustar_avg(jx,jy)**2) /u_avg(jx,jy)
       txz(jx,jy,1)=const *u1(jx,jy)
       tyz(jx,jy,1)=const *v1(jx,jy)
       ! this is as in Moeng 84
       dudz(jx,jy,1)=ustar_avg(jx,jy)/(0.5_rprec*dz*vonK)*u(jx,jy,1)/u_avg(jx,jy)&
           *phi_m(jx,jy)
       dvdz(jx,jy,1)=ustar_avg(jx,jy)/(0.5_rprec*dz*vonK)*v(jx,jy,1)/u_avg(jx,jy)&
           *phi_m(jx,jy)
       dudz(jx,jy,1)=merge(0._rprec,dudz(jx,jy,1),u(jx,jy,1).eq.0._rprec)
       dvdz(jx,jy,1)=merge(0._rprec,dvdz(jx,jy,1),v(jx,jy,1).eq.0._rprec)
    end do
    end do

  case ('stress free')

    txz(:, :, 1) = 0._rprec
    tyz(:, :, 1) = 0._rprec
    dudz(:, :, 1) = 0._rprec
    dvdz(:, :, 1) = 0._rprec

  case default

    write (*, *) 'invalid lbc_mom'
    stop

end select

end subroutine wallstress
