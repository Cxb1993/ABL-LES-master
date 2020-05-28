subroutine initial()
use types,only:rprec
use param
use sim_param,only:path,u,v,w,RHSx,RHSy,RHSz,theta,q
use sgsmodule,only:Cs_opt2,Cs_opt2_avg,F_LM,F_MM,F_QN,F_NN, &
                   G_LM,G_MM,G_QN,G_NN,Pr_t
! use bottombc,only:psi_m ! added by VK
use scalars_module,only:RHS_T,sgs_t3,psi_m ! added by VK
use scalars_module2,only:ic_scal,ic_scal_gabls,ic_scal_GABLS_diurnal ! added by VK
! VK -label 122 assigned to vel_sc.out for reading input files in case of
! scalars
!!!!XXXXXXXXXX--------Added by Vijayant----XXXXXXX!!!!!
use immersedbc,only:fx,fy,fz,u_des,v_des,w_des,n_bldg,bldg_pts
use io,only:mean_u,mean_u2,mean_v,mean_v2,mean_w,mean_w2
implicit none

$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif

logical, parameter :: use_add_random = .false.

character (200) :: fname, temp
real(kind=rprec),dimension(ld,ny,nz)::crap
real(kind=rprec)::z
integer::i,jz

Cs_opt2_avg=0._rprec
fx=0._rprec;fy=0._rprec;fz=0._rprec
u_des=0._rprec;v_des=0._rprec;w_des=0._rprec
mean_u=0._rprec;mean_u2=0._rprec;mean_v=0._rprec;mean_v2=0._rprec
mean_w=0._rprec;mean_w2=0._rprec

!VK Modified so that when running scalars, the file is named vel_sc.out
!VK else called vel.out

if (S_FLAG) then
    fname = path // 'vel_sc.out'
else
    fname = path // 'vel.out'
end if

$if ($MPI)
  write (temp, '(".c",i0)') coord
  fname = trim (fname) // temp
$endif

open(11,file=fname,form='unformatted')

if(initu)then
  if(initsc) then
    print *,'Reading initial velocity and temperature from file'
    select case (model)
      case (1)
        read (11) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz),theta(:,:,1:nz), &
                  RHSx(:, :, 1:nz), RHSy(:, :, 1:nz),RHSz(:, :, 1:nz),         &
                  RHS_T(:,:,1:nz), sgs_t3(:,:,1), psi_m
      case (2:3)
        read (11) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),theta(:,:,1:nz),   &
                  RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),        &
                  RHS_T(:,:,1:nz), sgs_t3(:,:,1), psi_m, Cs_opt2
                 
      case (4)
        if (inilag) then  !--not sure if Cs_opt2 should be there, just quickie
          read (11) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),theta(:,:,1:nz), &
                    RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),      &
                    RHS_T(:,:,1:nz), sgs_t3(:,:,1), psi_m, Cs_opt2 
        else
          read (11) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),theta(:,:,1:nz),  &
                    RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),       &
                    RHS_T(:,:,1:nz),sgs_t3(:,:,1), psi_m, Cs_opt2, F_LM, F_MM,  &
                    crap, crap
        end if
      case (5)
        if (inilag) then  !--not sure if Cs_opt2 should be there, just quickie
          read (11) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),theta(:,:,1:nz), &
                    RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),      &
                    RHS_T(:,:,1:nz), sgs_t3(:,:,1), psi_m, Cs_opt2
        else
          read (11) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),theta(:,:,1:nz),  &
                    RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),       &
                    RHS_T(:,:,1:nz), sgs_t3(:,:,1), psi_m, Cs_opt2, F_LM, F_MM, &
                    F_QN, F_NN
        end if
      ! SKS
      case (6)
        read (11) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),theta(:,:,1:nz),   &
                  RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),        &
                  RHS_T(:,:,1:nz), sgs_t3(:,:,1), psi_m, Cs_opt2
      case (7)
        if (inilag) then  ! Not sure what all are needed when inilag = .true.
          read (11) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),theta(:,:,1:nz), &
                    RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),      &
                    RHS_T(:,:,1:nz), sgs_t3(:,:,1), psi_m, Cs_opt2
        else
          read (11) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),theta(:,:,1:nz), &
                    RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),      &
                    RHS_T(:,:,1:nz), sgs_t3(:,:,1), psi_m, Cs_opt2, F_LM, F_MM,&
                    F_QN, &
                    F_NN, &
                    G_LM, &
                    G_MM, &
                    G_QN, &
                    G_NN, &
                    Pr_t(:, :, 1:nz)
        end if
      ! SKS
      case default
        write (*,*) 'initial: invalid model number'
    end select
    
!TS INITIALIZE THE ZERO CONCENTRATION FIELD IF jt_total=0
    if(passive_scalar)then
       open(1,file=path//'run')
       read(1,*)i
       close(1)
       if(i.eq.0)then
          theta=0._rprec;RHS_T=0._rprec
       endif
    endif

open(unit=44,file=path//'output/init_profiles.dat',status="unknown",position="append")
    do jz=1,nz
         $if ($MPI)
            z = (coord*(nz-1) + jz - 0.5_rprec) * dz
         $else
            z=(real(jz)-0.5_rprec)*dz
         $endif
     write(6,7781) jz,z,(sum(u(:,:,jz))/float(nx*ny))*u_star,(sum(v(:,:,jz))/&
     float(nx*ny))*u_star,(sum(w(:,:,jz))/float(nx*ny))*u_star,&
     (sum(theta(:,:,jz))/float(nx*ny))*T_scale

     write(44,7781) jz,z,(sum(u(:,:,jz))/float(nx*ny))*u_star,(sum(v(:,:,jz))/&
     float(nx*ny))*u_star,(sum(w(:,:,jz))/float(nx*ny))*u_star,&
     (sum(theta(:,:,jz))/float(nx*ny))*T_scale

!    write(6,7781) jz,sum(u(:,:,jz))/float(nx*ny),sum(v(:,:,jz))/float(nx*ny)&
!         ,sum(w(:,:,jz))/float(nx*ny),sum(theta(:,:,jz))/float(nx*ny)
!    write(6,7782) jz,sum(Cs_opt2(:,:,jz))/float(nx*ny),sum(F_LM(:,:,jz))/float(nx*ny)&
!         ,sum(F_MM(:,:,jz))/float(nx*ny),sum(RHS_T(1:nx,1:ny,jz))/float(nx*ny)
    end do
 7781 format('jz, z, ubar, vbar, wbar,T_bar:',(1x,I3,1x,F9.4,1x,F9.4,1x,F9.4,1x,F9.4,1x,F9.4))
!7781 format('jz, ubar, vbar, wbar, Tbar:',(1x,I3,1x,F9.4,1x,F9.4,1x,F9.4,1x,F9.4))
!7782 format('jz, Cs_bar, F_LM_bar, F_MM_bar, RHS_Tbar:',(1x,I3,1x,F9.4,1x,F9.4,1x,F9.4,1x,F9.4))

!     do jz = 1, nz
!       print *,'jt,Cs_opt2,F_LM,F_MM',jt,sum(Cs_opt2(1:nx,1:ny,jz))/float(nx*ny)&
!      ,sum(F_LM(1:nx,1:ny,jz))/float(nx*ny),sum(F_MM(1:nx,1:ny,jz))/float(nx*ny)
!     end do

close(44)
  else

    print *,'Reading initial velocity field from file'

    select case (model)
      case (1)
        read (11) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz),       &
                  RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz)
      case (2:4)
        read(11) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),             &
                 RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),  &
                 Cs_opt2
      case (5)
        if (inilag) then  !--not sure if Cs_opt2 should be there, just quickie
          read (11) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),             &
                    RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),  &
                    Cs_opt2
        else
          read(11) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),             &
                   RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),  &
                   Cs_opt2, F_LM, F_MM, F_QN, F_NN
        end if
      ! SKS
      case (6)
        read(11) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),             &
                 RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),  &
                 Cs_opt2
      case (7)
        if (inilag) then  ! not sure what inilag does
          read (11) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),             &
                    RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),  &
                    Cs_opt2
        else
          read(11) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),             &
                   RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),  &
                   Cs_opt2, F_LM, F_MM, F_QN, F_NN
        end if
      ! SKS
      case default
        write (*, *) 'initial: invalid model number'
    end select

    do jz=1,nz
      write(6,7780) jz, sum (u(1:nx, :, jz)) / (nx * ny),  &
                        sum (v(1:nx, :, jz)) / (nx * ny),  &
                        sum (w(1:nx, :, jz)) / (nx * ny)
    end do
7780 format('jz, ubar, vbar, wbar:',(1x,I3,1x,F9.4,1x,F9.4,1x,F9.4))
  end if
else
  if (dns_bc) then
     print*, 'Creating initial velocity field with DNS BCs'
     call ic_dns()
  else
    print*, 'Creating initial fields'
    if (S_FLAG) then
       if (GABLS_test) then
         print *, 'Creating initial velocity & scalar fields for the GABLS test case'
         call ic_scal_GABLS()
       else if (GABLS_diurnal_test) then
         print *, 'Creating initial velocity & scalar fields for the GABLS diurnal case'
         call ic_scal_GABLS_diurnal()
       else
         print *, 'Creating initial velocity & scalar fields'
         call ic_scal()
       end if
    else
       ! SKS
       ! If S_FLAG == 0 then ic() is called, else the above ones are called
       print*, 'Creating initial velocity field'
       call ic()
       ! SKS
    end if
  end if
end if

! bldg stuff
if (use_bldg) then
   open(1,file=path//'bldg.dat')
   read(1,*) n_bldg
   allocate(bldg_pts(5,n_bldg))
   do i=1,n_bldg
      read(1,*) bldg_pts(1:5,i)
      if(bldg_pts(5,i).ge.nz)then
         write(*,*)"lz should be less than nz"
         stop
      end if
   end do
   close(1)
end if

$if ($MPI)
  !--synchronize the overlapping parts 0 <-> nz-1 and 1 <-> nz 
  call mpi_sendrecv (u(1, 1, nz-1), ld*ny, MPI_RPREC, up, 1,  &
                     u(1, 1, 0), ld*ny, MPI_RPREC, down, 1,   &
                     comm, status, ierr)
  call mpi_sendrecv (v(1, 1, nz-1), ld*ny, MPI_RPREC, up, 2,  &
                     v(1, 1, 0), ld*ny, MPI_RPREC, down, 2,   &
                     comm, status, ierr)
  call mpi_sendrecv (w(1, 1, nz-1), ld*ny, MPI_RPREC, up, 3,  &
                     w(1, 1, 0), ld*ny, MPI_RPREC, down, 3,   &
                     comm, status, ierr)
  call mpi_sendrecv (u(1, 1, 1), ld*ny, MPI_RPREC, down, 4,  &
                     u(1, 1, nz), ld*ny, MPI_RPREC, up, 4,   &
                     comm, status, ierr)
  call mpi_sendrecv (v(1, 1, 1), ld*ny, MPI_RPREC, down, 5,  &
                     v(1, 1, nz), ld*ny, MPI_RPREC, up, 5,   &
                     comm, status, ierr)
  call mpi_sendrecv (w(1, 1, 1), ld*ny, MPI_RPREC, down, 6,  &
                     w(1, 1, nz), ld*ny, MPI_RPREC, up, 6,   &
                     comm, status, ierr)
  call mpi_sendrecv (theta(1, 1, nz-1), ld*ny, MPI_RPREC, up, 7,  &
                     theta(1, 1, 0), ld*ny, MPI_RPREC, down, 7,   &
                     comm, status, ierr)
  call mpi_sendrecv (theta(1, 1, 1), ld*ny, MPI_RPREC, down, 8,  &
                     theta(1, 1, nz), ld*ny, MPI_RPREC, up, 8,   &
                     comm, status, ierr)
$endif

if (USE_MPI .and. coord == 0) then
  !--set 0-level velocities to BOGUS
  u(:, :, $lbz) = BOGUS
  v(:, :, $lbz) = BOGUS
  w(:, :, $lbz) = BOGUS
  theta(:, :, $lbz) = BOGUS
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine add_random ()
implicit none

real (rprec), parameter :: rms = 0.2_rprec
integer :: i, j, k
integer :: seed
real (rprec) :: noise
real (rprec) :: ran3

!---------------------------------------------------------------------
seed = -80

do k = 1, nz
  do j = 1, ny
    do i = 1, nx
      noise=rms/.289_rprec*(ran3(seed)-0.5_rprec)
      u(i, j, k) = u(i, j, k) + noise
      noise=rms/.289_rprec*(ran3(seed)-0.5_rprec)
      v(i, j, k) = v(i, j, k) + noise
      noise=rms/.289_rprec*(ran3(seed)-0.5_rprec)
      w(i, j, k) = w(i, j, k) + noise
    end do
  end do
end do

end subroutine add_random
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine initial
