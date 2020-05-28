module io
use types,only:rprec
use param, only : ld, nx, ny, nz, write_inflow_file, path,  &
                  use_avgslice, USE_MPI, coord, rank, nproc,      &
                  average_dim_num, nz_tot,jt_total,p_count,dt,z_i,u_star

!use sim_param,only:L11t,L22t,L33t,Q11t,Q22t,Q33t
implicit none
private
public :: jt_total
public :: openfiles,output_loop,output_final,            &
          mean_u,mean_u2,mean_v,mean_v2,mean_w,mean_w2,  &
          inflow_read, inflow_write, avg_stats
public :: average_dim_select_flag, dim1_size, dim2_size,        &
          dim1_global, dim2_global, collocate_MPI_averages,     &
          compute_avg_var 

integer,parameter::num_hour_out=1
integer,parameter::base=50000,nwrite=base
! SKS
! Shouldnt be hard coded..base corresponds to the 
! time interval after which you want to write file
! So can be a factor of nsteps.
! SKS

logical,parameter:: io_spec=.false.,output_fields_3d_flag=.false.
integer,parameter::spec_write_freqz=600, fields_3d_write_freqz=p_count*6
integer,parameter::spec_write_start=1,spec_write_end=24*base
!! --------------------------------------------------------------------
!! The following block defines parameters for instantaneous slice output
!! slice_inst sets the value of the y node for the chosen x-z inst slice
!! inst_slice_freqz controls the output frequency
!! The 5 slices outputted every inst_slice_freqz (i.e. u,v,w,T,Cs in this order) ...
!! ... are saved in the 3rd dimension of inst_slice_array for inst_array_lim times ...
!! ... following which this array is outputted as a binary file and the process 
!! ... starts over
!! Therefore, each file will contain inst_array_lim x-z slices of 5 variables
!! This information is important for post-processing of these binary files
!! --------------------------------------------------------------------

logical,parameter:: inst_slice_flag=.false.
integer,parameter:: num_vars=4 ! works currently only for u,v,w,T due to the size difference in Cs
integer,parameter:: slice_inst=(nz_tot-1)/2, inst_slice_freqz=5, inst_array_lim=200

logical, parameter :: cumulative_time = .true.
character (*), parameter :: fcumulative_time = path // 'total_time.dat'

integer, parameter :: n_avg_stats = 5000000 !--interval for updates in avg_stats
character (*), parameter :: end_hdr_avg = '# end header'

!! --------------------------------------------------------------------
!! The following block defines parameters for use in avgslice and scalar_slice
!! --------------------------------------------------------------------
integer,parameter :: average_dim_select_flag=1-(average_dim_num/2) 
! The variable average_dim_select_flag generates the following values based
! on the value of average_dim_num in param.f90 :- 
! a) average_dim_num = 2 : 
! 	average_dim_select_flag = 0 ==> average the 3d array over x and y and output the z profile
! b) average_dim_num = 1 : 
! 	average_dim_select_flag = 1 ==> average the 3d array over y and output the (x,z) profile
integer, parameter :: dim1_size=average_dim_select_flag*(nx-nz+1)+nz-1
integer, parameter :: dim2_size=average_dim_select_flag*(nz-2)+1
integer, parameter :: dim1_global=average_dim_select_flag*(nx-nz_tot+1)+nz_tot-1
integer, parameter :: dim2_global=average_dim_select_flag*(nz_tot-2)+1
!! --------------------------------------------------------------------
!! --------------------------------------------------------------------

!!!!  time_spec>0 output time series spectrum (need additional calcu.)
integer,parameter::time_spec=0
integer::n_obs, jt_total_init
integer,allocatable::obs_pt(:,:)

!!!!  io_mean=.true. output small domain time-averaged velocity
logical,parameter::io_mean=.true.
integer,parameter::jx_pls=1,jx_ple=nx,width=0
integer,parameter::jy_pls=ny/2-width,jy_ple=ny/2+width
real(kind=rprec),dimension(jx_pls:jx_ple,jy_pls:jy_ple,nz)::&
     mean_u,mean_v,mean_w,mean_u2,mean_v2,mean_w2

!!!!  io_lambda2
logical,parameter::io_lambda2=.false.
real(kind=rprec),dimension(nx,ny,nz)::lam2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! file number 1-10 are used for temporary use
! 11-19 are basic output
! 20-40 are avgslice
! use >50 for debugging
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine openfiles()
use sim_param,only:path
implicit none

!--to hold file names
character (64) :: temp
character (64) :: fCS1plan, fCS2plan, fCS4plan, fVISCplan,  &
                  fDISSplan, fCS1Vplan, fCS2Vplan, fCS4Vplan

integer::i

logical :: exst

if (cumulative_time) then

  inquire (file=fcumulative_time, exist=exst)
  if (exst) then
    open (1, file=fcumulative_time)
    read (1, *) jt_total
    jt_total_init=jt_total 
    close (1)
  else  !--assume this is the first run on cumulative time
    write (*, *) 'file ', fcumulative_time, ' not found'
    write (*, *) 'assuming jt_total = 0'
    jt_total = 0
    jt_total_init=jt_total 
  end if

end if

if ((.not. USE_MPI) .or. (USE_MPI .and. rank == 0)) then
  open(13,file=path//'output/check_ke.out',status="unknown",position="append")
  ! SKS
  open(14,file=path//'output/Cs2byPr_t.out',status="unknown",position="append")
  open(84,file=path//'output/spec_u.out',status="unknown",position="append")
  open(85,file=path//'output/spec_v.out',status="unknown",position="append")
  open(86,file=path//'output/spec_w.out',status="unknown",position="append")
  open(87,file=path//'output/spec_T.out',status="unknown",position="append")
  open(88,file=path//'output/spec_freq.out',status="unknown",position="append")
  open(111,file=path//'output/tkeVsz.dat',status="unknown",position="append")
  ! SKS
end if

if(time_spec.gt.0)then
  open(15,file=path//'output/velspec.out',form='unformatted',position='append')
  if(jt_total.eq.0)rewind(15)
endif

if(io_mean)then
  open(51,file=path//'output/mean_u.out',form='unformatted',position='append')
  if(jt_total.eq.0)then
    rewind(51)
    write(51)jx_pls,jx_ple,jy_pls,jy_ple
  endif
endif

fCS1plan = path // 'output/CS1plan.out'
fCS2plan = path // 'output/CS2plan.out'
fCS4plan = path // 'output/CS4plan.out'
fVISCplan = path // 'output/VISCplan.out'
fDISSplan = path // 'output/DISSplan.out'
fCS1Vplan = path // 'output/CS1Vplan.out'
fCS2Vplan = path // 'output/CS2Vplan.out'
fCS4Vplan = path // 'output/CS4Vplan.out'

$if ($MPI)
  !--append coordinate identifiers
  write (temp, '(".c",i0)') coord
  fCS1plan = trim (fCS1plan) // temp
  fCS2plan = trim (fCS2plan) // temp
  fCS4plan = trim (fCS4plan) // temp
  fVISCplan = trim (fVISCplan) // temp
  fDISSplan = trim (fDISSplan) // temp
  fCS1Vplan = trim (fCS1Vplan) // temp
  fCS2Vplan = trim (fCS2Vplan) // temp
  fCS4Vplan = trim (fCS4Vplan) // temp
$endif

if(time_spec.gt.0)then
open(1,file=path//'obs.pt')
read(1,*)n_obs
allocate(obs_pt(1:2,n_obs))
do i=1,n_obs
read(1,*)obs_pt(1:2,i)
enddo
close(1)
endif

end subroutine openfiles
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_loop
use param,only:output,dt,c_count,S_FLAG,SCAL_init,jt,jan_diurnal_run
use sim_param,only:path,u,v,w,dudz,dudx,p,&
     RHSx,RHSy,RHSz,theta, txx, txy, txz, tyy, tyz, tzz
use sgsmodule,only:Cs_opt2
use scalars_module,only:sgs_t3
use scalars_module2,only:scalar_slice,budget_TKE_scalar
implicit none
real(kind=rprec),dimension(nz)::u_ndim
character(len=20)::req

! SKS
character (100) :: fname, temp  ! With 64 the code was giving an error !
! character (64) :: fname, temp ! sort of segmentation fault i guess.
! SKS

integer::jx,jy,jz
integer:: fields_3d_write_freqz_temp

$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif

real(kind=rprec),dimension(ld,$lbz:nz,num_vars*inst_array_lim),save::inst_slice_array
integer,save:: inst_slice_counter

!jt_total=jt_total+1  !--moved into main
call calculate_mean

if ((use_avgslice) .and. (mod(jt,c_count)==0)) then
       if ((S_FLAG) .and. (jt.GE.SCAL_init)) then ! Output scalar variables
         !call avgslice()
         call MM_budget_slice()
         call scalar_slice() ! Uses file unit numbers (35-75)
         ! SKS
         call budget_TKE_scalar()
     !    call MM_XYZ_Out()
         ! SKS
       elseif (.not. S_FLAG) then
         !call avgslice()
         call MM_budget_slice()
       end if
end if

if (output) then
  if ((mod(jt_total,base)==0).and.(jt_total.ge.1)) then
    if (S_FLAG) then
       write (fname, '(a,i6.6,a)') path // 'output/vel_sc', jt_total, '.out'
    else
       write (fname, '(a,i6.6,a)') path // 'output/vel', jt_total, '.out'
    end if
    $if ($MPI)
       write (temp, '(".c",i0)') coord
       fname = trim (fname) // temp
    $endif

    open(1,file=fname,form='unformatted')
    call checkpoint (1)
    close(1)
    if (io_mean) call io_mean_out
  end if
end if

 if (S_FLAG) then ! If Block 1
      if ((inst_slice_flag) .AND. mod(jt_total, inst_slice_freqz)==0) then !If Block 2
        if (jt .eq. inst_slice_freqz) inst_slice_counter=1
        if (jt .eq. inst_slice_freqz) print *,'inst_slice_counter = ',inst_slice_counter
            
            if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
                print *,'inst_slice_counter = ',inst_slice_counter
            end if

       inst_slice_array(:,:,(inst_slice_counter-1)*num_vars+1) = u(:,slice_inst,:)
       inst_slice_array(:,:,(inst_slice_counter-1)*num_vars+2) = v(:,slice_inst,:)
       inst_slice_array(:,:,(inst_slice_counter-1)*num_vars+3) = w(:,slice_inst,:)
       inst_slice_array(:,:,(inst_slice_counter-1)*num_vars+4) = theta(:,slice_inst,:)

         if (mod(inst_slice_counter,inst_array_lim) .eq. 0) then !If Block 3 begins
            if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then !If Block 4
                print *,'INSIDE:inst_slice_counter = ',inst_slice_counter
            end if !If Block 4 ends
          write(fname,'(A,i6.6,A)')path//'output/fields_3d/inst_slices_uvwT_till_',jt_total,'.bin'
          $if ($MPI)
            write (temp, '(".c",i0)') coord
            fname = trim (fname) // temp
          $endif
           open(1,file=fname,form='unformatted')
           write(1) real(inst_slice_array)
           close(1)
           inst_slice_array=0._rprec; inst_slice_counter=0;
         end if ! If Block 3 ends

       inst_slice_counter = inst_slice_counter+1 !increment the slice counter by 1
            if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then !If Block 4
                print *,'inst_slice_counter = ',inst_slice_counter
            end if !If Block 4 ends
     end if ! If Block 2 ends
       
       fields_3d_write_freqz_temp=50

    if ((output_fields_3d_flag) .and. mod(jt_total,fields_3d_write_freqz_temp)==0) then !If Block 5 begins

    write(fname,'(A,i6.6,A)')path//'output/fields_3d/o_uvwT_',jt_total,'.bin'
    $if ($MPI)
      write (temp, '(".c",i0)') coord
      fname = trim (fname) // temp
    $endif
    open(1,file=fname,form='unformatted')
    write(1) real(u),real(v),real(w),real(theta); close(1)
    
    end if ! If Block 5 ends
 end if ! If Block 1 ends

  if ((io_spec) .and. (jt_total .gt. spec_write_start .and. jt_total .le. spec_write_end)) then 
   if (mod(jt_total,spec_write_freqz)==0) call post_spec(jt_total)
  end if

  if (time_spec.gt.0) call timeseries_spec

end subroutine output_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine compute_avg_var(avg_var,raw_var,output_dim)
use param,only:nx,ny,nz,c_count,p_count
integer :: output_dim
real(kind=rprec),dimension(dim1_size,dim2_size)::avg_var
real(kind=rprec),dimension(1:nx,1:ny,1:nz-1):: raw_var
real(kind=rprec):: fr
     
fr=(1._rprec/real(p_count,kind=rprec))*real(c_count,kind=rprec)
select case (output_dim)
  case(1) !average over y and t
      avg_var=avg_var+fr*sum(raw_var(1:nx,1:ny,1:nz-1),dim=2)/ny
  case(2) ! average over x,y and t
      avg_var(:,1)=avg_var(:,1)+fr*sum(sum(raw_var(1:nx,1:ny,1:nz-1),dim=1),dim=2)/(nx*ny)
end select
end subroutine compute_avg_var

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--The following subroutine does the collocation of the MPI arrays for
! averaging in avgslice and scalar_slice (in scalars_module2.f90)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine collocate_MPI_averages(avg_var_proc,avg_var_tot_domain,file_ind)
use param
$if ($MPI)
  integer :: recvcounts(nproc)
  integer :: displs(nproc)
$endif
integer :: ind1,ind2,file_ind
real(kind=rprec),dimension(dim1_size,dim2_size)::avg_var_proc
real(kind=rprec),dimension(dim1_global,dim2_global)::avg_var_tot_domain

  avg_var_tot_domain=0._rprec
$if ($MPI)
  recvcounts = size (avg_var_proc)
  displs = coord_of_rank * recvcounts 
  call mpi_gatherv (avg_var_proc(1,1), size (avg_var_proc), MPI_RPREC,                &
                    avg_var_tot_domain(1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)
$else
  avg_var_tot_domain=avg_var_proc
$endif

  if((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then  
      if (average_dim_num .eq. 1) then
         do ind2=1,nz_tot-1
          write(file_ind,5168) jt*dt,(avg_var_tot_domain(ind1,ind2),ind1=1,nx)
         end do
      else if (average_dim_num .eq. 2) then
         write(file_ind,5168) jt*dt,(avg_var_tot_domain(ind1,1),ind1=1,nz_tot-1)
      end if
         close(file_ind)
  end if

5168     format(1400(E14.5))
end subroutine collocate_MPI_averages

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--The following subroutine does the collocation of the MPI arrays for
! averaging in avgslice and scalar_slice (in scalars_module2.f90)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine collocate_MPI_averages_N(avg_var_proc,avg_var_tot_domain,file_ind,filename_str)
!subroutine collocate_MPI_averages(avg_var_proc,avg_var_tot_domain,file_ind)
use param
$if ($MPI)
  integer :: recvcounts(nproc)
  integer :: displs(nproc)
$endif
integer :: ind1,ind2,file_ind
character (*),intent(in) :: filename_str
character (len=256) :: local_filename
real(kind=rprec),dimension(dim1_size,dim2_size)::avg_var_proc
real(kind=rprec),dimension(dim1_global,dim2_global)::avg_var_tot_domain

local_filename=path//'output/aver_'//trim(filename_str)//'.out'

avg_var_tot_domain=0._rprec
$if ($MPI)
  recvcounts = size (avg_var_proc)
  displs = coord_of_rank * recvcounts 
  call mpi_gatherv (avg_var_proc(1,1), size (avg_var_proc), MPI_RPREC,       &
                    avg_var_tot_domain(1, 1), recvcounts, displs, MPI_RPREC, &
                    rank_of_coord(0), comm, ierr)
$else
  avg_var_tot_domain=avg_var_proc
$endif

  if((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then  
        open(file_ind,file=trim(local_filename),status="unknown",position="append")
        if (average_dim_num .eq. 1) then
           do ind2=1,nz_tot-1
            write(file_ind,5168) jt*dt,(avg_var_tot_domain(ind1,ind2),ind1=1,nx)
           end do
        else if (average_dim_num .eq. 2) then
           write(file_ind,5168) jt*dt,(avg_var_tot_domain(ind1,1),ind1=1,nz_tot-1)
        end if
        close(file_ind)
  end if
5168     format(1400(E14.5))
end subroutine collocate_MPI_averages_N

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!subroutine avgslice
!use sim_param,only:path,u,v,w,dudz,dvdz,txx,txz,tyy,tyz,tzz,p
!use param,only:dz,p_count,c_count,jt
!use sgsmodule,only:Cs_opt2,Cs_Ssim,Beta_avg,Betaclip_avg
!implicit none
!integer::i,j,k
!real(kind=rprec),dimension(nx,nz-1),save::ap,au,av,aw,p2,u2,v2,w2,auw,avw,acs
!real(kind=rprec),dimension(nx,nz-1),save::adudz,advdz,aCs_Ssim,abeta_sgs,abetaclip_sgs
!real(kind=rprec),dimension(nx,nz-1),save::atxx,atxz,atyy,atyz,atzz
!real(kind=rprec),dimension(nx,nz-1),save::u3,v3,w3
!real(kind=rprec)::tu1,tv1,tw1,ttxx,ttxz,ttyy,ttyz,ttzz,tdudz,tdvdz,&
!     tu2,tv2,tw2,tp1,tp2,tuw,tvw,tCs,fr,arg1,arg2,tu3,tv3,tw3
!real(kind=rprec)::tCs_Ssim
!real(kind=rprec),dimension(:,:),allocatable::avg_out
!
!fr=(1._rprec/real(p_count,kind=rprec))*real(c_count,kind=rprec)
!do k=1,Nz-1
!do i=1,Nx
!   tu1=0._rprec;tv1=0._rprec;tw1=0._rprec;tp1=0._rprec
!   ttxx=0._rprec;ttxz=0._rprec;ttyy=0._rprec;ttyz=0._rprec
!   ttzz=0._rprec;tdudz=0._rprec;tdvdz=0._rprec;tu2=0._rprec
!   tv2=0._rprec;tw2=0._rprec;tp2=0._rprec;tuw=0._rprec;tvw=0._rprec
!   tCs=0._rprec;tCs_Ssim=0._rprec;tu3=0._rprec;tv3=0._rprec;tw3=0._rprec
!
!   do j=1,Ny
!      tu1=tu1+u(i,j,k)
!      tv1=tv1+v(i,j,k)
!      tw1=tw1+w(i,j,k)
!      tp1=tp1+p(i,j,k)
!      ttxx=ttxx+txx(i,j,k)
!      ttxz=ttxz+txz(i,j,k)
!      ttyy=ttyy+tyy(i,j,k)
!      ttyz=ttyz+tyz(i,j,k)
!      ttzz=ttzz+tzz(i,j,k)
!      tdudz=tdudz+dudz(i,j,k)
!      tdvdz=tdvdz+dvdz(i,j,k)
!      tu2=tu2+u(i,j,k)*u(i,j,k)
!      tv2=tv2+v(i,j,k)*v(i,j,k)
!      tw2=tw2+w(i,j,k)*w(i,j,k)
!      tp2=tp2+p(i,j,k)*p(i,j,k)
!      tCs=tCs+sqrt(Cs_opt2(i,j,k))
!      tCs_Ssim=tCs_Ssim+sqrt(Cs_Ssim(i,j,k))
!      if((k .eq. 1) .AND. ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0))) then
!         arg1=0._rprec
!         arg2=0._rprec
!      else
!         arg1=(u(i,j,k)+u(i,j,k-1))/2.
!         arg2=(v(i,j,k)+v(i,j,k-1))/2.
!      end if
!      tuw=tuw+w(i,j,k)*arg1
!      tvw=tvw+w(i,j,k)*arg2
!      tu3=tu3+u(i,j,k)*u(i,j,k)*u(i,j,k)
!      tv3=tv3+v(i,j,k)*v(i,j,k)*v(i,j,k)
!      tw3=tw3+w(i,j,k)*w(i,j,k)*w(i,j,k)
!   end do
!   au(i,k)=au(i,k)+(fr)*tu1/Ny
!   av(i,k)=av(i,k)+(fr)*tv1/Ny
!   aw(i,k)=aw(i,k)+(fr)*tw1/Ny
!   ap(i,k)=ap(i,k)+(fr)*tp1/Ny
!   adudz(i,k)=adudz(i,k)+(fr)*tdudz/Ny
!   advdz(i,k)=advdz(i,k)+(fr)*tdvdz/Ny
!   u2(i,k)=u2(i,k)+(fr)*tu2/Ny
!   v2(i,k)=v2(i,k)+(fr)*tv2/Ny
!   w2(i,k)=w2(i,k)+(fr)*tw2/Ny
!   atxx(i,k)=atxx(i,k)+(fr)*ttxx/Ny
!   atxz(i,k)=atxz(i,k)+(fr)*ttxz/Ny
!   atyy(i,k)=atyy(i,k)+(fr)*ttyy/Ny
!   atyz(i,k)=atyz(i,k)+(fr)*ttyz/Ny
!   atzz(i,k)=atzz(i,k)+(fr)*ttzz/Ny
!   p2(i,k)=p2(i,k)+fr*tp2/Ny
!   aCs(i,k)=aCs(i,k)+(fr)*tCs/Ny
!   auw(i,k)=auw(i,k)+(fr)*tuw/Ny
!   avw(i,k)=avw(i,k)+(fr)*tvw/Ny
!   aCs_Ssim(i,k)=aCs_Ssim(i,k)+fr*tCs_Ssim/Ny
!   abeta_sgs(i,k)=abeta_sgs(i,k)+fr*Beta_avg(k)
!   abetaclip_sgs(i,k)=abetaclip_sgs(i,k)+fr*Betaclip_avg(k)
!   u3(i,k)=u3(i,k)+(fr)*tu3/Ny
!   v3(i,k)=v3(i,k)+(fr)*tv3/Ny
!   w3(i,k)=w3(i,k)+(fr)*tw3/Ny
!end do
!end do
!
!if (mod(jt,p_count)==0) then
!        allocate(avg_out(1:nx,1:(nz_tot-1)));
!        call collocate_MPI_averages_N(au,avg_out,20,'u')
!        call collocate_MPI_averages_N(av,avg_out,21,'v')
!        call collocate_MPI_averages_N(aw,avg_out,22,'w')
!        call collocate_MPI_averages_N(ap,avg_out,23,'p')
!        call collocate_MPI_averages_N(u2,avg_out,24,'u2')
!        call collocate_MPI_averages_N(v2,avg_out,25,'v2')
!        call collocate_MPI_averages_N(w2,avg_out,26,'w2')
!        call collocate_MPI_averages_N(p2,avg_out,32,'p2')
!        call collocate_MPI_averages_N(atxx,avg_out,27,'txx')
!        call collocate_MPI_averages_N(atxz,avg_out,28,'txz')
!        call collocate_MPI_averages_N(atyy,avg_out,29,'tyy')
!        call collocate_MPI_averages_N(atyz,avg_out,30,'tyz')
!        call collocate_MPI_averages_N(atzz,avg_out,31,'tzz')
!        call collocate_MPI_averages_N(auw,avg_out,33,'uw')
!        call collocate_MPI_averages_N(avw,avg_out,34,'vw')
!        call collocate_MPI_averages_N(aCs,avg_out,35,'Cs')
!        call collocate_MPI_averages_N(adudz,avg_out,36,'dudz')
!        call collocate_MPI_averages_N(advdz,avg_out,37,'dvdz')
!        call collocate_MPI_averages_N(aCs_Ssim,avg_out,38,'Cs_Ssim')
!        call collocate_MPI_averages_N(abeta_sgs,avg_out,39,'beta_sgs')
!        call collocate_MPI_averages_N(abetaclip_sgs,avg_out,40,'betaclip_sgs');
!        call collocate_MPI_averages_N(u3,avg_out,41,'u3')
!        call collocate_MPI_averages_N(v3,avg_out,42,'v3')
!        call collocate_MPI_averages_N(w3,avg_out,43,'w3');deallocate(avg_out)
!
!!VK Zero out the outputted averages !!
!        au=0._rprec;av=0._rprec;aw=0._rprec;ap=0._rprec;u2=0._rprec;v2=0._rprec
!        w2=0._rprec;atxx=0._rprec;atxz=0._rprec;atyy=0._rprec;atyz=0._rprec
!        atzz=0._rprec;p2=0._rprec;auw=0._rprec;avw=0._rprec;aCs=0._rprec
!        adudz=0._rprec;advdz=0._rprec;aCs_Ssim=0._rprec;abeta_sgs=0._rprec
!        abetaclip_sgs=0._rprec;u3=0._rprec;v3=0._rprec;w3=0._rprec;
!end if
!end subroutine avgslice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!xxxxxxxxxx----------MM----------XXXXXXXXXXXXXXXXXXXXX
!!!xxxxxxxxxx----------MM----------XXXXXXXXXXXXXXXXXXXXX
subroutine MM_budget_slice
use sim_param,only:path,u,v,w,dudz,dvdz,txx,txz,tyy,tyz,tzz,p,txy,&
         dudx,dudy,dvdx,dvdy,dwdx,dwdy,dwdz,L11t,L22t,L33t,&
         Q11t,Q22t,Q33t
use param,only:dz,p_count,c_count,jt
use sgsmodule,only:Cs_opt2,Cs_Ssim,Beta_avg,Betaclip_avg,dissip,tau !,SS,SS11,SS12,SS13,SS22,SS23,SS33,CS_O,Sl
implicit none
integer::i,j,k
real(kind=rprec),dimension(nx,nz-1),save::ap,au,av,aw,p2,u2,v2,w2,auw,avw,auv,atke,au2w,av2w,awwu,awwv,awe,adissip,u3,v3,w3
real(kind=rprec),dimension(nx,nz-1),save::aCs_Ssim,aCs,asig,aL11f,aL22f,aL33f,asigQ,aQ11f,aQ22f,aQ33f
real(kind=rprec),dimension(nx,nz-1),save::atxx,atxy,atxz,atyy,atyz,atzz,autxx,autxy,autxz,avtxy,avtyy,avtyz,awtxz,awtzz,awtyz
real(kind=rprec),dimension(nx,nz-1),save::aup,avp,awp,ausig,avsig,awsig,asigux,asiguy,asiguz,asigvx,asigvy,asigvz,asigwx,asigwy,asigwz,at11ux,at12uy,at13uz
real(kind=rprec),dimension(nx,nz-1),save::at21vx,at22vy,at23vz,at31wx,at32wy,at33wz
real(kind=rprec),dimension(nx,nz-1),save::apux,apuy,apuz,apvx,apvy,apvz,apwx,apwy,apwz,apux2,apuy2,apuz2,apvx2,apvy2,apvz2,apwx2,apwy2,apwz2
real(kind=rprec),dimension(nx,nz-1),save::at11s11,at12s12,at13s13,at22s22,at23s23,at33s33
real(kind=rprec),dimension(nx,nz-1),save::aupf,awpf,avpf,ap_c,ap_c2,auf,avf,awf
real(kind=rprec),dimension(nx,nz-1),save::ausigQ,avsigQ,awsigQ,asigQux,asigQuy,asigQuz,asigQvx,asigQvy,asigQvz,asigQwx,asigQwy,asigQwz
real(kind=rprec),dimension(nx,nz-1),save::autzz,at11wx,at12wy,at13wz,at31ux,at32uy,at33uz,as11,as12,as13,as22,as23,as33
real(kind=rprec),dimension(nx,nz-1),save::adudx,adudy,adudz,advdx,advdy,advdz,adwdx,adwdy,adwdz,adivu,apdivu,apfdivu
real(kind=rprec),dimension(nx,ny,nz-1)::SS11,SS12,SS22,SS13,SS23,SS33


real(kind=rprec)::tp,tu,tv,tw,tp2,tu2,tv2,tw2,tuw,tvw,tuv,ttke,tu2w,tv2w,twwu,twwv,twe,tdissip,tu3,tv3,tw3
real(kind=rprec)::tsig,tL11f,tL22f,tL33f,tsigQ,tQ11f,tQ22f,tQ33f
real(kind=rprec)::ttxx,ttxz,ttxy,ttyy,ttyz,ttzz,tutxx,tutxy,tutxz,tvtxy,tvtyy,tvtyz,twtxz,twtzz,twtyz
real(kind=rprec)::tup,tvp,twp,tusig,tvsig,twsig,tsigux,tsiguy,tsiguz,tsigvx,tsigvy,tsigvz,tsigwx,tsigwy,tsigwz,tt11ux,tt12uy,tt13uz
real(kind=rprec)::tt21vx,tt22vy,tt23vz,tt31wx,tt32wy,tt33wz
real(kind=rprec)::tpux,tpuy,tpuz,tpvx,tpvy,tpvz,tpwx,tpwy,tpwz,tpux2,tpuy2,tpuz2,tpvx2,tpvy2,tpvz2,tpwx2,tpwy2,tpwz2
real(kind=rprec)::tt11s11,tt12s12,tt13s13,tt22s22,tt23s23,tt33s33
real(kind=rprec)::tupf,twpf,tvpf,tp_c,tp_c2
real(kind=rprec)::tusigQ,tvsigQ,twsigQ,tsigQux,tsigQuy,tsigQuz,tsigQvx,tsigQvy,tsigQvz,tsigQwx,tsigQwy,tsigQwz
real(kind=rprec)::arg1,arg2,arg3,arg3f,arg4,arg5,arg6,arg7,arg8,arg9,arg10,sig1,sig1q
real(kind=rprec)::tutzz,tt11wx,tt12wy,tt13wz,tt31ux,tt32uy,tt33uz,ts11,ts12,ts13,ts22,ts23,ts33  
real(kind=rprec)::tCs_Ssim,tCs,fr,tuf,tvf,twf
real(kind=rprec)::tdudx,tdudy,tdudz,tdvdx,tdvdy,tdvdz,tdwdx,tdwdy,tdwdz,tdivu,tpdivu,tpfdivu


real(kind=rprec),dimension(:,:),allocatable::avg_out

real(kind=rprec),dimension(nx,ny),save::SD
real(kind=rprec)::S11,S22,S33,S12,S13,S23,&
     ux,uy,uz,vx,vy,vz,wx,wy,wz

real(rprec),parameter::delbar=2._rprec
real(rprec)::denomsig=0.5_rprec/((delbar**(2.0_rprec/3.0_rprec))-1._rprec)
real(rprec)::denomsigq=0.5_rprec/(((2._rprec*delbar)**(2.0_rprec/3.0_rprec))-1._rprec)



fr=(1._rprec/real(p_count,kind=rprec))*real(c_count,kind=rprec)


do k=1,Nz-1
do i=1,Nx

   tu=0._rprec;tv=0._rprec;tw=0._rprec;tp=0._rprec;tp2=0._rprec;tu2=0._rprec;tv2=0._rprec;tw2=0._rprec;
   tuw=0._rprec;tvw=0._rprec;tuv=0._rprec;ttke=0._rprec;tu2w=0._rprec;tv2w=0._rprec;
   twe=0._rprec;tdissip=0._rprec;tu3=0._rprec;tv3=0._rprec;tw3=0._rprec;
   tCs_Ssim=0._rprec;tsig=0._rprec;
   tsigQ=0._rprec;tCs=0._rprec;tQ11f=0._rprec;tQ22f=0._rprec;tQ33f=0._rprec;
   tL11f=0._rprec;tL22f=0._rprec;tL33f=0._rprec;ttxx=0._rprec;ttxy=0._rprec;ttxz=0._rprec;ttyy=0._rprec;ttyz=0._rprec;
   ttzz=0._rprec;tutxx=0._rprec;tutxy=0._rprec;tutxz=0._rprec;
   tvtxy=0._rprec;tvtyy=0._rprec;tvtyz=0._rprec;
   twtxz=0._rprec;twtzz=0._rprec;twtyz=0._rprec;
   tup=0._rprec;tvp=0._rprec;twp=0._rprec;tusig=0._rprec;tvsig=0._rprec;twsig=0._rprec;tsigux=0._rprec;
   tsiguy=0._rprec;tsiguz=0._rprec;tsigvx=0._rprec;tsigvy=0._rprec;tsigvz=0._rprec;tsigwx=0._rprec;
   tsigwy=0._rprec;tsigwz=0._rprec;tt11ux=0._rprec;tt12uy=0._rprec;tt13uz=0._rprec;tt21vx=0._rprec;
   tt22vy=0._rprec;tt23vz=0._rprec;tt31wx=0._rprec;tt32wy=0._rprec;tt33wz=0._rprec;
   tpux=0._rprec;tpuy=0._rprec;tpuz=0._rprec;tpvx=0._rprec;tpvy=0._rprec;tpvz=0._rprec;
   tpwx=0._rprec;tpwy=0._rprec;
   tpwz=0._rprec;tpux2=0._rprec;tpuy2=0._rprec;tpuz2=0._rprec;tpvx2=0._rprec;tpvy2=0._rprec;tpvz2=0._rprec;
   tpwx2=0._rprec;tpwy2=0._rprec;tpwz2=0._rprec;
   tt11s11=0._rprec;tt12s12=0._rprec;tt13s13=0._rprec;tt22s22=0._rprec;
   tt23s23=0._rprec;tt33s33=0._rprec;
   tupf=0._rprec;twpf=0._rprec;tvpf=0._rprec;
   tCs_Ssim=0._rprec;
   tusigQ=0._rprec;tvsigQ=0._rprec;twsigQ=0._rprec;tsigQux=0._rprec;
   tsigQuy=0._rprec;tsigQuz=0._rprec;tsigQvx=0._rprec;tsigQvy=0._rprec;tsigQvz=0._rprec;tsigQwx=0._rprec;
   tsigQwy=0._rprec;tsigQwz=0._rprec;tp_c=0._rprec;tp_c2=0._rprec;
   twwu=0._rprec;twwv=0._rprec;tutzz=0._rprec;tt11wx=0._rprec;tt12wy=0._rprec;tt13wz=0._rprec;
   tt31ux=0._rprec;tt32uy=0._rprec;tt33uz=0._rprec; 
   ts11=0._rprec;ts12=0._rprec;ts13=0._rprec;ts22=0._rprec;ts23=0._rprec;ts33=0._rprec; 
   tdudx=0._rprec;tdudy=0._rprec;tdudz=0._rprec;tdvdx=0._rprec;tdvdy=0._rprec;
   tdvdz=0._rprec;tdwdx=0._rprec;tdwdy=0._rprec;tdwdz=0._rprec;
   tdivu=0._rprec;tpdivu=0._rprec;tpfdivu=0._rprec;tuf=0._rprec;tvf=0._rprec;twf=0._rprec;


  do j=1,Ny

! Put all variables on cs nodes (i.e. at dz/2 [uvp node] for the 1st point and at dz,2dz,et.[w nodes] everywhere else)

        if((k .eq. 1) .AND. ((.not. USE_MPI) .or. (USE_MPI .and. coord .eq. 0))) then
        
         arg1=u(i,j,k)
         arg2=v(i,j,k)
         arg4=0.25_rprec*w(i,j,k+1)
         arg3=p(i,j,k)-0.5_rprec*(arg1*arg1+arg2*arg2+arg4*arg4)
! KMG: Lii and Qii are already on cs nodes 
         arg3f= p(i,j,k)- 0.5_rprec*(arg1*arg1+arg2*arg2+arg4*arg4)-(1._rprec/3._rprec)*(denomsig)*(L11t(i,j,k)+L22t(i,j,k)+L33t(i,j,k))
!         arg5=txx(i,j,k)
!         arg6=txy(i,j,k)
!         arg7=tyy(i,j,k)
!         arg8=tzz(i,j,k)
!         arg9=(txz(i,j,k)+txz(i,j,k+1))/2.0_rprec
!         arg10=(tyz(i,j,k)+tyz(i,j,k+1))/2.0_rprec
         arg9=txz(i,j,k)
         arg10=tyz(i,j,k)
         tp=tp+p(i,j,k)
         ux=dudx(i,j,k)
         uy=dudy(i,j,k)
         vx=dvdx(i,j,k)
         vy=dvdy(i,j,k)
         wz=dwdz(i,j,k)
         wx=(dwdx(i,j,k)+dwdx(i,j,k+1))/2.0_rprec
         wy=(dwdy(i,j,k)+dwdy(i,j,k+1))/2.0_rprec
 
         else
        
         arg1=(u(i,j,k)+u(i,j,k-1))/2.0_rprec
         arg2=(v(i,j,k)+v(i,j,k-1))/2.0_rprec
         arg4=w(i,j,k)
         arg3=(p(i,j,k)+p(i,j,k-1))/2.0_rprec - 0.5_rprec*(arg1*arg1+arg2*arg2+arg4*arg4)
         arg3f=(p(i,j,k)+p(i,j,k-1))/2.0_rprec - 0.5_rprec*(arg1*arg1+arg2*arg2+arg4*arg4)-(1._rprec/3._rprec)*(denomsig)*(L11t(i,j,k)+L22t(i,j,k)+L33t(i,j,k))
!         arg5=(txx(i,j,k)+txx(i,j,k-1))/2.0_rprec
!         arg6=(txy(i,j,k)+txy(i,j,k-1))/2.0_rprec
!         arg7=(tyy(i,j,k)+tyy(i,j,k-1))/2.0_rprec
!         arg8=(tzz(i,j,k)+tzz(i,j,k-1))/2.0_rprec
         arg9=txz(i,j,k)
         arg10=tyz(i,j,k)
         tp=tp+(p(i,j,k)+p(i,j,k-1))/2.0_rprec
         ux=(dudx(i,j,k)+dudx(i,j,k-1))/2.0_rprec
         uy=(dudy(i,j,k)+dudy(i,j,k-1))/2.0_rprec
         vx=(dvdx(i,j,k)+dvdx(i,j,k-1))/2.0_rprec
         vy=(dvdy(i,j,k)+dvdy(i,j,k-1))/2.0_rprec
         wz=(dwdz(i,j,k)+dwdz(i,j,k-1))/2.0_rprec
         wx=dwdx(i,j,k)
         wy=dwdy(i,j,k) 


        end if

      arg5=txx(i,j,k)
      arg6=txy(i,j,k)
      arg7=tyy(i,j,k)
      arg8=tzz(i,j,k)
 
      tu=tu+arg1
      tv=tv+arg2
      tw=tw+arg4
      ttxx=ttxx+arg5
      ttyy=ttyy+arg7
      ttzz=ttzz+arg8
      ttxy=ttxy+arg6
      

      tuf=tuf+u(i,j,k)
      tvf=tvf+v(i,j,k)
      twf=twf+w(i,j,k)



     ttxz=ttxz+arg9
     ttyz=ttyz+arg10
     uz=dudz(i,j,k)
     vz=dvdz(i,j,k) 

 
      tu2=tu2+arg1*arg1
      tv2=tv2+arg2*arg2
      tw2=tw2+arg4*arg4
      tp_c=tp_c + arg3
      tp_c2=tp_c2 + arg3f
      tp2=tp2+arg3f*arg3f
      tCs=tCs+sqrt(Cs_opt2(i,j,k))
      tCs_Ssim=tCs_Ssim+sqrt(Cs_Ssim(i,j,k))
      
        S11=ux 
        S12=0.5_rprec*(uy+vx) 
        S13=0.5_rprec*(uz+wx)
        S22=vy
        S23=0.5_rprec*(vz+wy) 
        S33=wz 
       
      SS11(i,j,k)=S11
      SS12(i,j,k)=S12
      SS13(i,j,k)=S13
      SS22(i,j,k)=S22
      SS23(i,j,k)=S23
      SS33(i,j,k)=S33  
      
      tdudx=tdudx+ux
      tdudy=tdudy+uy
      tdudz=tdudz+uz
      tdvdx=tdvdx+vx
      tdvdy=tdvdy+vy
      tdvdz=tdvdz+vz
      tdwdx=tdwdx+wx
      tdwdy=tdwdy+wy
      tdwdz=tdwdz+wz

   
      ts11=ts11+SS11(i,j,k)
      ts12=ts12+SS12(i,j,k)
      ts13=ts13+SS13(i,j,k)
      ts22=ts22+SS22(i,j,k)
      ts23=ts23+SS23(i,j,k)
      ts33=ts33+SS33(i,j,k)

      tuw=tuw+arg1*arg4
      tvw=tvw+arg2*arg4
      tuv=tuv+arg1*arg2

      tu3=tu3+arg1*arg1*arg1
      tv3=tv3+arg2*arg2*arg2
      tw3=tw3+arg4*arg4*arg4
     
      tu2w=tu2w+arg1*arg1*arg4
      tv2w=tv2w+arg2*arg2*arg4
      twwu=twwu+arg1*arg4*arg4
      twwv=twwv+arg2*arg4*arg4 
     

      tpux=tpux +  arg3*ux
      tpuy=tpuy +  arg3*uy
      tpuz=tpuz +  arg3*uz
      
      tpvx=tpvx + arg3*vx
      tpvy=tpvy + arg3*vy
      tpvz=tpvz + arg3*vz


      tpwx=tpwx + arg3*wx
      tpwy=tpwy + arg3*wy
      tpwz=tpwz + arg3*wz

      tpux2=tpux2 +  arg3f*ux
      tpuy2=tpuy2 +  arg3f*uy
      tpuz2=tpuz2 +  arg3f*uz

      tpvx2=tpvx2 + arg3f*vx
      tpvy2=tpvy2 + arg3f*vy
      tpvz2=tpvz2 + arg3f*vz


      tpwx2=tpwx2 + arg3f*wx
      tpwy2=tpwy2 + arg3f*wy
      tpwz2=tpwz2 + arg3f*wz


    tdivu=tdivu+(ux+vy+wz)
    tpdivu=tpdivu+arg3*(ux+vy+wz)
    tpfdivu=tpfdivu+arg3f*(ux+vy+wz)


      tup=tup+arg3*arg1
      tupf=tupf+arg3f*arg1

      tvp=tvp+arg3*arg2
      tvpf=tvpf+arg3f*arg2


      twp=twp+arg4*arg3
      twpf=twpf+arg4*arg3f
   
      
      tdissip=tdissip+dissip(i,j,k)
      ttke=ttke+0.5_rprec*(arg1*arg1+arg2*arg2+arg4*arg4)
      twe=twe+0.5_rprec*arg4*(arg1*arg1+arg2*arg2+arg4*arg4) 
      

      sig1=denomsig*(L11t(i,j,k)+L22t(i,j,k)+L33t(i,j,k))
      sig1q=denomsigq*(Q11t(i,j,k)+Q22t(i,j,k)+Q33t(i,j,k))

      tsig=tsig+sig1
      tsigQ=tsigQ+sig1q
      tL11f=tL11f+L11t(i,j,k)
      tL22f=tL22f+L22t(i,j,k)
      tL33f=tL33f+L33t(i,j,k)
      tQ11f=tQ11f+Q11t(i,j,k)
      tQ22f=tQ22f+Q22t(i,j,k)
      tQ33f=tQ33f+Q33t(i,j,k)    

      tusig=tusig+arg1*sig1
      tvsig=tvsig+arg2*sig1
      twsig=twsig+arg4*sig1


     tsigux=tsigux+ux*sig1
     tsiguy=tsiguy+uy*sig1
     tsiguz=tsiguz+uz*sig1

     tsigvx=tsigvx+vx*sig1
     tsigvy=tsigvy+vy*sig1
     tsigvz=tsigvz+vz*sig1


     tsigwx=tsigwx+wx*sig1
     tsigwy=tsigwy+wy*sig1
     tsigwz=tsigwz+wz*sig1

     tusigQ=tusigQ+arg1*sig1q
     tvsigQ=tvsigQ+arg2*sig1q
     twsigQ=twsigQ+arg4*sig1q


     tsigQux=tsigQux+ux*sig1q
     tsigQuy=tsigQuy+uy*sig1q
     tsigQuz=tsigQuz+uz*sig1q

     tsigQvx=tsigQvx+vx*sig1q
     tsigQvy=tsigQvy+vy*sig1q
     tsigQvz=tsigQvz+vz*sig1q


     tsigQwx=tsigQwx+wx*sig1q
     tsigQwy=tsigQwy+wy*sig1q
     tsigQwz=tsigQwz+wz*sig1q


     tt11ux=tt11ux+ux*arg5  ! txx*du/dx
     tt12uy=tt12uy+uy*arg6  ! txy*du/dy
     tt13uz=tt13uz+uz*arg9
     tt21vx=tt21vx+vx*arg6
     tt22vy=tt22vy+vy*arg7
     tt23vz=tt23vz+vz*arg10
     tt31wx=tt31wx+wx*arg9
     tt32wy=tt32wy+wy*arg10
     tt33wz=tt33wz+wz*arg8
     tt11wx=tt11wx+wx*arg5
     tt12wy=tt12wy+wy*arg6
     tt13wz=tt13wz+wz*arg9
     tt31ux=tt31ux+ux*arg9
     tt32uy=tt32uy+uy*arg10
     tt33uz=tt33uz+uz*arg8

      tutxx=tutxx+arg1*arg5
      tutxy=tutxy+arg1*arg6
      tutxz=tutxz+arg1*arg9
      tvtxy=tvtxy+arg2*arg6
      tvtyy=tvtyy+arg2*arg7
      tvtyz=tvtyz+arg2*arg10
      twtxz=twtxz+arg4*arg9
      twtyz=twtyz+arg4*arg10
      twtzz=twtzz+arg4*arg8
      tutzz=tutzz+arg1*arg8
      
      tt11s11=tt11s11+ arg5*SS11(i,j,k)
      tt12s12=tt12s12+ arg6*SS12(i,j,k)
      tt13s13=tt13s13+ arg9*SS13(i,j,k)
      tt22s22=tt22s22+ arg7*SS22(i,j,k)
      tt23s23=tt23s23+ arg10*SS23(i,j,k)  
      tt33s33=tt33s33+ arg8*SS33(i,j,k)

 end do

   
   au(i,k)=au(i,k)+(fr)*tu/Ny
   av(i,k)=av(i,k)+(fr)*tv/Ny
   aw(i,k)=aw(i,k)+(fr)*tw/Ny
   auf(i,k)=auf(i,k)+(fr)*tuf/Ny
   avf(i,k)=avf(i,k)+(fr)*tvf/Ny
   awf(i,k)=awf(i,k)+(fr)*twf/Ny
   ap(i,k)=ap(i,k)+(fr)*tp/Ny
   u2(i,k)=u2(i,k)+(fr)*tu2/Ny 
   v2(i,k)=v2(i,k)+(fr)*tv2/Ny
   w2(i,k)=w2(i,k)+(fr)*tw2/Ny
   p2(i,k)=p2(i,k)+(fr)*tp2/Ny
   auw(i,k)=auw(i,k)+(fr)*tuw/Ny
   avw(i,k)=avw(i,k)+(fr)*tvw/Ny
   auv(i,k)=auv(i,k)+(fr)*tuv/Ny
   atke(i,k)=atke(i,k)+(fr)*ttke/Ny
   au2w(i,k)=au2w(i,k)+(fr)*tu2w/Ny
   av2w(i,k)=av2w(i,k)+(fr)*tv2w/Ny
   awwu(i,k)=awwu(i,k)+(fr)*twwu/Ny
   awwv(i,k)=awwv(i,k)+(fr)*twwv/Ny
   awe(i,k)=awe(i,k)+(fr)*twe/Ny
   adissip(i,k)=adissip(i,k)+(fr)*tdissip/Ny
   u3(i,k)=u3(i,k)+(fr)*tu3/Ny
   v3(i,k)=v3(i,k)+(fr)*tv3/Ny
   w3(i,k)=w3(i,k)+(fr)*tw3/Ny
   aCs_Ssim(i,k)=aCs_Ssim(i,k)+(fr)*tCs_Ssim/Ny
   aCs(i,k)=aCs(i,k)+(fr)*tCs/Ny
   asig(i,k)=asig(i,k)+(fr)*tsig/Ny
   aL11f(i,k)=aL11f(i,k)+(fr)*tL11f/Ny
   aL22f(i,k)=aL22f(i,k)+(fr)*tL22f/Ny  
   aL33f(i,k)=aL33f(i,k)+(fr)*tL33f/Ny
   aQ11f(i,k)=aQ11f(i,k)+(fr)*tQ11f/Ny
   aQ22f(i,k)=aQ22f(i,k)+(fr)*tQ22f/Ny
   aQ33f(i,k)=aQ33f(i,k)+(fr)*tQ33f/Ny
   asigQ(i,k)=asigQ(i,k)+(fr)*tsigQ/Ny
   atxx(i,k)=atxx(i,k)+(fr)*ttxx/Ny
   atxy(i,k)=atxy(i,k)+(fr)*ttxy/Ny
   atxz(i,k)=atxz(i,k)+(fr)*ttxz/Ny
   atyy(i,k)=atyy(i,k)+(fr)*ttyy/Ny
   atyz(i,k)=atyz(i,k)+(fr)*ttyz/Ny
   atzz(i,k)=atzz(i,k)+(fr)*ttzz/Ny
   autxx(i,k)=autxx(i,k)+(fr)*tutxx/Ny
   autxy(i,k)=autxy(i,k)+(fr)*tutxy/Ny
   autxz(i,k)=autxz(i,k)+(fr)*tutxz/Ny  
   avtxy(i,k)=avtxy(i,k)+(fr)*tvtxy/Ny
   avtyy(i,k)=avtyy(i,k)+(fr)*tvtyy/Ny
   avtyz(i,k)=avtyz(i,k)+(fr)*tvtyz/Ny
   awtxz(i,k)=awtxz(i,k)+(fr)*twtxz/Ny
   awtyz(i,k)=awtyz(i,k)+(fr)*twtyz/Ny
   awtzz(i,k)=awtzz(i,k)+(fr)*twtzz/Ny
   aup(i,k)=aup(i,k)+(fr)*tup/Ny
   avp(i,k)=avp(i,k)+(fr)*tvp/Ny
   awp(i,k)=awp(i,k)+(fr)*twp/Ny
   ausig(i,k)=ausig(i,k)+(fr)*tusig/Ny
   avsig(i,k)=avsig(i,k)+(fr)*tvsig/Ny
   awsig(i,k)=awsig(i,k)+(fr)*twsig/Ny
   asigux(i,k)=asigux(i,k)+(fr)*tsigux/Ny
   asiguy(i,k)=asiguy(i,k)+(fr)*tsiguy/Ny
   asiguz(i,k)=asiguz(i,k)+(fr)*tsiguz/Ny
   asigvx(i,k)=asigvx(i,k)+(fr)*tsigvx/Ny
   asigvy(i,k)=asigvy(i,k)+(fr)*tsigvy/Ny
   asigvz(i,k)=asigvz(i,k)+(fr)*tsigvz/Ny
   asigwx(i,k)=asigwx(i,k)+(fr)*tsigwx/Ny
   asigwy(i,k)=asigwy(i,k)+(fr)*tsigwy/Ny
   asigwz(i,k)=asigwz(i,k)+(fr)*tsigwz/Ny
   at11ux(i,k)=at11ux(i,k)+(fr)*tt11ux/Ny
   at12uy(i,k)=at12uy(i,k)+(fr)*tt12uy/Ny
   at13uz(i,k)=at13uz(i,k)+(fr)*tt13uz/Ny
   at21vx(i,k)=at21vx(i,k)+(fr)*tt21vx/Ny
   at22vy(i,k)=at22vy(i,k)+(fr)*tt22vy/Ny
   at23vz(i,k)=at23vz(i,k)+(fr)*tt23vz/Ny
   at31wx(i,k)=at31wx(i,k)+(fr)*tt31wx/Ny
   at32wy(i,k)=at32wy(i,k)+(fr)*tt32wy/Ny
   at33wz(i,k)=at33wz(i,k)+(fr)*tt33wz/Ny
   apux(i,k)=apux(i,k)+(fr)*tpux/Ny
   apuy(i,k)=apuy(i,k)+(fr)*tpuy/Ny
   apuz(i,k)=apuz(i,k)+(fr)*tpuz/Ny
   apvx(i,k)=apvx(i,k)+(fr)*tpvx/Ny
   apvy(i,k)=apvy(i,k)+(fr)*tpvy/Ny
   apvz(i,k)=apvz(i,k)+(fr)*tpvz/Ny
   apwx(i,k)=apwx(i,k)+(fr)*tpwx/Ny
   apwy(i,k)=apwy(i,k)+(fr)*tpwy/Ny
   apwz(i,k)=apwz(i,k)+(fr)*tpwz/Ny
   apux2(i,k)=apux2(i,k)+(fr)*tpux2/Ny
   apuy2(i,k)=apuy2(i,k)+(fr)*tpuy2/Ny
   apuz2(i,k)=apuz2(i,k)+(fr)*tpuz2/Ny
   apvx2(i,k)=apvx2(i,k)+(fr)*tpvx2/Ny
   apvy2(i,k)=apvy2(i,k)+(fr)*tpvy2/Ny
   apvz2(i,k)=apvz2(i,k)+(fr)*tpvz2/Ny
   apwx2(i,k)=apwx2(i,k)+(fr)*tpwx2/Ny
   apwy2(i,k)=apwy2(i,k)+(fr)*tpwy2/Ny
   apwz2(i,k)=apwz2(i,k)+(fr)*tpwz2/Ny
   at11s11(i,k)=at11s11(i,k)+(fr)*tt11s11/Ny
   at12s12(i,k)=at12s12(i,k)+(fr)*tt12s12/Ny
   at13s13(i,k)=at13s13(i,k)+(fr)*tt13s13/Ny
   at22s22(i,k)=at22s22(i,k)+(fr)*tt22s22/Ny
   at23s23(i,k)=at23s23(i,k)+(fr)*tt23s23/Ny
   at33s33(i,k)=at33s33(i,k)+(fr)*tt33s33/Ny
   aupf(i,k)=aupf(i,k)+(fr)*tupf/Ny
   avpf(i,k)=avpf(i,k)+(fr)*tvpf/Ny
   awpf(i,k)=awpf(i,k)+(fr)*twpf/Ny
   ap_c(i,k)=ap_c(i,k)+(fr)*tp_c/Ny
   ap_c2(i,k)=ap_c2(i,k)+(fr)*tp_c2/Ny
   ausigQ(i,k)=ausigQ(i,k)+(fr)*tusigQ/Ny
   avsigQ(i,k)=avsigQ(i,k)+(fr)*tvsigQ/Ny
   awsigQ(i,k)=awsigQ(i,k)+(fr)*twsigQ/Ny
   asigQux(i,k)=asigQux(i,k)+(fr)*tsigQux/Ny
   asigQuy(i,k)=asigQuy(i,k)+(fr)*tsigQuy/Ny
   asigQuz(i,k)=asigQuz(i,k)+(fr)*tsigQuz/Ny
   asigQvx(i,k)=asigQvx(i,k)+(fr)*tsigQvx/Ny
   asigQvy(i,k)=asigQvy(i,k)+(fr)*tsigQvy/Ny
   asigQvz(i,k)=asigQvz(i,k)+(fr)*tsigQvz/Ny
   asigQwx(i,k)=asigQwx(i,k)+(fr)*tsigQwx/Ny
   asigQwy(i,k)=asigQwy(i,k)+(fr)*tsigQwy/Ny
   asigQwz(i,k)=asigQwz(i,k)+(fr)*tsigQwz/Ny
   autzz(i,k)=autzz(i,k)+(fr)*tutzz/Ny
   at11wx(i,k)=at11wx(i,k)+(fr)*tt11wx/Ny
   at12wy(i,k)=at12wy(i,k)+(fr)*tt12wy/Ny
   at13wz(i,k)=at13wz(i,k)+(fr)*tt13wz/Ny
   at31ux(i,k)=at31ux(i,k)+(fr)*tt31ux/Ny
   at32uy(i,k)=at32uy(i,k)+(fr)*tt32uy/Ny
   at33uz(i,k)=at33uz(i,k)+(fr)*tt33uz/Ny


   as11(i,k)=as11(i,k)+(fr)*ts11/Ny
   as12(i,k)=as12(i,k)+(fr)*ts12/Ny
   as13(i,k)=as13(i,k)+(fr)*ts13/Ny  
   as22(i,k)=as22(i,k)+(fr)*ts22/Ny
   as23(i,k)=as23(i,k)+(fr)*ts23/Ny
   as33(i,k)=as33(i,k)+(fr)*ts33/Ny 

  adudx(i,k)=adudx(i,k)+(fr)*tdudx/Ny
  adudy(i,k)=adudy(i,k)+(fr)*tdudy/Ny
  adudz(i,k)=adudz(i,k)+(fr)*tdudz/Ny
  advdx(i,k)=advdx(i,k)+(fr)*tdvdx/Ny
  advdy(i,k)=advdy(i,k)+(fr)*tdvdy/Ny
  advdz(i,k)=advdz(i,k)+(fr)*tdvdz/Ny
  adwdx(i,k)=adwdx(i,k)+(fr)*tdwdx/Ny
  adwdy(i,k)=adwdy(i,k)+(fr)*tdwdy/Ny
  adwdz(i,k)=adwdz(i,k)+(fr)*tdwdz/Ny
  adivu(i,k)=adivu(i,k)+(fr)*tdivu/Ny
  apdivu(i,k)=apdivu(i,k)+(fr)*tpdivu/Ny
  apfdivu(i,k)=apfdivu(i,k)+(fr)*tpfdivu/Ny

end do
end do


if (mod(jt,p_count)==0) then
! file number 35 till 75 reserved for scalar_slice and scalar_TKE_budget in scalars_module2.f90
        allocate(avg_out(1:nx,1:(nz_tot-1)));
        call collocate_MPI_averages_N(au,avg_out,20,'u')
        call collocate_MPI_averages_N(av,avg_out,21,'v')
        call collocate_MPI_averages_N(aw,avg_out,22,'w')
        call collocate_MPI_averages_N(ap,avg_out,23,'p')
        call collocate_MPI_averages_N(u2,avg_out,24,'u2')
        call collocate_MPI_averages_N(v2,avg_out,25,'v2')
        call collocate_MPI_averages_N(w2,avg_out,26,'w2')
        call collocate_MPI_averages_N(p2,avg_out,27,'p2')
        call collocate_MPI_averages_N(atxx,avg_out,28,'txx')
        call collocate_MPI_averages_N(atxz,avg_out,29,'txz')
        call collocate_MPI_averages_N(atyy,avg_out,30,'tyy')
        call collocate_MPI_averages_N(atyz,avg_out,31,'tyz')
        call collocate_MPI_averages_N(atzz,avg_out,32,'tzz')
        call collocate_MPI_averages_N(atxy,avg_out,33,'txy')
        call collocate_MPI_averages_N(auw,avg_out,34,'uw')
        call collocate_MPI_averages_N(avw,avg_out,76,'vw')
        call collocate_MPI_averages_N(auv,avg_out,77,'uv')
        call collocate_MPI_averages_N(aCs,avg_out,78,'Cs')
        call collocate_MPI_averages_N(adudz,avg_out,79,'dudz')
        call collocate_MPI_averages_N(advdz,avg_out,80,'dvdz')
        call collocate_MPI_averages_N(aCs_Ssim,avg_out,81,'Cs_Ssim')
        call collocate_MPI_averages_N(u3,avg_out,82,'u3')
        call collocate_MPI_averages_N(v3,avg_out,83,'v3')
        call collocate_MPI_averages_N(w3,avg_out,84,'w3')
        
        call collocate_MPI_averages_N(atke,avg_out,181,'tke');
       
        call collocate_MPI_averages_N(au2w,avg_out,183,'au2w')
        call collocate_MPI_averages_N(av2w,avg_out,184,'av2w')
        call collocate_MPI_averages_N(awe,avg_out,185,'awe')
        call collocate_MPI_averages_N(adissip,avg_out,186,'dissip')
        call collocate_MPI_averages_N(asig,avg_out,187,'sigma')
        call collocate_MPI_averages_N(asigQ,avg_out,188,'sigmaQ')
        call collocate_MPI_averages_N(autxx,avg_out,189,'utxx')
        call collocate_MPI_averages_N(autxy,avg_out,190,'utxy')
        call collocate_MPI_averages_N(autxz,avg_out,191,'utxz')
        call collocate_MPI_averages_N(avtxy,avg_out,192,'vtxy')
        call collocate_MPI_averages_N(avtyy,avg_out,193,'vtyy')
        call collocate_MPI_averages_N(avtyz,avg_out,194,'vtyz')
        call collocate_MPI_averages_N(awtxz,avg_out,195,'wtxz')
        call collocate_MPI_averages_N(awtyz,avg_out,196,'wtyz')
        call collocate_MPI_averages_N(awtzz,avg_out,197,'wtzz')
        call collocate_MPI_averages_N(aup,avg_out,198,'up')
        call collocate_MPI_averages_N(avp,avg_out,199,'vp')
        call collocate_MPI_averages_N(awp,avg_out,200,'wp')
        call collocate_MPI_averages_N(ausig,avg_out,201,'usig')
        call collocate_MPI_averages_N(avsig,avg_out,202,'vsig')
        call collocate_MPI_averages_N(awsig,avg_out,203,'wsig')
        call collocate_MPI_averages_N(asigux,avg_out,204,'sigux')
        call collocate_MPI_averages_N(asiguy,avg_out,205,'siguy')
        call collocate_MPI_averages_N(asiguz,avg_out,206,'siguz')
        call collocate_MPI_averages_N(asigvx,avg_out,207,'sigvx')
        call collocate_MPI_averages_N(asigvy,avg_out,208,'sigvy')
        call collocate_MPI_averages_N(asigvz,avg_out,209,'sigvz')
        call collocate_MPI_averages_N(asigwx,avg_out,210,'sigwx')
        call collocate_MPI_averages_N(asigwy,avg_out,211,'sigwy')
        call collocate_MPI_averages_N(asigwz,avg_out,212,'sigwz')
        call collocate_MPI_averages_N(at11ux,avg_out,213,'t11ux')
        call collocate_MPI_averages_N(at12uy,avg_out,214,'t12uy')
        call collocate_MPI_averages_N(at13uz,avg_out,215,'t13uz')
        call collocate_MPI_averages_N(at21vx,avg_out,216,'t21vx')
        call collocate_MPI_averages_N(at22vy,avg_out,217,'t22vy')
        call collocate_MPI_averages_N(at23vz,avg_out,218,'t23vz')
        call collocate_MPI_averages_N(at31wx,avg_out,219,'t31wx')
        call collocate_MPI_averages_N(at32wy,avg_out,220,'t32wy')
        call collocate_MPI_averages_N(at33wz,avg_out,221,'t33wz')  
        call collocate_MPI_averages_N(apux,avg_out,222,'pux')
        call collocate_MPI_averages_N(apuy,avg_out,223,'puy')   
        call collocate_MPI_averages_N(apuz,avg_out,224,'puz')
        call collocate_MPI_averages_N(apvx,avg_out,225,'pvx')
        call collocate_MPI_averages_N(apvy,avg_out,226,'pvy')
        call collocate_MPI_averages_N(apvz,avg_out,227,'pvz')
        call collocate_MPI_averages_N(apwx,avg_out,228,'pwx')
        call collocate_MPI_averages_N(apwy,avg_out,229,'pwy')
        call collocate_MPI_averages_N(apwz,avg_out,230,'pwz')
        call collocate_MPI_averages_N(apux2,avg_out,231,'pux2')    
        call collocate_MPI_averages_N(apuy2,avg_out,232,'puy2')
        call collocate_MPI_averages_N(apuz2,avg_out,233,'puz2')
        call collocate_MPI_averages_N(apvx2,avg_out,234,'pvx2')
        call collocate_MPI_averages_N(apvy2,avg_out,235,'pvy2')
        call collocate_MPI_averages_N(apvz2,avg_out,236,'pvz2')
        call collocate_MPI_averages_N(apwx2,avg_out,237,'pwx2')
        call collocate_MPI_averages_N(apwy2,avg_out,238,'pwy2')
        call collocate_MPI_averages_N(apwz2,avg_out,239,'pwz2')
        call collocate_MPI_averages_N(at11s11,avg_out,240,'t11s11')
        call collocate_MPI_averages_N(at12s12,avg_out,241,'t12s12')  
        call collocate_MPI_averages_N(at13s13,avg_out,242,'t13s13')
        call collocate_MPI_averages_N(at22s22,avg_out,243,'t22s22')
        call collocate_MPI_averages_N(at23s23,avg_out,244,'t23s23')
        call collocate_MPI_averages_N(at33s33,avg_out,245,'t33s33')
        call collocate_MPI_averages_N(aupf,avg_out,246,'upf')
        call collocate_MPI_averages_N(avpf,avg_out,247,'vpf')
        call collocate_MPI_averages_N(awpf,avg_out,248,'wpf')
        call collocate_MPI_averages_N(ap_c,avg_out,249,'ap_c')
        call collocate_MPI_averages_N(ap_c2,avg_out,250,'ap_c2')
        call collocate_MPI_averages_N(ausigQ,avg_out,251,'usigQ')
        call collocate_MPI_averages_N(avsigQ,avg_out,252,'vsigQ')
        call collocate_MPI_averages_N(awsigQ,avg_out,253,'wsigQ')
        call collocate_MPI_averages_N(asigQux,avg_out,254,'sigQux')
        call collocate_MPI_averages_N(asigQuy,avg_out,255,'sigQuy')
        call collocate_MPI_averages_N(asigQuz,avg_out,256,'sigQuz')
        call collocate_MPI_averages_N(asigQvx,avg_out,257,'sigQvx')
        call collocate_MPI_averages_N(asigQvy,avg_out,258,'sigQvy')
        call collocate_MPI_averages_N(asigQvz,avg_out,259,'sigQvz')
        call collocate_MPI_averages_N(asigQwx,avg_out,260,'sigQwx')
        call collocate_MPI_averages_N(asigQwy,avg_out,261,'sigQwy')
        call collocate_MPI_averages_N(asigQwz,avg_out,262,'sigQwz')
        call collocate_MPI_averages_N(aL11f,avg_out,263,'L11')
        call collocate_MPI_averages_N(aL22f,avg_out,264,'L22')
        call collocate_MPI_averages_N(aL33f,avg_out,265,'L33')
        call collocate_MPI_averages_N(aQ11f,avg_out,266,'Q11')
        call collocate_MPI_averages_N(aQ22f,avg_out,267,'Q22')
        call collocate_MPI_averages_N(aQ33f,avg_out,268,'Q33')
        call collocate_MPI_averages_N(awwu,avg_out,269,'wwu')
        call collocate_MPI_averages_N(autzz,avg_out,270,'utzz')
        call collocate_MPI_averages_N(at11wx,avg_out,271,'t11wx')
        call collocate_MPI_averages_N(at12wy,avg_out,272,'t12wy')
        call collocate_MPI_averages_N(at13wz,avg_out,273,'t13wz')
        call collocate_MPI_averages_N(at31ux,avg_out,274,'t31ux')
        call collocate_MPI_averages_N(at32uy,avg_out,275,'t32uy')
        call collocate_MPI_averages_N(at33uz,avg_out,276,'t33uz')
        call collocate_MPI_averages_N(as11,avg_out,277,'s11')
        call collocate_MPI_averages_N(as12,avg_out,278,'s12')
        call collocate_MPI_averages_N(as13,avg_out,279,'s13')
        call collocate_MPI_averages_N(as22,avg_out,280,'s22')
        call collocate_MPI_averages_N(as23,avg_out,281,'s23')
        call collocate_MPI_averages_N(as33,avg_out,282,'s33')
        call collocate_MPI_averages_N(adudx,avg_out,283,'dudx')
        call collocate_MPI_averages_N(adudy,avg_out,284,'dudy')
        call collocate_MPI_averages_N(advdx,avg_out,285,'dvdx')
        call collocate_MPI_averages_N(advdy,avg_out,286,'dvdy')
        call collocate_MPI_averages_N(adwdx,avg_out,287,'dwdx')
        call collocate_MPI_averages_N(adwdy,avg_out,288,'dwdy')
        call collocate_MPI_averages_N(adwdz,avg_out,289,'dwdz')
        call collocate_MPI_averages_N(adivu,avg_out,290,'divu')
        call collocate_MPI_averages_N(apdivu,avg_out,291,'pdivu')
        call collocate_MPI_averages_N(apfdivu,avg_out,292,'pfdivu')
        call collocate_MPI_averages_N(awwv,avg_out,293,'wwv')
        call collocate_MPI_averages_N(auf,avg_out,294,'zzuf')
        call collocate_MPI_averages_N(avf,avg_out,295,'zzvf')
        call collocate_MPI_averages_N(awf,avg_out,296,'zzwf')
        deallocate(avg_out)
!Zero out the outputted averages !!
   au=0._rprec;av=0._rprec;aw=0._rprec;ap=0._rprec;p2=0._rprec;u2=0._rprec;v2=0._rprec;w2=0._rprec;
   auw=0._rprec;avw=0._rprec;auv=0._rprec;atke=0._rprec;au2w=0._rprec;av2w=0._rprec;
   awe=0._rprec;adissip=0._rprec;u3=0._rprec;v3=0._rprec;w3=0._rprec;
   adudz=0._rprec;advdz=0._rprec;aCs_Ssim=0._rprec;asig=0._rprec;
   asigQ=0._rprec;aCs=0._rprec;aQ11f=0._rprec;aQ22f=0._rprec;aQ33f=0._rprec;
   aL11f=0._rprec;aL22f=0._rprec;aL33f=0._rprec;atxx=0._rprec;atxy=0._rprec;atxz=0._rprec;atyy=0._rprec;atyz=0._rprec;
   atzz=0._rprec;autxx=0._rprec;autxy=0._rprec;autxz=0._rprec;
   avtxy=0._rprec;avtyy=0._rprec;avtyz=0._rprec;
   awtxz=0._rprec;awtzz=0._rprec;awtyz=0._rprec;
   aup=0._rprec;avp=0._rprec;awp=0._rprec;ausig=0._rprec;avsig=0._rprec;awsig=0._rprec;asigux=0._rprec;
   asiguy=0._rprec;asiguz=0._rprec;asigvx=0._rprec;asigvy=0._rprec;asigvz=0._rprec;asigwx=0._rprec;
   asigwy=0._rprec;asigwz=0._rprec;at11ux=0._rprec;at12uy=0._rprec;at13uz=0._rprec;at21vx=0._rprec;
   at22vy=0._rprec;at23vz=0._rprec;at31wx=0._rprec;at32wy=0._rprec;at33wz=0._rprec;
   apux=0._rprec;apuy=0._rprec;apuz=0._rprec;apvx=0._rprec;apvy=0._rprec;apvz=0._rprec;
   apwx=0._rprec;apwy=0._rprec;
   apwz=0._rprec;apux2=0._rprec;apuy2=0._rprec;apuz2=0._rprec;apvx2=0._rprec;apvy2=0._rprec;apvz2=0._rprec;
   apwx2=0._rprec;apwy2=0._rprec;apwz2=0._rprec;
   at11s11=0._rprec;at12s12=0._rprec;at13s13=0._rprec;at22s22=0._rprec;
   at23s23=0._rprec;at33s33=0._rprec;
   aupf=0._rprec;awpf=0._rprec;avpf=0._rprec;
   aCs_Ssim=0._rprec;
   ausigQ=0._rprec;avsigQ=0._rprec;awsigQ=0._rprec;asigQux=0._rprec;
   asigQuy=0._rprec;asigQuz=0._rprec;asigQvx=0._rprec;asigQvy=0._rprec;asigQvz=0._rprec;asigQwx=0._rprec;
   asigQwy=0._rprec;asigQwz=0._rprec;ap_c=0._rprec;ap_c2=0._rprec;awwu=0._rprec;awwv=0._rprec;
   autzz=0._rprec;at11wx=0._rprec;at12wy=0._rprec;at13wz=0._rprec;at31ux=0._rprec;
   at32uy=0._rprec;at33uz=0._rprec;
   as11=0._rprec;as12=0._rprec;as13=0._rprec;as22=0._rprec;as23=0._rprec;as33=0._rprec;
   adudx=0._rprec;adudy=0._rprec;advdx=0._rprec;advdy=0._rprec;
   adwdx=0._rprec;adwdy=0._rprec;adwdz=0._rprec;adivu=0._rprec;apdivu=0._rprec;apfdivu=0._rprec;
   auf=0._rprec;avf=0._rprec;awf=0._rprec;

   end if
end subroutine MM_budget_slice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--This Subroutine Give the output for whole domain (MM)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine MM_XYZ_Out
!use sim_param,only:path,u,v,w,dudz,dvdz,txx,txz,tyy,tyz,tzz,p,theta
!use param,only:dz,p_count,c_count,jt
!use sgsmodule,only:Cs_opt2,Cs_Ssim,Beta_avg,Betaclip_avg
!implicit none
!integer::i,j,k
!real(kind=rprec),dimension(nx,ny,nz-1),save::ap,au,av,aw,p2,u2,v2,w2,auw,avw,acs,atheta
!real(kind=rprec),dimension(nx,ny,nz-1),save::adudz,advdz,aCs_Ssim,abeta_sgs,abetaclip_sgs
!real(kind=rprec),dimension(nx,ny,nz-1),save::atxx,atxz,atyy,atyz,atzz
!real(kind=rprec),dimension(nx,ny,nz-1),save::u3,v3,w3
!real(kind=rprec),dimension(:,:,:),allocatable::avg_out
!real(kind=rprec)::fr,arg1, arg2
!character (len=256) :: local_filename
!
!fr=(1._rprec/real(p_count,kind=rprec))*real(c_count,kind=rprec)
!do k=1,Nz-1
!do i=1,Nx
!do j=1,Ny
!      au(i,j,k)=au(i,j,k)+fr*u(i,j,k)
!      av(i,j,k)=av(i,j,k)+fr*v(i,j,k)
!      aw(i,j,k)=aw(i,j,k)+fr*w(i,j,k)
!      atheta(i,j,k)=atheta(i,j,k)+fr*theta(i,j,k)
!      ap(i,j,k)=ap(i,j,k)+fr*p(i,j,k)
!      atxx(i,j,k)=atxx(i,j,k)+fr*txx(i,j,k)
!      atxz(i,j,k)=atxz(i,j,k)+fr*txz(i,j,k)
!      atyy(i,j,k)=atyy(i,j,k)+fr*tyy(i,j,k)
!      atyz(i,j,k)=atyz(i,j,k)+fr*tyz(i,j,k)
!      atzz(i,j,k)=atzz(i,j,k)+fr*tzz(i,j,k)
!      adudz(i,j,k)=adudz(i,j,k)+fr*dudz(i,j,k)
!      advdz(i,j,k)=advdz(i,j,k)+fr*dvdz(i,j,k)
!      u2(i,j,k)=u2(i,j,k)+fr*u(i,j,k)*u(i,j,k)
!      v2(i,j,k)=v2(i,j,k)+fr*v(i,j,k)*v(i,j,k)
!      w2(i,j,k)=w2(i,j,k)+fr*w(i,j,k)*w(i,j,k)
!      p2(i,j,k)=p2(i,j,k)+fr*p(i,j,k)*p(i,j,k)
!      aCs(i,j,k)=sqrt(Cs_opt2(i,j,k))
!      aCs_Ssim(i,j,k)=sqrt(Cs_Ssim(i,j,k))
!      if((k .eq. 1) .AND. ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0))) then
!         arg1=0._rprec
!         arg2=0._rprec
!      else
!         arg1=(u(i,j,k)+u(i,j,k-1))/2.
!         arg2=(v(i,j,k)+v(i,j,k-1))/2.
!      end if
!      auw(i,j,k)=auw(i,j,k)+fr*w(i,j,k)*arg1
!      avw(i,j,k)=avw(i,j,k)+fr*w(i,j,k)*arg2
!      u3(i,j,k)=u3(i,j,k)+fr*u(i,j,k)*u(i,j,k)*u(i,j,k)
!      v3(i,j,k)=v3(i,j,k)+fr*v(i,j,k)*v(i,j,k)*v(i,j,k)
!      w3(i,j,k)=w3(i,j,k)+fr*w(i,j,k)*w(i,j,k)*w(i,j,k)
!      abeta_sgs(i,j,k)=abeta_sgs(i,j,k)+fr*Beta_avg(k)
!      abetaclip_sgs(i,j,k)=abetaclip_sgs(i,j,k)+fr*Betaclip_avg(k)
!end do
!end do
!end do
!
!if (mod(jt,p_count)==0) then
!        allocate(avg_out(1:nx,1:ny,1:(nz_tot-1)));
!        call collocate_MPI_averages_SHH(au,avg_out,720,'u_all')
!        call collocate_MPI_averages_SHH(av,avg_out,721,'v_all')
!        call collocate_MPI_averages_SHH(aw,avg_out,722,'w_all')
!        call collocate_MPI_averages_SHH(ap,avg_out,723,'p_all')
!        call collocate_MPI_averages_SHH(u2,avg_out,724,'u2_all')
!        call collocate_MPI_averages_SHH(v2,avg_out,725,'v2_all')
!        call collocate_MPI_averages_SHH(w2,avg_out,726,'w2_all')
!        call collocate_MPI_averages_SHH(p2,avg_out,732,'p2_all')
!        call collocate_MPI_averages_SHH(atxx,avg_out,727,'txx_all')
!        call collocate_MPI_averages_SHH(atxz,avg_out,728,'txz_all')
!        call collocate_MPI_averages_SHH(atyy,avg_out,729,'tyy_all')
!        call collocate_MPI_averages_SHH(atyz,avg_out,730,'tyz_all')
!        call collocate_MPI_averages_SHH(atzz,avg_out,731,'tzz_all')
!        call collocate_MPI_averages_SHH(auw,avg_out,733,'uw_all')
!        call collocate_MPI_averages_SHH(avw,avg_out,734,'vw_all')
!        call collocate_MPI_averages_SHH(aCs,avg_out,735,'Cs_all')
!        call collocate_MPI_averages_SHH(adudz,avg_out,736,'dudz_all')
!        call collocate_MPI_averages_SHH(advdz,avg_out,737,'dvdz_all')
!        call collocate_MPI_averages_SHH(aCs_Ssim,avg_out,738,'Cs_Ssim_all')
!        call collocate_MPI_averages_SHH(abeta_sgs,avg_out,739,'beta_sgs_all')
!        call collocate_MPI_averages_SHH(abetaclip_sgs,avg_out,740,'betaclip_sgs_all');
!        call collocate_MPI_averages_SHH(u3,avg_out,741,'u3_all')
!        call collocate_MPI_averages_SHH(v3,avg_out,742,'v3_all')
!        call collocate_MPI_averages_SHH(w3,avg_out,743,'w3_all');
!        call collocate_MPI_averages_SHH(atheta,avg_out,744,'theta_all')
!
!        deallocate(avg_out)
!
!!VK Zero out the outputted averages !!
!        au=0._rprec;av=0._rprec;aw=0._rprec;ap=0._rprec;u2=0._rprec;v2=0._rprec
!        w2=0._rprec;atxx=0._rprec;atxz=0._rprec;atyy=0._rprec;atyz=0._rprec
!        atzz=0._rprec;p2=0._rprec;auw=0._rprec;avw=0._rprec;aCs=0._rprec;atheta=0._rprec
!        adudz=0._rprec;advdz=0._rprec;aCs_Ssim=0._rprec;abeta_sgs=0._rprec
!        abetaclip_sgs=0._rprec;u3=0._rprec;v3=0._rprec;w3=0._rprec;
!end if
!5168     format(1400(E14.5))
!end subroutine MM_XYZ_Out
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--The following subroutine does the collocation of the MPI arrays for
! SHHOutput Subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine collocate_MPI_averages_SHH(avg_var_proc,avg_var_tot_domain,file_ind,filename_str)
!subroutine collocate_MPI_averages(avg_var_proc,avg_var_tot_domain,file_ind)
use param
$if ($MPI)
  integer :: recvcounts(nproc)
  integer :: displs(nproc)
$endif
integer :: ind1,ind2,ind3,file_ind
character (*),intent(in) :: filename_str
character (len=256) :: local_filename
real(kind=rprec),dimension(nx,ny,nz-1)::avg_var_proc
real(kind=rprec),dimension(nx,ny,nz_tot-1)::avg_var_tot_domain

local_filename=path//'output/aver_'//trim(filename_str)//'.out'

avg_var_tot_domain=0._rprec
$if ($MPI)
  recvcounts = size (avg_var_proc)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (avg_var_proc(1,1,1), size (avg_var_proc), MPI_RPREC,       &
                    avg_var_tot_domain(1,1,1), recvcounts, displs, MPI_RPREC, &
                    rank_of_coord(0), comm, ierr)
$else
  avg_var_tot_domain=avg_var_proc
$endif

  if((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
        open(file_ind,file=trim(local_filename),status="unknown",position="append")
           do ind3=1,nz_tot-1
           do ind2=1,ny
            write(file_ind,5168)jt*dt,(avg_var_tot_domain(ind1,ind2,ind3),ind1=1,nx)
           end do
           end do
        close(file_ind)
  end if
5168     format(1400(E14.5))
end subroutine collocate_MPI_averages_SHH













!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--assumes lun is open and positioned correctly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine checkpoint (lun)

use param,only:nz,S_FLAG
use sim_param,only:u,v,w,RHSx,RHSy,RHSz,theta
use sgsmodule,only:Cs_opt2,F_LM,F_MM,F_QN,F_NN,G_LM,G_MM,G_QN,G_NN,Pr_t
use scalars_module,only:RHS_T,sgs_t3,psi_m

implicit none
integer,intent(in)::lun

if (S_FLAG) then ! WITH SCALARS
   write (lun) u(:,:,1:nz),v(:,:,1:nz),w(:,:,1:nz),theta(:,:,1:nz),   &
               RHSx(:,:,1:nz),RHSy(:,:,1:nz),RHSz(:,:,1:nz),          &
               RHS_T(:,:,1:nz),sgs_t3(:,:,1),psi_m,Cs_opt2,F_LM,F_MM, &
               F_QN,F_NN,G_LM,G_MM,G_QN,G_NN,Pr_t(:,:,1:nz)
else ! No SCALARS
   write (lun) u(:,:,1:nz),v(:,:,1:nz),w(:,:,1:nz),          &
               RHSx(:,:,1:nz),RHSy(:,:,1:nz),RHSz(:,:,1:nz), &
               Cs_opt2,F_LM,F_MM,F_QN,F_NN
end if

end subroutine checkpoint

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--lun_opt gives option for writing to different unit, and is used by inflow_write
!--assumes the unit to write to (lun_default or lun_opt is already
!  open for sequential unformatted writing
!--this routine also closes the unit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_final(jt, lun_opt)
implicit none

integer,intent(in)::jt
integer, intent (in), optional :: lun_opt  !--if present, write to unit lun
integer, parameter :: lun_default = 11
integer::jx,jy,jz
integer :: lun
logical :: opn

if (present (lun_opt)) then
  lun = lun_opt
else
  lun = lun_default
end if

inquire (unit=lun, opened=opn)

if (.not. opn) then
  write (*, *) 'output_final: lun=', lun, ' is not open'
  stop
end if

rewind (lun)

call checkpoint (lun)

close (lun)

if ((cumulative_time) .and. (lun == lun_default)) then
  !--only do this for true final output, not intermediate recording
  open (1, file=fcumulative_time)
  write (1, *) jt_total
  close (1)
end if

end subroutine output_final

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine io_lambda2_out
use sim_param,only:path
implicit none
character(len=24)::fname
call lambda2()
write(fname,'(A13,i6.6,A4)')path//'output/lam-',jt_total,'.out'
open(1,file=fname,form='unformatted')
write(1)nx,ny,nz
write(1)real(lam2)
close(1)
end subroutine io_lambda2_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine lambda2()
use types,only:rprec
use sim_param,only:u,v,w,&
     dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
use param,only:dx,dy,dz
implicit none
real(kind=rprec)::S11,S22,S33,S12,O12,S13,O13,S23,O23,&
     ux,uy,uz,vx,vy,vz,wx,wy,wz
integer::jx,jy,jz
! following used for eispack call...
integer::neis,nmeis,matzeis,ierreis,iv1eis(3)
double precision,dimension(3,3)::aeis,zeis
double precision,dimension(3)::wreis,wieis,fv1eis
double precision::ave

! assignments for eispack call
neis=3
nmeis=3
matzeis=0
ierreis=0
lam2=0._rprec

! at level z=dz/2.  (Where dudz and dvdz are on UVP nodes)
jz=1
do jy=1,ny
do jx=1,nx              
   ux=dudx(jx,jy,1)  ! uvp-node
   uy=dudy(jx,jy,1)  ! uvp-node
   uz=dudz(jx,jy,1)  ! uvp-node
   vx=dvdx(jx,jy,1)  ! uvp-node
   vy=dvdy(jx,jy,1)  ! uvp-node
   vz=dvdz(jx,jy,1)  ! uvp-node 
! special case
   wx=0.5_rprec*(dwdx(jx,jy,1)+dwdx(jx,jy,2))  ! uvp-node
   wy=0.5_rprec*(dwdy(jx,jy,1)+dwdy(jx,jy,2))  ! uvp-node
   wz=dwdz(jx,jy,1)  ! uvp-node
   S11=ux          ! uvp-node
   S12=0.5_rprec*(uy+vx) ! uvp-node
! taken care of with wall stress routine
   S13=0.5_rprec*(uz+wx) ! uvp
   O12=0.5_rprec*(uy-vx) ! w-node
   O13=0.5_rprec*(uz-wx) ! w-node
   S22=vy          ! uvp-node
! taken care of with wall stress routine 
   S23=0.5_rprec*(vz+wy) ! uvp
   O23=0.5_rprec*(vz-wy) ! w-node
   S33=wz          ! uvp-node
   aeis(1,1)=s11*s11+s12*s12+s13*s13-O12*O12-O13*O13
   aeis(1,2)=s11*s12+s12*s22+s13*s23-O13*O23
   aeis(1,3)=s11*s13+s12*s23+s13*s33+O12*O23
   aeis(2,2)=s12*s12+s22*s22+s23*s23-O12*O12-O23*O23
   aeis(2,3)=s12*s13+s22*s23+s23*s33-O12*O13
   aeis(3,3)=s13*s13+s23*s23+s33*s33-O13*O13-O23*O23
   aeis(2,1)=aeis(1,2)
   aeis(3,1)=aeis(1,3)
   aeis(3,2)=aeis(2,3)
  write (*, *) 'rg temporarily removed, sorry'; stop
   if(wreis(1).ge.wreis(2).and.wreis(1).le.wreis(3)) then
      lam2(jx,jy,jz)=real(wreis(1),kind=rprec)
   elseif(wreis(2).ge.wreis(1).and.wreis(2).le.wreis(3)) then
      lam2(jx,jy,jz)=real(wreis(2),kind=rprec)
   else
      lam2(jx,jy,jz)=real(wreis(3),kind=rprec)
   endif
end do
end do
! calculate derivatives/strain on w-nodes
do jz=2,nz-1  
do jy=1,ny
do jx=1,nx              
   ux=0.5_rprec*(dudx(jx,jy,jz) + dudx(jx,jy,jz-1))  ! w-node
   uy=0.5_rprec*(dudy(jx,jy,jz) + dudy(jx,jy,jz-1))  ! w-node
   uz=dudz(jx,jy,jz)  ! w-node
   vx=0.5_rprec*(dvdx(jx,jy,jz) + dvdx(jx,jy,jz-1))  ! w-node
   vy=0.5_rprec*(dvdy(jx,jy,jz) + dvdy(jx,jy,jz-1))  ! w-node
   vz=dvdz(jx,jy,jz)  ! w-node
   wx=dwdx(jx,jy,jz)  ! w-node
   wy=dwdy(jx,jy,jz)  ! w-node
   wz=0.5_rprec*(dwdz(jx,jy,jz) + dwdz(jx,jy,jz-1))  ! w-node
   S11=ux          ! w-node
   S12=0.5_rprec*(uy+vx) ! w-node
   S13=0.5_rprec*(uz+wx) ! w-node
   O12=0.5_rprec*(uy-vx) ! w-node
   O13=0.5_rprec*(uz-wx) ! w-node
   S22=vy          ! w-node
   S23=0.5_rprec*(vz+wy) ! w-node
   O23=0.5_rprec*(vz-wy) ! w-node
   S33=wz          ! w-node
   aeis(1,1)=s11*s11+s12*s12+s13*s13-O12*O12-O13*O13
   aeis(1,2)=s11*s12+s12*s22+s13*s23-O13*O23
   aeis(1,3)=s11*s13+s12*s23+s13*s33+O12*O23
   aeis(2,2)=s12*s12+s22*s22+s23*s23-O12*O12-O23*O23
   aeis(2,3)=s12*s13+s22*s23+s23*s33-O12*O13
   aeis(3,3)=s13*s13+s23*s23+s33*s33-O13*O13-O23*O23
   aeis(2,1)=aeis(1,2)
   aeis(3,1)=aeis(1,3)
   aeis(3,2)=aeis(2,3)
   if(wreis(1).ge.wreis(2).and.wreis(1).le.wreis(3)) then
      lam2(jx,jy,jz)=real(wreis(1),kind=rprec)
   elseif(wreis(2).ge.wreis(1).and.wreis(2).le.wreis(3)) then
      lam2(jx,jy,jz)=real(wreis(2),kind=rprec)
   else
      lam2(jx,jy,jz)=real(wreis(3),kind=rprec)
   endif
end do
end do
end do

print*,'minmax',minval(lam2),maxval(lam2)
end subroutine lambda2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine io_mean_out
implicit none
write(51)real(mean_u),real(mean_u2),real(mean_v),real(mean_v2),&
     real(mean_w),real(mean_w2)
!--the nz/4*3 stuff has to go
mean_u(jx_pls:jx_ple,jy_pls:jy_ple,1:nz/4*3)=0._rprec
mean_u2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz/4*3)=0._rprec
mean_v(jx_pls:jx_ple,jy_pls:jy_ple,1:nz/4*3)=0._rprec
mean_v2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz/4*3)=0._rprec
mean_w(jx_pls:jx_ple,jy_pls:jy_ple,1:nz/4*3)=0._rprec
mean_w2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz/4*3)=0._rprec
end subroutine io_mean_out

subroutine calculate_mean
use sim_param,only:u,v,w
use sgsmodule,only:Cs_opt2,Cs_opt2_avg
implicit none
Cs_opt2_avg(:,:,:)=Cs_opt2_avg(:,:,:)+Cs_opt2(:,:,:)/nwrite
!TS
mean_u(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)=&
     mean_u(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)+&
     u(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)/nwrite
mean_u2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)=&
     mean_u2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)+&
     u(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)**2/nwrite
mean_v(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)=&
     mean_v(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)+&
     v(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)/nwrite
mean_v2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)=&
     mean_v2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)+&
     v(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)**2/nwrite
mean_w(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)=&
     mean_w(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)+&
     w(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)/nwrite
mean_w2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)=&
     mean_w2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)+&
     w(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)**2/nwrite
end subroutine calculate_mean

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine timeseries_spec
use sim_param,only:u,v,w,theta
implicit none
integer::jx,jy,jz,i
if(mod(jt_total,time_spec)==0.and.jt_total.gt.2000)then
jx=NX/8
jy=NY/2+1
jz=NZ/2
endif
end subroutine timeseries_spec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine post_spec(jt_local)
use sim_param,only:path,u,v,w,theta
use param
use fft
implicit none
real(kind=rprec),dimension(nx/2,nz)::spectra_u,spectra_v,spectra_w,&
     spectra_theta
real(kind=rprec),dimension(4,nx/2,nz-1)::spectra_uvwT
real(kind=rprec),dimension(4,nx/2,nz_tot-1)::spectra_uvwT_tot
integer,intent(in)::jt_local
integer::k,jz,z
character(len=64)::fname1,fname2,fname3,fname4
$if ($MPI)
  $define $lbz 0
  integer :: recvcounts(nproc)
  integer :: displs(nproc)
$else
  $define $lbz 1
$endif

write(fname1,'(a,a)') path//'output/spec_x','.dat'
open(82,file=fname1,form='formatted')
do jz=1,nz-1
   z=(jz-0.5_rprec)*dz*z_i
   write(82,*) (real(kx(k,1)/z_i*z),k=1,nx/2)
   call spectrum(u(:, :, jz), spectra_uvwT(1,:,jz))
   call spectrum(v(:, :, jz), spectra_uvwT(2,:,jz))
   call spectrum(w(:, :, jz), spectra_uvwT(3,:,jz))
   call spectrum(theta(:, :, jz), spectra_uvwT(4,:,jz))
enddo
   close(82)
$if ($MPI)
  recvcounts = size (spectra_uvwT)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (spectra_uvwT(1, 1,1), size (spectra_uvwT), MPI_RPREC,&
                    spectra_uvwT_tot(1, 1, 1), recvcounts, displs,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
$else
  spectra_uvwT_tot=spectra_uvwT
$endif

if((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then  
   write(fname1,'(A,i6.6,A)')path//'output/spec_uvwT_',jt_local,'.bin'
   open(83,file=fname1,form='unformatted')
   write(83) real(spectra_uvwT_tot(:,1:nx/2,:))
   close(83)
end if

end subroutine post_spec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine spectrum(u, spec)
use fft
implicit none      
real(kind=rprec),dimension(ld,ny),intent(in)::u
real(kind=rprec),dimension(nx/2),intent(out)::spec  !--assumes Nyquist is 0

integer::jy,jz,k
real(kind=rprec),dimension(nx)::vel_r,vel_c

integer*8, save :: plan
logical, save :: init = .false.

if (.not. init) then
  call rfftw_f77_create_plan(plan,nx,FFTW_REAL_TO_COMPLEX,FFTW_MEASURE)
  init = .true.
end if

! initialize
spec(:)=0._rprec
do jy=1,ny
   vel_r(:)= u(1:nx,jy)/real(nx,kind=rprec)
! check this normaliztion-part of forward; call the fft
   call rfftw_f77_one(plan,vel_r,vel_c)
! compute magnitudes the 0.5 is the 1/2, all others are taken care of! (except maybe Nyquist)
   spec(1)=spec(1)+0.5*vel_c(1)*vel_c(1)
   do k=2,nx/2
      spec(k)=spec(k)+vel_c(k)*vel_c(k)+vel_c(nx+2-k)*vel_c(nx+2-k)
   end do

   !--assume Nyquist is 0
   !spec(nx/2+1)=spec(nx/2+1)+vel_c(nx/2+1)*vel_c(nx/2+1)
end do
spec(:)=spec(:)/real(Ny,kind=rprec) ! for average over Ny
end subroutine spectrum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine avg_stats ()
use param
use sim_param, only : u, v, w, txz
use fft, only : kx
implicit none

!--choose naming convention that does not conflict with qpost
character (*), parameter :: fubar_avg = 'output/ubar-avg_stats.dat'
character (*), parameter :: fupr2bar_avg = 'output/upr2bar-avg_stats.dat'
character (*), parameter :: fstressbar_avg = 'output/stressbar-avg_stats.dat'
character (*), parameter :: fEozbar_avg = 'output/Eozbar-avg_stats.dat'

integer, parameter :: hdr_len = 256
logical, parameter :: DEBUG = .false.
character (hdr_len) :: Eozbar_hdr

$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif
$if ($MPI)
  integer :: recvcounts(nproc)
  integer :: displs(nproc)
$endif
integer, save :: n_ubar_avg
integer, save :: n_upr2bar_avg
integer, save :: n_stressbar_avg
integer, save :: n_Eozbar_avg
integer :: jz

logical, save :: init = .false.

real (rprec) :: z
real (rprec) :: zu(1, nz_tot-1)
real (rprec) :: kz_z(2, nx/2)
real (rprec), save :: ubar_avg(1, nz_tot-1)      !--1 is <u>
real (rprec), save :: upr2bar_avg(3, nz_tot-1)   !--<u'^2>, <v'^2>, <w'^2>
real (rprec), save :: stressbar_avg(3, nz_tot-1) !--1 is <u'w'>, 2 is <txz>, 3 is <u'w'> + <txz>
real (rprec), save :: Eozbar_avg(1, nx/2, nz_tot-1)  !--E11(k1,z)/z
!--tot is a temp for current stats at nz_tot size
real (rprec), save :: ubar_tot(1, nz_tot-1)      !--1 is <u>
real (rprec), save :: upr2bar_tot(3, nz_tot-1)   !--<u'^2>, <v'^2>, <w'^2>
real (rprec), save :: stressbar_tot(3, nz_tot-1) !--1 is <u'w'>, 2 is <txz>, 3 is <u'w'> + <txz>
real (rprec), save :: Eozbar_tot(1, nx/2, nz_tot-1)  !--E11(k1,z)/z
real (rprec) :: upr(nx, ny), vpr(nx, ny), wpr(nx, ny)
real (rprec) :: ubar(nz-1), vbar(nz-1), wbar(nz-1)
real (rprec) :: upr2bar(3, nz-1)
real (rprec) :: stressbar(3, nz-1)
real (rprec) :: Eozbar(nx/2, nz-1)
!---------------------------------------------------------------------

!--check whether or not to actually do anything
!--motivation for doing this way is that it cleans up interface in main
if (modulo (jt, n_avg_stats) /= 0) goto 001  !--do nothing, exit cleanly

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then  
  if (.not. init) then  !--initialization

    call init_avg (fubar_avg, 1, ubar_avg, n_ubar_avg)
    call init_avg (fupr2bar_avg, 1, upr2bar_avg, n_upr2bar_avg)
    call init_avg (fstressbar_avg, 1, stressbar_avg, n_stressbar_avg) 
    do jz = 1, nz-2
      call init_avg (fEozbar_avg, 2, Eozbar_avg(:, :, jz), n_Eozbar_avg,  &
                     leaveopn='yes')
    end do
    call init_avg (fEozbar_avg, 2, Eozbar_avg(:, :, nz-1), n_Eozbar_avg)

    init = .true.

  end if
end if

!--calculation of current stats
do jz = $lbz, nz-1

  ubar(jz) = sum (u(1:nx, 1:ny, jz)) / (nx * ny)
  vbar(jz) = sum (v(1:nx, 1:ny, jz)) / (nx * ny)
  wbar(jz) = sum (w(1:nx, 1:ny, jz)) / (nx * ny)

  if ( ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) .and.  &
       (jz == 1) ) then
    upr = 0._rprec
    vpr = 0._rprec
    wpr = 0._rprec
  else
    !--see qpost for u/w-node interpolation
    !--convention will be to put up, vp, wp on w-nodes
    upr = 0.5_rprec * (u(1:nx, 1:ny, jz) - ubar(jz) +  &
                       u(1:nx, 1:ny, jz-1) - ubar(jz-1))
    vpr = 0.5_rprec * (v(1:nx, 1:ny, jz) - vbar(jz) +  &
                       v(1:nx, 1:ny, jz-1) - vbar(jz-1))
    wpr = w(1:nx, 1:ny, jz) - wbar(jz)
  end if
 
  upr2bar(1, jz) = sum (upr**2) / (nx * ny)
  upr2bar(2, jz) = sum (vpr**2) / (nx * ny)
  upr2bar(3, jz) = sum (wpr**2) / (nx * ny)

  stressbar(1, jz) = sum (upr * wpr) / (nx * ny) 
  stressbar(2, jz) = sum (txz(1:nx, 1:ny, jz)) / (nx * ny)
  stressbar(3, jz) = sum (stressbar(1:2, jz))

  !--energy spectra
  call spectrum (u(:, :, jz), Eozbar(:, jz))  !--not /z yet
  z = (jz - 0.5_rprec) * dz
  Eozbar(:, jz) = Eozbar(:, jz) / z

end do

!--collect current stats into nz_tot sized arrays
$if ($MPI)

  if (DEBUG) then
    write (*, *) coord, ': ubar(1) = ', ubar(1)
  end if

  recvcounts = size (ubar)
  displs = coord_of_rank * recvcounts 
  call mpi_gatherv (ubar(1), size (ubar), MPI_RPREC,                &
                    ubar_tot(1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)

  recvcounts = size (upr2bar)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (upr2bar(1, 1), size (upr2bar), MPI_RPREC,          &
                    upr2bar_tot(1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)  

  recvcounts = size (stressbar)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (stressbar(1, 1), size (stressbar), MPI_RPREC,        &
                    stressbar_tot(1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)

  recvcounts = size (Eozbar)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (Eozbar(1, 1), size (Eozbar), MPI_RPREC,              &
                    Eozbar_tot(1, 1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)

$else

  ubar_tot(1, :) = ubar
  upr2bar_tot = upr2bar
  stressbar_tot = stressbar
  Eozbar_tot(1, :, :) = Eozbar

$endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  !--calculation of cumulative average stats
  ubar_avg = (n_ubar_avg * ubar_avg + ubar_tot) / (n_ubar_avg + 1)
  n_ubar_avg = n_ubar_avg + 1

  upr2bar_avg = (n_upr2bar_avg * upr2bar_avg + upr2bar_tot) /  &
                (n_upr2bar_avg + 1)
  n_upr2bar_avg = n_upr2bar_avg + 1

  stressbar_avg = (n_stressbar_avg * stressbar_avg + stressbar_tot) /  &
                  (n_stressbar_avg + 1)
  n_stressbar_avg = n_stressbar_avg + 1

  Eozbar_avg = (n_Eozbar_avg * Eozbar_avg + Eozbar_tot) / (n_Eozbar_avg + 1)
  n_Eozbar_avg = n_Eozbar_avg + 1

  !--prepare list of z-coordinates
  forall (jz=1:nz_tot-1) zu(1, jz) = (jz - 0.5_rprec) * dz
  !--prepare  header, optional

  !--write out to file
  call write_avg (fubar_avg, n_ubar_avg, zu, ubar_avg)
  call write_avg (fupr2bar_avg, n_upr2bar_avg, zu, upr2bar_avg)
  call write_avg (fstressbar_avg, n_stressbar_avg, zu, stressbar_avg)

  !--this is a bit awkward: maybe just have another routine to do it right
  Eozbar_hdr = 'zone' !--this is for tecplot... 
  kz_z(1, :) = kx(1:nx/2, 1) * zu(1, 1)
  kz_z(2, :) = zu(1, 1)
  call write_avg (fEozbar_avg, n_Eozbar_avg, kz_z, Eozbar_avg(:, :, 1),  &
                  hdr=Eozbar_hdr) 

  do jz = 2, nz_tot - 1
    kz_z(1, :) = kx(1:nx/2, 1) * zu(1, jz)
    kz_z(2, :) = zu(1, jz)
    call write_avg (fEozbar_avg, n_Eozbar_avg, kz_z, Eozbar_avg(:, :, jz),  &
                    hdr=Eozbar_hdr, position='append') 
  end do

end if

001 continue  !--exit cleanly

end subroutine avg_stats

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine init_avg (file_avg, n_ccol, a_avg, n_avg, leaveopn)
implicit none

character (*), intent (in) :: file_avg
integer, intent (in) :: n_ccol  !--num. coord columns: x, y, etc.
real (rprec), intent (out) :: a_avg(:, :)
integer, intent (out) :: n_avg
character (*), optional, intent (in) :: leaveopn
character (128) :: buff
logical :: exst, opn
integer :: j
real (rprec) :: z(n_ccol)

!---------------------------------------------------------------------
inquire (file=file_avg, exist=exst, opened=opn)

if (exst) then

  if (.not. opn) then
    open (1, file=file_avg)
    read (1, '(a)') buff

    if (buff(1:1) == '#') then
      read (buff(2:), *) n_avg
    else
      write (*, *) 'avg_stats: error'
      write (*, *) trim (file_avg), ' does not have expected format on line 1'
      stop  !--need to replace all stops with nice mpi exits
    end if
  end if

  !--skip data header lines here
  do
    read (1, '(a)') buff
    if (trim (buff) == trim (end_hdr_avg)) exit
  end do

  do j = 1, size (a_avg, 2)
    read (1, *) z, a_avg(:, j)  !--z is just placeholder here
  end do

  if (present (leaveopn)) then
    if (leaveopn /= 'yes') close (1)  !--case sensitive here
  else
    close (1)
  end if

else

  n_avg = 0
  a_avg = 0._rprec

end if

end subroutine init_avg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_avg (file_avg, n_avg, x, a_avg, hdr, position)
implicit none

character (*), intent (in) :: file_avg
integer, intent (in) :: n_avg
real (rprec), intent (in) :: x(:, :)  !--coord columns for x, y, etc
real (rprec), intent (in) :: a_avg(:, :)

character (*), optional, intent (in) :: hdr
character (*), optional, intent (in) :: position

character (64) :: r_fmt, fmt
character (32) :: posn

integer :: j

!---------------------------------------------------------------------

!--check sizes compatible
if (size (x, 2) /= size (a_avg, 2)) then
  write (*, *) 'write_avg: error with sizes of x, a_avg'
  stop
end if

if (present (position)) then
  posn = position
else
  posn = 'rewind'
end if

open (1, file=file_avg, position=posn)

if (trim (posn) /= 'append') then  !--case sensitive
  write (1, '(a,i0)') '# ', n_avg  !--not needed when appending
end if

if (present (hdr)) then
  !--write data header, if present
  write (1, '(a)') trim (hdr)
end if

!--write something to indicate end of header, always do this
write (1, '(a)') end_hdr_avg

!--define output format
write (r_fmt, '(2(a,i0))') 'es', precision (1._rprec) + 7,  &
                           '.', precision (1._rprec)
write (fmt, '(a,i0,3a)') '(', size (x, 1) + size (a_avg, 1),  &
                         '(1x,', trim (r_fmt), '))'

!--write to file
do j = 1, size (a_avg, 2)
  write (1, fmt) x(:, j), a_avg(:, j)
end do

close (1)

end subroutine write_avg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine inflow_write ()
use param, only : jt_total, model, jt_start_write, buff_end,  &
                  read_inflow_file, write_inflow_file
use sgsmodule, only : F_MM, F_LM, F_QN, F_NN
use sim_param, only : u, v, w
implicit none

character (*), parameter :: sub = 'inflow_write'
character (*), parameter :: inflow_file = 'output/inflow_BC.out'
character (*), parameter :: field_file = 'output/inflow.vel.out'
character (*), parameter :: MPI_suffix = '.c'

integer, parameter :: lun = 80
integer, parameter :: field_lun = 81

logical, parameter :: DEBUG = .false.

character (64) :: fname

integer, save :: rec = 0
integer :: nrec
integer :: iolen
integer :: iend, iend_w

logical, save :: initialized = .false.
logical :: opn, exst

!---------------------------------------------------------------------

!--option check
if ( read_inflow_file .and. write_inflow_file ) then
  write (*, *) sub // ': cannot have read_inflow_file and write_inflow_file'
  stop
end if

!--check consistency with inflow_cond
iend = floor (buff_end * nx + 1._rprec)
iend_w = modulo (iend - 1, nx) + 1

if (.not. initialized) then

  inquire ( unit=lun, exist=exst, opened=opn )
  if ( .not. exst ) then
    write (*, *) sub // ': lun = ', lun, ' does not exist'
    stop
  end if
  if (opn) then
    write (*, *) sub // ': lun = ', lun, ' is already open'
    stop
  end if

  if ( USE_MPI ) then
      write ( fname, '(a,a,i0)' ) trim (inflow_file), MPI_suffix, coord
  else
      write ( fname, '(a)' ) inflow_file
  end if
  
  inquire ( file=fname, exist=exst, opened=opn )
  if (exst .and. opn) then
    write (*, *) sub // ': file = ', trim (fname), ' is already open'
    stop
  end if
  
  !--figure out the record length
  if ( model.eq.4 ) then
    inquire (iolength=iolen) u(iend_w,:,:), v(iend_w,:,:), w(iend_w,:,:), F_MM(1,:,:), F_LM(1,:,:)
  else if (model.eq.5) then
    inquire (iolength=iolen) u(iend_w,:,:), v(iend_w,:,:), w(iend_w,:,:), F_MM(1,:,:), F_LM(1,:,:), F_QN(1,:,:), F_NN(1,:,:)
  else
    inquire (iolength=iolen) u(iend_w,:,:), v(iend_w,:,:), w(iend_w,:,:)
  end if

  !--always add to existing inflow_file
  !--inflow_file records always start at 1
  if ( exst ) then
      !--figure out the number of records already in file
      call len_da_file (fname, iolen, nrec)
      write (*, *) sub // ': #records in ' // trim (fname) // '= ', nrec
      rec = nrec
  else
      rec = 0
  end if
  
  !--using direct-access file to allow implementation of 'inflow recycling'
  !  more easily
  !--may want to put in some simple checks on ny, nz
  open (unit=lun, file=fname, access='direct', action='write',  &
        recl=iolen)

  initialized = .true.

end if

if (jt_total == jt_start_write) then  !--write entire flow field out
  inquire (unit=field_lun, exist=exst, opened=opn)
  if (exst .and. .not. opn) then
    open (unit=field_lun, file=field_file, form='unformatted')
    call output_final (jt_total, field_lun)
  else
    write (*, *) sub // ': problem opening ' // field_file
    stop
  end if
end if

if (jt_total >= jt_start_write) then
  rec = rec + 1
  if ( model.eq.4 ) then
    write (unit=lun, rec=rec) u(iend_w,:,:), v(iend_w,:,:), w(iend_w,:,:), F_MM(1,:,:), F_LM(1,:,:)
  else if ( model.eq.5) then 
    write (unit=lun, rec=rec) u(iend_w,:,:), v(iend_w,:,:), w(iend_w,:,:), F_MM(1,:,:), F_LM(1,:,:), F_QN(1,:,:), F_NN(1,:,:)
  else
    write (unit=lun, rec=rec) u(iend_w,:,:), v(iend_w,:,:), w(iend_w,:,:)
  end if
  if ( DEBUG ) write (*, *) sub // ': wrote record ', rec
end if

end subroutine inflow_write

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine inflow_read ()
use param, only : model, ny, nz, pi, nsteps, jt_total, buff_end
use sgsmodule, only : FMM_hold, FLM_hold, FQN_hold, FNN_hold
use sim_param, only : u, v, w
implicit none

character (*), parameter :: sub = 'inflow_read'
character (*), parameter :: inflow_file = 'output/inflow_BC.out'
character (*), parameter :: debug_file = 'inflow_read_debug.dat'
character (*), parameter :: MPI_suffix = '.c'

integer, parameter :: lun = 80  !--inflow_write for now
integer, parameter :: lun_DEBUG = 88

integer, parameter :: l_blend = 300  !--length of blending zone (recycling)
                                     !--should correspond to integral scale
                                     !--this is number of t-steps
logical, parameter :: recycle = .false.

logical, parameter :: DEBUG = .false.

character (32) :: fmt
character (64) :: fname

!--check for consistency with sim_param here
!--could define a fortran integer lbz in sim_param, and make it visible
!  here, however, this may have complications elsewhere where the name lbz
!  is used.
$if ( $MPI )
    $define $lbz 0
$else
    $define $lbz 1
$endif

integer :: jy, jz
integer :: iend, iend_w
integer :: i
integer :: iolen
integer, save :: rec
integer, save :: nrec
integer :: recp

logical, save :: init_DEBUG = .false.
logical, save :: initialized = .false.
logical :: exst, opn

real (rprec) :: wgt

real (rprec) :: u_tmp(ny, $lbz:nz), v_tmp(ny, $lbz:nz), w_tmp(ny, $lbz:nz)

!---------------------------------------------------------------------

iend = floor ( buff_end * nx + 1.0_rprec )
iend_w = modulo ( iend - 1, nx ) + 1

if ( .not. initialized ) then

    inquire ( unit=lun, exist=exst, opened=opn )
    if ( .not. exst ) then
        write (*, *) sub // ': lun = ', lun, ' does not exist'
        stop
    end if
    if ( opn ) then
        write (*, *) sub // ': lun = ', lun, ' is already open'
        stop
    end if

    if ( USE_MPI ) then
        write ( fname, '(a,a,i0)' ) trim (inflow_file), MPI_suffix, coord
    else
        write ( fname, '(a)' ) inflow_file
    end if
    
    inquire ( file=fname, exist=exst, opened=opn )
    if ( exst ) then
        if ( opn ) then
            write (*, *) sub // ': file = ', fname, ' is already open'
            stop
        end if
    else
        write (*, *) sub // ': file = ', fname, ' does not exist'
        stop
    end if

    !--can only reach this point if exst and .not. opn
  
    !--figure out the record length
    if ( model.eq.4 ) then
        inquire ( iolength=iolen ) u(iend_w,:,:), v(iend_w,:,:), w(iend_w,:,:), FMM_hold, FLM_hold
    else if ( model.eq.5 ) then
        inquire ( iolength=iolen ) u(iend_w,:,:), v(iend_w,:,:), w(iend_w,:,:), FMM_hold, FMM_hold, FQN_hold, FNN_hold 
    else
        inquire ( iolength=iolen ) u(iend_w, :, :), v(iend_w, :, :), w(iend_w, :, :)
    endif
 
    !--figure out the number of records
    call len_da_file ( fname, iolen, nrec )

    write (*, *) sub // ': number of records = ', nrec

    if ( recycle ) then
        !--check minimum length requirement
        !  checks that there are some points that will be non-blended
        
        if ( 2 * (l_blend - 1) > nrec ) then
            write (*, *) sub // ': ', fname, 'is too short to recycle'
            stop
        end if
    end if

    open ( unit=lun, file=fname, access='direct', action='read',  &
           recl=iolen )

    !--file always starts a record 1, but in continued runs, we may need to
    !  access a record that is not 1 to begin with
    !--actually, with wrap-around of records, it means the reading can start
    !  at any point in the file and be OK
    !--intended use: jt_total = 1 here at start of set of runs reading
    !  from the inflow_file, so the first record read will be record 1
    rec = jt_total - 1

    initialized = .true.

end if
rec = rec + 1
if ( recycle ) then
    rec = modulo ( rec - 1, nrec - l_blend + 1 ) + 1
else
    rec = modulo ( rec - 1, nrec ) + 1
end if

if ( model.eq.4 ) then
    read ( unit=lun, rec=rec ) u(iend_w,:,:), v(iend_w,:,:), w(iend_w,:,:), FMM_hold, FLM_hold
else if ( model.eq.5 ) then
    read ( unit=lun, rec=rec ) u(iend_w,:,:), v(iend_w,:,:), w(iend_w,:,:), FMM_hold, FLM_hold, FQN_hold, FNN_hold
else
    read ( unit=lun, rec=rec ) u(iend_w, :, :), v(iend_w, :, :), w(iend_w, :, :) 
endif
if ( DEBUG ) write (*, *) sub // ' : read record ', rec
    
if ( recycle ) then
    if ( rec < l_blend ) then
        recp = nrec - l_blend + 1 + rec
        wgt = 0.5_rprec * ( 1.0_rprec -                              &
                            cos ( pi * real (rec, rprec) / l_blend ) )
            !--wgt = 0+ when rec = 1
            !  wgt = 1- when rec = l_blend
        read ( unit=lun, rec=recp ) u_tmp, v_tmp, w_tmp
        u(iend_w, :, :) = wgt * u(iend_w, :, :) + (1.0_rprec - wgt) * u_tmp
        v(iend_w, :, :) = wgt * v(iend_w, :, :) + (1.0_rprec - wgt) * v_tmp
        w(iend_w, :, :) = wgt * w(iend_w, :, :) + (1.0_rprec - wgt) * w_tmp
    end if
end if

if ( DEBUG ) then  !--write out slices as an ascii time series
    if ( .not. init_DEBUG ) then
        inquire ( unit=lun_DEBUG, exist=exst, opened=opn )
        if ( exst .and. (.not. opn) ) then
            if ( USE_MPI ) then
                open ( unit=lun_DEBUG, file=debug_file // MPI_suffix )
            else
                open ( unit=lun_DEBUG, file=debug_file )
            end if
        
            write ( lun_DEBUG, '(a)' ) 'variables = "y" "z" "t" "u" "v" "w"'
            write ( lun_DEBUG, '(3(a,i0))' ) 'zone, f=point, i= ', ny,  &
                                             ', j= ', nz,               &
                                             ', k= ', nsteps
        else
            write (*, *) sub // ': problem opening debug file'
            stop
        end if
        init_DEBUG = .true.
    end if

    fmt = '(3(1x,i0),3(1x,es12.5))'
    do jz = 1, nz
        do jy = 1, ny
            write ( lun_DEBUG, fmt ) jy, jz, jt_total, u(iend_w, jy, jz),  &
                                     v(iend_w, jy, jz), w(iend_w, jy, jz)
        end do
    end do
end if

end subroutine inflow_read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--finds number of records on existing direct-access unformatted file
!--taken from Clive Page's comp.lang.fortran posting (12/16/2003), 
!  under the topic counting number of records in a Fortran direct file
!--minor changes/renaming
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine len_da_file(fname, lenrec, length)
implicit none
character (*), intent(in) :: fname  ! name of existing direct-access file
integer, intent(in)       :: lenrec ! record length (O/S dependent units)
integer, intent(out) :: length      ! number of records.
!
character (1) :: cdummy
integer :: lunit, nlo, nhi, mid, kode
logical :: exists, open
!
! find a free unit on which to open the file
!
do lunit = 99, 1, -1
  !--units to skip (compiler dependent)
  select case (lunit)
    case (5:6)
      !--do nothing
    case default
      inquire(unit=lunit, exist=exists, opened=open)
      if(exists .and. .not. open) exit
  end select
end do
open(unit=lunit, file=fname, access="direct", recl=lenrec, iostat=kode)
if(kode /= 0) then
  print *, 'error in len_da_file: ', trim(fname), ' does not exist'
  return
end if
!
! expansion phase
!
mid = 1
do
  read(lunit, rec=mid, iostat=kode) cdummy
  if(kode /= 0) exit
  mid = 2 * mid
end do
!
! length is between mid/2 and mid, do binary search to refine
!
nlo = mid/2
nhi = mid
do while(nhi - nlo > 1)
  mid = (nlo + nhi) / 2
  read(lunit, rec=mid, iostat=kode) cdummy
  if(kode == 0) then
     nlo = mid
  else
     nhi = mid
  end if
end do
length = nlo
close(unit=lunit)
return
end subroutine len_da_file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module io
