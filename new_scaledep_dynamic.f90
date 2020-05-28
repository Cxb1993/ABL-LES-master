subroutine new_scaledep_dynamic(Cs_opt2,S11,S12,S13,S22,S23,S33)
! Subroutine written by Stimit Shah on 7th Sept 2010
! Returns scale dependent planar averaged Cs 

use types,only:rprec
use param ! ,only:ld,nx,ny,nz,dx,dy,dz,jt,USE_MPI,coord
use sim_param,only:path,u,v,w
use test_filtermodule

implicit none
integer::jz,tagCs=400
real(kind=rprec),dimension(ld,ny,nz)::S11,S12,S13,S22,S23,S33
real(kind=rprec),dimension(ld,ny)::L11,L12,L13,L22,L23,L33
real(kind=rprec),dimension(ld,ny)::M11,M12,M13,M22,M23,M33
real(kind=rprec),dimension(ld,ny)::Q11,Q12,Q13,Q22,Q23,Q33
real(kind=rprec),dimension(ld,ny)::N11,N12,N13,N22,N23,N33
real(kind=rprec)::LM,MM,QN,NN
real(kind=rprec),dimension(nz),intent(out)::Cs_opt2

real(kind=rprec),dimension(ld,ny)::S_bar,S11_bar,S12_bar, &
     S13_bar,S22_bar,S23_bar,S33_bar,S_S11_bar,S_S12_bar, &
     S_S13_bar,S_S22_bar,S_S23_bar,S_S33_bar
real(kind=rprec),dimension(ld,ny)::S_hat,S11_hat,S12_hat, &
     S13_hat,S22_hat,S23_hat,S33_hat,S_S11_hat,S_S12_hat, &
     S_S13_hat,S_S22_hat,S_S23_hat,S_S33_hat

real(kind=rprec),dimension(ld,ny)::u_bar,v_bar,w_bar
real(kind=rprec),dimension(ld,ny)::u_hat,v_hat,w_hat
real(kind=rprec),dimension(ld,ny)::S
real(kind=rprec),dimension(nz)::Cs_opt2_2d,Cs_opt2_4d,beta
real(kind=rprec)::delta,const,betaclip
real(kind=rprec)::tf1,tf2,tf1_2,tf2_2,z
character(len=24)::fname

$if ($MPI)
  $define $lbz 0
  integer::recvcounts(nproc)
  integer::displs(nproc)
$else
  $define $lbz 1
$endif

delta = filter_size*(dx*dy*dz)**(1._rprec/3._rprec)
tf1=2._rprec
tf2=4._rprec
tf1_2=tf1**2
tf2_2=tf2**2

const = 2._rprec*(delta**2)

do jz=1,nz
   if ( ((.not. USE_MPI) .or. (USE_MPI .and. coord==0)) .and. (jz==1) ) then
     ! watch the 0.25's: w = c*z**z near wall, so get 0.25
     L11(:,:) = u(:,:,1)*u(:,:,1) ! uv-node
     L12(:,:) = u(:,:,1)*v(:,:,1) ! uv-node
     L13(:,:) = u(:,:,1)*0.25_rprec*w(:,:,2)  ! assume parabolic near wall
     L22(:,:) = v(:,:,1)*v(:,:,1) ! uv-node
     L23(:,:) = v(:,:,1)*0.25_rprec*w(:,:,2)  ! uv-node
     L33(:,:) = (0.25_rprec*w(:,:,2))**2      ! uv-node
     u_bar(:,:) = u(:,:,1)
     v_bar(:,:) = v(:,:,1)
     w_bar(:,:) = 0.25_rprec*w(:,:,2)
   else  ! w-nodes
     L11(:,:) = 0.5_rprec*(u(:,:,jz) + u(:,:,jz-1))*&
                0.5_rprec*(u(:,:,jz) + u(:,:,jz-1))
     L12(:,:) = 0.5_rprec*(u(:,:,jz) + u(:,:,jz-1))*&
                0.5_rprec*(v(:,:,jz) + v(:,:,jz-1))
     L13(:,:) = 0.5_rprec*(u(:,:,jz) + u(:,:,jz-1))*w(:,:,jz)
     L22(:,:) = 0.5_rprec*(v(:,:,jz) + v(:,:,jz-1))*&
                0.5_rprec*(v(:,:,jz) + v(:,:,jz-1))
     L23(:,:) = 0.5_rprec*(v(:,:,jz) + v(:,:,jz-1))*w(:,:,jz)
     L33(:,:) = w(:,:,jz)*w(:,:,jz)
     u_bar(:,:) = 0.5_rprec*(u(:,:,jz) + u(:,:,jz-1))
     v_bar(:,:) = 0.5_rprec*(v(:,:,jz) + v(:,:,jz-1))
     w_bar(:,:) = w(:,:,jz)
   end if
   u_hat = u_bar
   v_hat = v_bar
   w_hat = w_bar
   
   call test_filter(u_bar,G_test)
   call test_filter(v_bar,G_test)
   call test_filter(w_bar,G_test)
   call test_filter(L11,G_test)  ! in-place filtering
   L11 = L11 - u_bar*u_bar
   call test_filter(L12,G_test)
   L12 = L12 - u_bar*v_bar
   call test_filter(L13,G_test)
   L13 = L13 - u_bar*w_bar
   call test_filter(L22,G_test)
   L22 = L22 - v_bar*v_bar
   call test_filter(L23,G_test)
   L23 = L23 - v_bar*w_bar
   call test_filter(L33,G_test)
   L33 = L33 - w_bar*w_bar
   
   Q11 = u_bar*u_bar
   Q12 = u_bar*v_bar
   Q13 = u_bar*w_bar
   Q22 = v_bar*v_bar
   Q23 = v_bar*w_bar
   Q33 = w_bar*w_bar

   call test_filter(u_hat,G_test_test)
   call test_filter(v_hat,G_test_test)
   call test_filter(w_hat,G_test_test)
   call test_filter(Q11,G_test_test)
   Q11 = Q11 - u_hat*u_hat
   call test_filter(Q12,G_test_test)
   Q12 = Q12 - u_hat*v_hat
   call test_filter(Q13,G_test_test)
   Q13 = Q13 - u_hat*w_hat
   call test_filter(Q22,G_test_test)
   Q22 = Q22 - v_hat*v_hat
   call test_filter(Q23,G_test_test)
   Q23 = Q23 - v_hat*w_hat
   call test_filter(Q33,G_test_test)
   Q33 = Q33 - w_hat*w_hat
   
! calculate |S|
   S(:,:) = sqrt(2._rprec*(S11(:,:,jz)**2 + S22(:,:,jz)**2 + & 
                 S33(:,:,jz)**2 + 2._rprec*(S12(:,:,jz)**2 + &
                 S13(:,:,jz)**2 + S23(:,:,jz)**2)))

! already on w-nodes
   S11_bar(:,:) = S11(:,:,jz)  
   S12_bar(:,:) = S12(:,:,jz)  
   S13_bar(:,:) = S13(:,:,jz)  
   S22_bar(:,:) = S22(:,:,jz)  
   S23_bar(:,:) = S23(:,:,jz)  
   S33_bar(:,:) = S33(:,:,jz)  
   
   S11_hat = S11_bar
   S12_hat = S12_bar
   S13_hat = S13_bar
   S22_hat = S22_bar
   S23_hat = S23_bar
   S33_hat = S33_bar
   
   call test_filter(S11_bar,G_test)
   call test_filter(S12_bar,G_test)
   call test_filter(S13_bar,G_test)
   call test_filter(S22_bar,G_test)
   call test_filter(S23_bar,G_test)
   call test_filter(S33_bar,G_test)

   call test_filter(S11_hat,G_test_test)
   call test_filter(S12_hat,G_test_test)
   call test_filter(S13_hat,G_test_test)
   call test_filter(S22_hat,G_test_test)
   call test_filter(S23_hat,G_test_test)
   call test_filter(S33_hat,G_test_test)

   S_bar = sqrt(2._rprec*(S11_bar**2 + S22_bar**2 + S33_bar**2 + &
                2._rprec*(S12_bar**2 + S13_bar**2 + S23_bar**2)))

   S_hat = sqrt(2._rprec*(S11_hat**2 + S22_hat**2 + S33_hat**2 + &
                2._rprec*(S12_hat**2 + S13_hat**2 + S23_hat**2)))

   S_S11_bar = S*S11(:,:,jz)
   S_S12_bar = S*S12(:,:,jz)
   S_S13_bar = S*S13(:,:,jz)
   S_S22_bar = S*S22(:,:,jz)
   S_S23_bar = S*S23(:,:,jz)
   S_S33_bar = S*S33(:,:,jz)

   S_S11_hat = S_S11_bar
   S_S12_hat = S_S12_bar
   S_S13_hat = S_S13_bar
   S_S22_hat = S_S22_bar
   S_S23_hat = S_S23_bar
   S_S33_hat = S_S33_bar

   call test_filter(S_S11_bar,G_test)
   call test_filter(S_S12_bar,G_test)
   call test_filter(S_S13_bar,G_test)
   call test_filter(S_S22_bar,G_test)
   call test_filter(S_S23_bar,G_test)
   call test_filter(S_S33_bar,G_test)     
   
   call test_filter(S_S11_hat,G_test_test)
   call test_filter(S_S12_hat,G_test_test)
   call test_filter(S_S13_hat,G_test_test)
   call test_filter(S_S22_hat,G_test_test)
   call test_filter(S_S23_hat,G_test_test)
   call test_filter(S_S33_hat,G_test_test)     

   M11 = const*(S_S11_bar - tf1_2*S_bar*S11_bar)
   M12 = const*(S_S12_bar - tf1_2*S_bar*S12_bar)
   M13 = const*(S_S13_bar - tf1_2*S_bar*S13_bar)
   M22 = const*(S_S22_bar - tf1_2*S_bar*S22_bar)
   M23 = const*(S_S23_bar - tf1_2*S_bar*S23_bar)
   M33 = const*(S_S33_bar - tf1_2*S_bar*S33_bar)

   N11 = const*(S_S11_hat - tf2_2*S_hat*S11_hat)
   N12 = const*(S_S12_hat - tf2_2*S_hat*S12_hat)
   N13 = const*(S_S13_hat - tf2_2*S_hat*S13_hat)
   N22 = const*(S_S22_hat - tf2_2*S_hat*S22_hat)
   N23 = const*(S_S23_hat - tf2_2*S_hat*S23_hat)
   N33 = const*(S_S33_hat - tf2_2*S_hat*S33_hat)

   LM = sum(L11*M11 + L22*M22 + L33*M33 + 2._rprec*(L12*M12 + L13*M13 + L23*M23))
   MM = sum(M11*M11 + M22*M22 + M33*M33 + 2._rprec*(M12*M12 + M13*M13 + M23*M23))
   QN = sum(Q11*N11 + Q22*N22 + Q33*N33 + 2._rprec*(Q12*N12 + Q13*N13 + Q23*N23))
   NN = sum(N11*N11 + N22*N22 + N33*N33 + 2._rprec*(N12*N12 + N13*N13 + N23*N23))

   ! LM = max(real(1E-24),real(LM))
   ! QN = max(real(1E-24),real(QN))

   Cs_opt2_2d(jz) = LM/MM
   Cs_opt2_4d(jz) = QN/NN

   Cs_opt2_2d(jz) = max(real(1E-24),real(Cs_opt2_2d(jz)))
   Cs_opt2_4d(jz) = max(real(1E-24),real(Cs_opt2_4d(jz)))

   beta(jz)    = Cs_opt2_4d(jz)/Cs_opt2_2d(jz)
   betaclip    = max(real(beta(jz)),real(1._rprec/8._rprec))

   Cs_opt2(jz) = Cs_opt2_2d(jz)/betaclip
   Cs_opt2(jz) = max(0._rprec,real(Cs_opt2(jz),kind=rprec))
end do

end subroutine new_scaledep_dynamic
