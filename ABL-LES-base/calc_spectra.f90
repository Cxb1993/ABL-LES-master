subroutine calc_spectra(u_c,spec)
! Calculates fft of u and then calculates 
! fft(u) * complex conjugate of fft(u)

use types,only:rprec
use fft
use param,only:lh,nx,ny,nz,ld
implicit none

real(kind=rprec),dimension(ld,ny),intent(in)::u_c
real(kind=rprec),dimension(nx/2),intent(out)::spec
! assumes Nyquist is 0

integer::jy,jz,k
real(kind=rprec),dimension(nx)::vel_r,vel_c

integer*8,save::plan
logical,save::init = .false.

if (.not. init) then
  call rfftw_f77_create_plan(plan,nx,FFTW_REAL_TO_COMPLEX,FFTW_MEASURE)
  init = .true.
end if

spec(:) = 0._rprec ! initialize
do jy=1,ny
   vel_r(:) = u_c(1:nx,jy)/real(nx,kind=rprec)

   ! check this normaliztion-part of forward; call the fft
   call rfftw_f77_one(plan,vel_r,vel_c)

   ! compute magnitudes the 0.5 is the 1/2, all others are taken care of
   spec(1) = spec(1) + 0.5_rprec*vel_c(1)*vel_c(1)
   do k=2,nx/2
      spec(k) = spec(k) + vel_c(k)*vel_c(k) + vel_c(nx+2-k)*vel_c(nx+2-k)
   end do

   ! spec(nx/2+1) = spec(nx/2+1) + vel_c(nx/2+1)*vel_c(nx/2+1)
   ! Assume Nyquist is zero
end do
spec(:) = 2._rprec*spec(:)/real(ny,kind=rprec) ! for average over Ny

end subroutine calc_spectra
