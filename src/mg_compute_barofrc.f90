module mg_compute_barofrc

  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids
  use mg_mpi_exchange
  use mg_netcdf_out

  implicit none

contains
  !-------------------------------------------------------------------------     
  subroutine compute_barofrc(dt,ru,rv)

    real(kind=rp),                            intent(in)  :: dt
    real(kind=rp), dimension(:,:)  , pointer, intent(out) :: ru,rv

    integer(kind=ip):: k, j, i
    integer(kind=ip):: nx, ny, nz

    real(kind=rp), dimension(:,:)  , pointer :: dx,dy
    real(kind=rp), dimension(:,:,:), pointer :: zr,zw

    real(kind=rp), dimension(:,:,:), pointer :: dzw
    real(kind=rp), dimension(:,:,:), pointer :: Arx,Ary
    real(kind=rp), dimension(:,:,:), pointer :: zxdy,zydx

    real(kind=rp), dimension(:,:,:), pointer :: p

    real(kind=rp), parameter :: two  = 2._rp
    real(kind=rp), parameter :: one  = 1._rp
    real(kind=rp), parameter :: hlf  = 0.5_rp
    real(kind=rp), parameter :: qrt  = 0.25_rp
    real(kind=rp), parameter :: zero = 0._rp

    if (myrank==0) write(*,*)'- compute barofrc:'

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz

    dx => grid(1)%dx
    dy => grid(1)%dy
    zr => grid(1)%zr
    zw => grid(1)%zw

    Arx => grid(1)%Arx
    Ary => grid(1)%Ary
    zxdy => grid(1)%zxdy
    zydx => grid(1)%zydx

    !! Compute
    p => grid(1)%p

    ru(:,:) = 0._8
    rv(:,:) = 0._8

    do i = 1,nx+1
       do j = 0,ny+1 
!       do j = 1,ny

          k = 1

             ru(i,j) = ru(i,j) &

                     - Arx(k,j,i) * (p(k,j,i)-p(k,j,i-1)) &

                     + qrt * (  & ! zero in our flat bottom case
                             + dx(j,i  )*zxdy(k,j,i  ) * (p(k+1,j,i)-p(k  ,j,i)) & 
                             +  & ! zero in our flat bottom case
                             + dx(j,i-1)*zxdy(k,j,i-1) * (p(k+1,j,i-1)-p(k  ,j,i-1)) &
                             )

          do k = 2,nz-1

             ru(i,j) = ru(i,j) &

                     - Arx(k,j,i) * (p(k,j,i)-p(k,j,i-1)) &

                     + qrt * ( dx(j,i  )*zxdy(k,j,i  ) * (p(k  ,j,i)-p(k-1,j,i)) & 
                             + dx(j,i  )*zxdy(k,j,i  ) * (p(k+1,j,i)-p(k  ,j,i)) & 
                             + dx(j,i-1)*zxdy(k,j,i-1) * (p(k  ,j,i-1)-p(k-1,j,i-1)) & 
                             + dx(j,i-1)*zxdy(k,j,i-1) * (p(k+1,j,i-1)-p(k  ,j,i-1)) &
                             )

          enddo

          k = nz

             ru(i,j) = ru(i,j) &

                     - Arx(k,j,i) * (p(k,j,i)-p(k,j,i-1)) &

                     + qrt * ( dx(j,i  )*zxdy(k,j,i  ) * (p(k  ,j,i)-p(k-1,j,i)) & 
                             + dx(j,i  )*zxdy(k,j,i  ) * two * (-p(k  ,j,i)) & 
                             + dx(j,i-1)*zxdy(k,j,i-1) * (p(k  ,j,i-1)-p(k-1,j,i-1)) & 
                             + dx(j,i-1)*zxdy(k,j,i-1) * two * (-p(k  ,j,i-1)) &
                             )

          ru(i,j) = ru(i,j)/dt

       enddo
    enddo
 
    do i = 0,nx+1
!    do i = 1,nx
       do j = 1,ny+1 

          k = 1

             rv(i,j) = rv(i,j) &

                     - Ary(k,j,i) * (p(k,j,i)-p(k,j-1,i)) &

                     + qrt * (  & ! zero in our flat bottom case
                             + dy(j  ,i)*zydx(k,j  ,i) * (p(k+1,j,i)-p(k  ,j,i)) & 
                             +  & ! zero in our flat bottom case
                             + dy(j-1,i)*zydx(k,j-1,i) * (p(k+1,j-1,i)-p(k  ,j-1,i)) &
                             )

          do k = 2,nz-1

             rv(i,j) = rv(i,j) &

                     - Ary(k,j,i) * (p(k,j,i)-p(k,j-1,i)) &

                     + qrt * ( dy(j  ,i)*zydx(k,j  ,i) * (p(k  ,j,i)-p(k-1,j,i)) & 
                             + dy(j  ,i)*zydx(k,j  ,i) * (p(k+1,j,i)-p(k  ,j,i)) & 
                             + dy(j-1,i)*zydx(k,j-1,i) * (p(k  ,j-1,i)-p(k-1,j-1,i)) & 
                             + dy(j-1,i)*zydx(k,j-1,i) * (p(k+1,j-1,i)-p(k  ,j-1,i)) &
                             )

          enddo

          k = nz

             rv(i,j) = rv(i,j) &

                     - Ary(k,j,i) * (p(k,j,i)-p(k,j-1,i)) &

                     + qrt * ( dy(j  ,i)*zydx(k,j  ,i) * (p(k  ,j,i)-p(k-1,j,i)) & 
                             + dy(j  ,i)*zydx(k,j  ,i) * two * (-p(k  ,j,i)) & 
                             + dy(j-1,i)*zydx(k,j-1,i) * (p(k  ,j-1,i)-p(k-1,j-1,i)) & 
                             + dy(j-1,i)*zydx(k,j-1,i) * two * (-p(k  ,j-1,i)) &
                             )

          rv(i,j) = rv(i,j)/dt

       enddo
    enddo

  end subroutine compute_barofrc

end module mg_compute_barofrc
