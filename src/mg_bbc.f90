module mg_bbc

  use mg_cst
  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids
  use mg_mpi_exchange
  use mg_netcdf_out

  implicit none

contains
  !-------------------------------------------------------------------------     
  subroutine set_bbc(u,v,w)

    real(kind=rp), dimension(:,:,:), pointer, intent(in)     :: u,v
    real(kind=rp), dimension(:,:,:), pointer, intent(inout)  :: w

    integer(kind=ip):: k,j,i
    integer(kind=ip):: nx,ny,nz

    real(kind=rp), dimension(:,:)  , pointer :: dx,dy
    real(kind=rp), dimension(:,:,:), pointer :: zr, dz
    real(kind=rp), dimension(:,:,:), pointer :: dzw
    real(kind=rp), dimension(:,:,:), pointer :: zxdy,zydx
    real(kind=rp), dimension(:,:,:), pointer :: alpha

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz

    dx    => grid(1)%dx
    dy    => grid(1)%dy
    zr    => grid(1)%zr
    dz    => grid(1)%dz
    dzw   => grid(1)%dzw
    alpha => grid(1)%alpha
    zxdy  => grid(1)%zxdy
    zydx  => grid(1)%zydx

    !- bottom vertical momentum -!

!! TODO !! 23 nov 2016
!!$    do i = 1,nx  
!!$       do j = 1,ny
!!$
!!$          k = 1 ! bottom
!!$          w(k,j,i) = ( &
!!$               + hlf * hlf * ( &
!!$               + zxdy(k,j,i) * ( dx(j,i  ) + dx(j,i-1) ) * u(k,j,i  ) &
!!$               + zxdy(k,j,i) * ( dx(j,i+1) + dx(j,i  ) ) * u(k,j,i+1) ) &
!!$               + hlf * hlf * ( &
!!$               + zydx(k,j,i) * ( dy(j  ,i) + dy(j-1,i) ) * v(k,j  ,i) &
!!$               + zydx(k,j,i) * ( dy(j+1,i) + dy(j  ,i) ) * v(k,j+1,i) ) &
!!$                      ) / (cw(k,j,i)*dzw(k,j,i))
!!$
!!$       enddo
!!$    enddo
!! TODO !! 23 nov 2016

    do i = 1,nx  
       do j = 1,ny

          k = 1 ! bottom
          w(k,j,i) = 99999999999999999._rp

       enddo
    enddo


  end subroutine set_bbc
 
end module mg_bbc
