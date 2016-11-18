module mg_correct_uvw

  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids
  use mg_mpi_exchange
  use mg_netcdf_out

  implicit none

contains
  !-------------------------------------------------------------------------     
  subroutine correct_uvw(u,v,w)

    real(kind=rp), dimension(:,:,:), pointer, intent(inout) :: u,v,w

    integer(kind=ip):: k, j, i
    integer(kind=ip):: nx, ny, nz

    real(kind=rp), dimension(:,:)  , pointer :: dx,dy
    real(kind=rp), dimension(:,:,:), pointer :: zr,zw
    real(kind=rp), dimension(:,:,:), pointer :: p
    real(kind=rp) :: dxu,dyv
    real(kind=rp) :: dzw

    real(kind=rp), parameter :: two  = 2._rp
    real(kind=rp), parameter :: one  = 1._rp
    real(kind=rp), parameter :: hlf  = 0.5_rp
    real(kind=rp), parameter :: qrt  = 0.25_rp
    real(kind=rp), parameter :: zero  = 0._rp

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz

    dx => grid(1)%dx
    dy => grid(1)%dy
    zr => grid(1)%zr
    zw => grid(1)%zw

    if (myrank==0) write(*,*)'- correct u,v,w:'

    !! Correct
    p => grid(1)%p

    do i = 1,nx+1
!    do i = 2,nx
!       do j = 0,ny+1
        do j = 1,ny
          do k = 1,nz

             dxu = hlf * (dx(j,i)+dx(j,i-1))

             u(k,j,i) = u(k,j,i) - one / dxu * (p(k,j,i)-p(k,j,i-1))

          enddo
       enddo
    enddo

!    do i = 0,nx+1
    do i = 1,nx
       do j = 1,ny+1 
!       do j = 2,ny 
          do k = 1,nz

             dyv = hlf * (dy(j,i)+dy(j-1,i))

             v(k,j,i) = v(k,j,i) - one / dyv * (p(k,j,i)-p(k,j-1,i  ))

          enddo
       enddo
    enddo

!    do i = 0,nx+1
    do i = 1,nx
!       do j = 0,ny+1
       do j = 1,ny

          k = 1 !bottom

          do k = 2,nz !interior levels
             dzw = zr(k,j,i)-zr(k-1,j,i)
             w(k,j,i) = w(k,j,i) - one / dzw * (p(k,j,i)-p(k-1,j,i))
          enddo

          k = nz+1 !surface
          dzw = zw(nz+1,j,i)-zr(nz,j,i)
          w(k,j,i) = w(k,j,i) - one / dzw * (-p(k-1,j,i))

       enddo
    enddo

  end subroutine correct_uvw

end module mg_correct_uvw
