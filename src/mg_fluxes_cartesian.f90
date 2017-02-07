module mg_fluxes

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
  subroutine set_fluxes(u,v,w,uf,vf)

    real(kind=rp), dimension(:,:,:), pointer, intent(in)  :: u,v,w
    real(kind=rp), dimension(:,:,:), pointer, intent(out) :: uf,vf

    integer(kind=ip) :: k,j,i
    integer(kind=ip) :: nx,ny,nz

    real(kind=rp), dimension(:,:,:), pointer :: Arx,Ary

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz

    Arx   => grid(1)%Arx
    Ary   => grid(1)%Ary

    !- uf -!

    do i = 1,nx+1  
       do j = 0,ny+1

          do k = 1,nz
             uf(k,j,i) = Arx(k,j,i) * u(k,j,i)
          enddo

       enddo
    enddo

    !- vf -!

    do i = 0,nx+1
       do j = 1,ny+1

          do k = 1,nz
             vf(k,j,i) = Ary(k,j,i) * v(k,j,i)
          enddo

       enddo
    enddo

  end subroutine set_fluxes
 
end module mg_fluxes

