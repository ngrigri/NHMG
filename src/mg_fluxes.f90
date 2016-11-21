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

    integer(kind=ip):: k,j,i
    integer(kind=ip):: nx,ny,nz

    real(kind=rp), dimension(:,:)  , pointer :: dx,dy
    real(kind=rp), dimension(:,:),   pointer :: dxu,dyv
    real(kind=rp), dimension(:,:,:), pointer :: zr,zw
    real(kind=rp), dimension(:,:,:), pointer :: dzw
    real(kind=rp), dimension(:,:,:), pointer :: Arx,Ary
    real(kind=rp), dimension(:,:,:), pointer :: zxdy,zydx
    real(kind=rp), dimension(:,:,:), pointer :: cw

    real(kind=rp), dimension(:,:,:), pointer :: uf_tmp,vf_tmp !!dirty

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz

    dx => grid(1)%dx
    dy => grid(1)%dy
    dxu => grid(1)%dxu
    dyv => grid(1)%dyv
    zr => grid(1)%zr
    zw => grid(1)%zw
    dzw => grid(1)%dzw
    Arx => grid(1)%Arx
    Ary => grid(1)%Ary
    cw  => grid(1)%cw
    zxdy => grid(1)%zxdy
    zydx => grid(1)%zydx

    !- UF -!

    ! lower level
    k = 1
    do i = 1,nx+1  
       do j = 1,ny 
!       do j = 0,ny+1

          uf(k,j,i) = Arx(k,j,i)/dxu(j,i) *dxu(j,i)*u(k,j,i) &
              - qrt * ( &
              + zxdy(k,j,i  )* two *dzw(k  ,j,i  )*w(k  ,j,i  ) & 
              + zxdy(k,j,i  )      *dzw(k+1,j,i  )*w(k+1,j,i  ) &
              + zxdy(k,j,i-1)* two *dzw(k  ,j,i-1)*w(k  ,j,i-1) &
              + zxdy(k,j,i-1)      *dzw(k+1,j,i-1)*w(k+1,j,i-1) )  
       enddo
    enddo

    !interior levels
    do i = 1,nx+1  
       do j = 1,ny 
          do k = 2,nz-1 
 
             uf(k,j,i) = Arx(k,j,i)/dxu(j,i) *dxu(j,i)*u(k,j,i) &
                  - qrt * ( &
                  + zxdy(k,j,i  ) *dzw(k  ,j,i  )*w(k  ,j,i  ) &
                  + zxdy(k,j,i  ) *dzw(k+1,j,i  )*w(k+1,j,i  ) &
                  + zxdy(k,j,i-1) *dzw(k  ,j,i-1)*w(k  ,j,i-1) &
                  + zxdy(k,j,i-1) *dzw(k+1,j,i-1)*w(k+1,j,i-1) )
          enddo
       enddo
    enddo

    k = nz !upper level
    do i = 1,nx+1  
       do j = 1,ny 

           uf(k,j,i) = Arx(k,j,i)/dxu(j,i) *dxu(j,i)*u(k,j,i) &

                  - qrt * ( zxdy(k,j,i  ) *      dzw(k  ,j,i  )*w(k  ,j,i  ) &
                           +zxdy(k,j,i  ) * two *dzw(k+1,j,i  )*w(k+1,j,i  ) &
                           +zxdy(k,j,i-1) *      dzw(k  ,j,i-1)*w(k  ,j,i-1) &
                           +zxdy(k,j,i-1) * two *dzw(k+1,j,i-1)*w(k+1,j,i-1) )
       enddo
    enddo

    allocate(uf_tmp(1:nz,0:ny+1,0:nx+1))
    uf_tmp(:,1:ny,1:nx+1)=uf(:,1:ny,1:nx+1)
    call fill_halo(1,uf_tmp,lbc_null='u')
    uf(:,0:ny+1,1:nx+1)=uf_tmp(:,0:ny+1,1:nx+1)
    deallocate(uf_tmp)

    !- VF -!

    !lower level
    k = 1
    do i = 1,nx
!    do i = 0,nx+1
       do j = 1,ny+1

          vf(k,j,i) = Ary(k,j,i)/dyv(j,i) *dyv(j,i)*v(k,j,i) &
               - qrt * ( &
               + zydx(k,j  ,i)* two *dzw(k  ,j  ,i)*w(k  ,j  ,i) &  
               + zydx(k,j  ,i)      *dzw(k+1,j  ,i)*w(k+1,j  ,i) &
               + zydx(k,j-1,i)* two *dzw(k  ,j-1,i)*w(k  ,j-1,i) &
               + zydx(k,j-1,i)      *dzw(k+1,j-1,i)*w(k+1,j-1,i) )
       enddo
    enddo

    !interior levels
    do i = 1,nx
       do j = 1,ny+1
          do k = 2,nz-1

             vf(k,j,i) = Ary(k,j,i)/dyv(j,i) *dyv(j,i)*v(k,j,i) &
                  - qrt * ( &
                  + zydx(k,j  ,i) *dzw(k  ,j  ,i)*w(k  ,j  ,i) &
                  + zydx(k,j  ,i) *dzw(k+1,j  ,i)*w(k+1,j  ,i) &
                  + zydx(k,j-1,i) *dzw(k  ,j-1,i)*w(k  ,j-1,i) &
                  + zydx(k,j-1,i) *dzw(k+1,j-1,i)*w(k+1,j-1,i) )
          enddo
       enddo
    enddo

    k = nz !upper level
    do i = 1,nx
       do j = 1,ny+1

            vf(k,j,i) = Ary(k,j,i)/dyv(j,i) *dyv(j,i)*v(k,j,i) &
                  - qrt * ( &
                  + zydx(k,j  ,i)       *dzw(k  ,j  ,i)*w(k  ,j  ,i) &
                  + zydx(k,j  ,i) * two *dzw(k+1,j  ,i)*w(k+1,j  ,i) &
                  + zydx(k,j-1,i)       *dzw(k  ,j-1,i)*w(k  ,j-1,i) &
                  + zydx(k,j-1,i) * two *dzw(k+1,j-1,i)*w(k+1,j-1,i) ) 
       enddo
    enddo

    allocate(vf_tmp(1:nz,0:ny+1,0:nx+1))
    vf_tmp(:,1:ny+1,1:nx)=vf(:,1:ny+1,1:nx)
    call fill_halo(1,vf_tmp,lbc_null='v')
    vf(:,1:ny+1,0:nx+1)=vf_tmp(:,1:ny+1,0:nx+1)
    deallocate(vf_tmp)

  end subroutine set_fluxes
 
end module mg_fluxes

