module mg_coupling

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
  subroutine bc2bt_coupling(dt,ru,rv)

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

    if (myrank==0) write(*,*)'- bc2bt_coupling :'

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

    ru(:,:) = zero
    rv(:,:) = zero

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

  end subroutine bc2bt_coupling

  !-------------------------------------------------------------------------     
  subroutine bt2bc_coupling(uf_bar,vf_bar,u,v,w,uf,vf)

    real(kind=rp), dimension(:,:),   pointer, intent(in)  :: uf_bar,vf_bar
    real(kind=rp), dimension(:,:,:), pointer, intent(inout)  :: u,v,w
    real(kind=rp), dimension(:,:,:), pointer, optional, intent(out):: uf,vf

    integer(kind=ip):: k,j,i
    integer(kind=ip):: nx,ny,nz

    real(kind=rp), dimension(:,:)  , pointer :: dx,dy
    real(kind=rp), dimension(:,:,:), pointer :: zr,zw

    real(kind=rp), dimension(:,:,:), pointer :: dzw
    real(kind=rp), dimension(:,:,:), pointer :: Arx,Ary
    real(kind=rp), dimension(:,:,:), pointer :: zydx,zxdy
    real(kind=rp), dimension(:,:,:), pointer :: cw

    real(kind=rp), dimension(:,:,:),   pointer :: uf_tmp,vf_tmp !!dirty

    real(kind=rp), dimension(:,:)  ,   pointer :: su_integr,sv_integr
    real(kind=rp), dimension(:,:)  ,   pointer :: uf_integr,vf_integr
    real(kind=rp), dimension(:,:,:),   pointer :: uf_integr_tmp,vf_integr_tmp !!dirty
    real(kind=rp), dimension(:,:,:),   pointer :: wc

    integer(kind=ip), save :: iter_coupling_in=0
    iter_coupling_in = iter_coupling_in + 1

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz

    dx => grid(1)%dx
    dy => grid(1)%dy
    zr => grid(1)%zr
    zw => grid(1)%zw

    dzw => grid(1)%dzw
    Arx => grid(1)%Arx
    Ary => grid(1)%Ary
    cw  => grid(1)%cw
    zxdy => grid(1)%zxdy
    zydx => grid(1)%zydx

    !------------------------------------------
    ! compute integrated transport anomalies

    allocate(su_integr(0:ny+1,0:nx+1))
    allocate(sv_integr(0:ny+1,0:nx+1))
    allocate(uf_integr(0:ny+1,0:nx+1))
    allocate(vf_integr(0:ny+1,0:nx+1))
    su_integr(:,:) = zero
    sv_integr(:,:) = zero
    uf_integr(:,:) = zero
    vf_integr(:,:) = zero

    do i = 1,nx+1  
       do j = 1,ny 
!       do j = 0,ny +1

          k = 1 ! lower level
          su_integr(j,i) =  Arx(k,j,i)
          uf_integr(j,i) =  Arx(k,j,i) * u(k,j,i) &
               - qrt * ( &
               + zxdy(k,j,i  )* two * dzw(k  ,j,i  )*w(k  ,j,i  ) & 
               + zxdy(k,j,i  )*       dzw(k+1,j,i  )*w(k+1,j,i  ) & 
               + zxdy(k,j,i-1)* two * dzw(k  ,j,i-1)*w(k  ,j,i-1) &
               + zxdy(k,j,i-1)*       dzw(k+1,j,i-1)*w(k+1,j,i-1) &
               )  ! umask

!!$               - qrt * ( &
!!$               + zxdy(k,j,i  )*dzw(k+1,j,i  )*w(k+1-1,j,i)   & !note the index trick for w
!!$               + zxdy(k,j,i-1)*dzw(k+1,j,i-1)*w(k+1-1,j,i-1) ) & !note the index trick for w
!!$               -( &
!!$               + zxdy(k,j,i  )*zxdy(k,j,i  )/(cw(k,j,i  )+cw(k+1,j,i  )) &
!!$               + zxdy(k,j,i-1)*zxdy(k,j,i-1)/(cw(k,j,i-1)+cw(k+1,j,i-1)) ) * &
!!$               (hlf * ( dx(j,i) + dx(j,i-1) )) * u(k,j,i) &
!!$               -( & 
!!$               + zxdy(k,j,i) * zydx(k,j,i)/(cw(k,j,i  )+cw(k+1,j,i  ))   &
!!$               * hlf * ( &
!!$               hlf * ( dy(j  ,i) + dy(j-1,i) )*v(k,j,i ) + &
!!$               hlf * ( dy(j+1,i) + dy(j  ,i) )*v(k,j+1,i ))   & 
!!$               + zxdy(k,j,i-1) * zydx(k,j,i-1)/(cw(k,j,i-1)+cw(k+1,j,i-1))   &
!!$               * hlf * ( &
!!$               hlf * ( dy(j  ,i-1) + dy(j-1,i-1) )*v(k,j,i-1) + &
!!$               hlf * ( dy(j+1,i-1) + dy(j  ,i-1) )*v(k,j+1,i-1)) )

          do k = 2,nz-1 ! interior levels
             su_integr(j,i) = su_integr(j,i) + Arx(k,j,i)
             uf_integr(j,i) = uf_integr(j,i) + Arx(k,j,i) * u(k,j,i) &
                  - qrt * ( &
                  + zxdy(k,j,i  ) * dzw(k  ,j,i  ) * w(k  ,j,i  ) &
                  + zxdy(k,j,i  ) * dzw(k+1,j,i  ) * w(k+1,j,i  ) &
                  + zxdy(k,j,i-1) * dzw(k  ,j,i-1) * w(k  ,j,i-1) & 
                  + zxdy(k,j,i-1) * dzw(k+1,j,i-1) * w(k+1,j,i-1) &
                  )  ! umask
          enddo

          k = nz ! upper level
          su_integr(j,i) = su_integr(j,i) + Arx(k,j,i)
          uf_integr(j,i) = uf_integr(j,i) + Arx(k,j,i) * u(k,j,i) &
               - qrt * ( &
               + zxdy(k,j,i  )*       dzw(k  ,j,i  )*w(k  ,j,i  ) & 
               + zxdy(k,j,i  )* two * dzw(k+1,j,i  )*w(k+1,j,i  ) & 
               + zxdy(k,j,i-1)*       dzw(k  ,j,i-1)*w(k  ,j,i-1) &
               + zxdy(k,j,i-1)* two * dzw(k+1,j,i-1)*w(k+1,j,i-1) & 
               )  ! umask

       enddo
    enddo

    allocate(uf_integr_tmp(1:nz,0:ny+1,0:nx+1))
    uf_integr_tmp(:,:,:) = zero
    uf_integr_tmp(1,1:ny,1:nx+1)=uf_integr(1:ny,1:nx+1)
    call fill_halo(1,uf_integr_tmp,lbc_null='u')
    uf_integr(0:ny+1,0:nx+1)=uf_integr_tmp(1,0:ny+1,0:nx+1)
    deallocate(uf_integr_tmp)

    if (check_output) then
       call write_netcdf(uf_integr,vname='ufin',netcdf_file_name='coin.nc',rank=myrank,iter=iter_coupling_in)
    endif

    do i = 1,nx+1  
       do j = 1,ny
!       do j = 0,ny+1
          uf_integr(j,i) = uf_integr(j,i) - uf_bar(j,i)
       enddo
    enddo

    do i = 1,nx
!    do i = 0,nx+1
       do j = 1,ny+1

          k = 1 ! lower level
          sv_integr(j,i) = Ary(k,j,i)
          vf_integr(j,i) = Ary(k,j,i) * v(k,j,i) &
               - qrt * ( &
               + zydx(k,j  ,i)* two * dzw(k  ,j  ,i) * w(k  ,j  ,i) & 
               + zydx(k,j  ,i)*       dzw(k+1,j  ,i) * w(k+1,j  ,i) & 
               + zydx(k,j-1,i)* two * dzw(k  ,j-1,i) * w(k  ,j-1,i) & 
               + zydx(k,j-1,i)*       dzw(k+1,j-1,i) * w(k+1,j-1,i) &  
               ) !* vmask(j,i)

!!$               - qrt * ( &
!!$               + zydx(k,j  ,i) * dzw(k+1,j  ,i)*w(k+1-1,j,i) & !note the index trick for w
!!$               + zydx(k,j-1,i) * dzw(k+1,j-1,i)*w(k+1-1,j-1,i) ) & !note the index trick for w        
!!$               -( &
!!$               + zydx(k,j  ,i) * zydx(k,j  ,i)/(cw(k,j  ,i)+cw(k+1,j  ,i)) &
!!$               + zydx(k,j-1,i) * zydx(k,j-1,i)/(cw(k,j-1,i)+cw(k+1,j-1,i)) ) * &
!!$               hlf * ( dy(j,i) + dy(j-1,i) ) * v(k,j,i) &
!!$               - ( &
!!$               + zxdy(k,j  ,i) * zydx(k,j  ,i)/(cw(k,j  ,i)+cw(k+1,j  ,i))  &
!!$               * hlf * ( &
!!$               hlf * ( dx(j,i  ) + dx(j,i-1) ) * u(k,j,i) + &
!!$               hlf * ( dx(j,i+1) + dx(j,i  ) ) * u(k,j,i+1)) &
!!$               + zxdy(k,j-1,i) * zydx(k,j-1,i)/(cw(k,j-1,i)+cw(k+1,j-1,i))   &
!!$               * hlf * ( &
!!$               hlf * ( dx(j-1,i  ) + dx(j-1,i-1) ) * u(k,j-1,i) + &
!!$               hlf * ( dx(j-1,i+1) + dx(j-1,i  ) ) * u(k,j-1,i+1)) &
!!$               ) 

          do k = 2,nz-1 ! interior levels
             sv_integr(j,i) = sv_integr(j,i) + Ary(k,j,i)
             vf_integr(j,i) = vf_integr(j,i) + Ary(k,j,i) * v(k,j,i) &
                  - qrt * ( &
                  + zydx(k,j  ,i) * dzw(k  ,j  ,i) * w(k  ,j  ,i) & 
                  + zydx(k,j  ,i) * dzw(k+1,j  ,i) * w(k+1,j  ,i) & 
                  + zydx(k,j-1,i) * dzw(k  ,j-1,i) * w(k  ,j-1,i) & 
                  + zydx(k,j-1,i) * dzw(k+1,j-1,i) * w(k+1,j-1,i) & 
                  )  !* vmask(j,i)
          enddo

          k = nz ! upper level
          sv_integr(j,i) = sv_integr(j,i) + Ary(k,j,i)
          vf_integr(j,i) = vf_integr(j,i) + Ary(k,j,i) * v(k,j,i) &
               - qrt * ( &
               + zydx(k,j  ,i)*       dzw(k  ,j  ,i) * w(k  ,j  ,i) & 
               + zydx(k,j  ,i)* two * dzw(k+1,j  ,i) * w(k+1,j  ,i) &
               + zydx(k,j-1,i)*       dzw(k  ,j-1,i) * w(k  ,j-1,i) & 
               + zydx(k,j-1,i)* two * dzw(k+1,j-1,i) * w(k+1,j-1,i) & 
               ) !* vmask(j,i)

       enddo
    enddo

    allocate(vf_integr_tmp(1:nz,0:ny+1,0:nx+1))
    vf_integr_tmp(:,:,:) = zero
    vf_integr_tmp(1,1:ny+1,1:nx)=vf_integr(1:ny+1,1:nx)
    call fill_halo(1,vf_integr_tmp,lbc_null='v')
    vf_integr(0:ny+1,0:nx+1)=vf_integr_tmp(1,0:ny+1,0:nx+1)
    deallocate(vf_integr_tmp)

    if (check_output) then
       call write_netcdf(vf_integr,vname='vfin',netcdf_file_name='coin.nc',rank=myrank,iter=iter_coupling_in)
    endif

    do i = 1,nx
!    do i = 0,nx+1
       do j = 1,ny+1
          vf_integr(j,i) = vf_integr(j,i) - vf_bar(j,i)
       enddo
    enddo

    !-------------------------------
    ! correct u,v,w at each depth

    allocate(wc(1:nz+1,0:ny+1,0:nx+1))
    wc(:,:,:) = zero

    do i = 1,nx
       do j = 1,ny
          
          k = 1 ! bottom

          wc(k,j,i) = ( &

               + hlf *( &
                 (zxdy(k,j,i)/dy(j,i)) * &
                 ( uf_integr(j,i  )/su_integr(j,i  ) & 
                  +uf_integr(j,i+1)/su_integr(j,i+1) ) ) &
            
               + hlf *( &
                 (zydx(k,j,i)/dx(j,i)) * &
                 ( vf_integr(j  ,i)/sv_integr(j  ,i) & 
                  +vf_integr(j+1,i)/sv_integr(j+1,i) ) ) &

                      )

          do k = 2,nz ! interior levels

             wc(k,j,i) = ( &

               - (zw(k,j,i)-zw(1,j,i))/(zw(nz+1,j,i)-zw(1,j,i)) &
                *( uf_integr(j,i+1)-uf_integr(j,i)+vf_integr(j+1,i)-vf_integr(j,i) ) / (dx(j,i)*dy(j,i)) &

               + qrt *( &
                 (zxdy(k-1,j,i)/dy(j,i)+zxdy(k,j,i)/dy(j,i)) * &
                 ( uf_integr(j,i  )/su_integr(j,i  ) & 
                  +uf_integr(j,i+1)/su_integr(j,i+1) ) ) &
            
               + qrt *( &
                 (zydx(k-1,j,i)/dx(j,i)+zydx(k,j,i)/dx(j,i)) * &
                 ( vf_integr(j  ,i)/sv_integr(j  ,i) & 
                  +vf_integr(j+1,i)/sv_integr(j+1,i) ) ) &
            
                          )
          enddo

          k = nz+1 ! surface

          wc(k,j,i) = ( &

               - ( uf_integr(j,i+1)-uf_integr(j,i)+vf_integr(j+1,i)-vf_integr(j,i) ) / (dx(j,i)*dy(j,i)) &

               + hlf *( &
                 (zxdy(k-1,j,i)/dy(j,i)) * &
                 ( uf_integr(j,i  )/su_integr(j,i  ) & 
                  +uf_integr(j,i+1)/su_integr(j,i+1) ) ) &

               + hlf *( &
                 (zydx(k-1,j,i)/dx(j,i)) * &
                 ( vf_integr(j  ,i)/sv_integr(j  ,i) & 
                  +vf_integr(j+1,i)/sv_integr(j+1,i) ) ) &

                                 )

       enddo
    enddo

    call fill_halo(1,wc)
    if (check_output) then
       call write_netcdf(wc,vname='wcin',netcdf_file_name='coin.nc',rank=myrank,iter=iter_coupling_in)
    endif

    do i = 1,nx+1  
       do j = 1,ny 
!       do j = 0,ny+1

          k = 1 ! lower level

             u(k,j,i) = u(k,j,i) - ( &
                 Arx(k,j,i)*uf_integr(j,i)/su_integr(j,i) &
               + qrt *( &
                 zxdy(k,j,i  )* two *dzw(k  ,j,i  )*wc(k  ,j,i  ) & 
               + zxdy(k,j,i  )      *dzw(k+1,j,i  )*wc(k+1,j,i  ) &
               + zxdy(k,j,i-1)* two *dzw(k  ,j,i-1)*wc(k  ,j,i-1) & 
               + zxdy(k,j,i-1)      *dzw(k+1,j,i-1)*wc(k+1,j,i-1) ) &        
                                   ) / Arx(k,j,i)
 
          do k = 2,nz-1 ! interior levels

             u(k,j,i) = u(k,j,i) - ( &
                 Arx(k,j,i)*uf_integr(j,i)/su_integr(j,i) &
               + qrt *( &
                 zxdy(k,j,i  ) *dzw(k  ,j,i  )*wc(k  ,j,i  ) & 
               + zxdy(k,j,i  ) *dzw(k+1,j,i  )*wc(k+1,j,i  ) &
               + zxdy(k,j,i-1) *dzw(k  ,j,i-1)*wc(k  ,j,i-1) & 
               + zxdy(k,j,i-1) *dzw(k+1,j,i-1)*wc(k+1,j,i-1) ) & 
                                   ) / Arx(k,j,i)
          enddo

          k = nz !upper level

             u(k,j,i) = u(k,j,i) - ( &
                 Arx(k,j,i)*uf_integr(j,i)/su_integr(j,i) &
               + qrt *( &
                 zxdy(k,j,i  )      *dzw(k  ,j,i  )*wc(k  ,j,i  ) & 
               + zxdy(k,j,i  )* two *dzw(k+1,j,i  )*wc(k+1,j,i  ) &
               + zxdy(k,j,i-1)      *dzw(k  ,j,i-1)*wc(k  ,j,i-1) & 
               + zxdy(k,j,i-1)* two *dzw(k+1,j,i-1)*wc(k+1,j,i-1) ) &        
                                   ) / Arx(k,j,i)

       enddo
    enddo

    allocate(uf_integr_tmp(1:nz,0:ny+1,0:nx+1))
    uf_integr_tmp(:,:,:) = zero
    uf_integr_tmp(:,1:ny,1:nx+1)=u(:,1:ny,1:nx+1)
    call fill_halo(1,uf_integr_tmp,lbc_null='u')
    u(:,0:ny+1,1:nx+1)=uf_integr_tmp(:,0:ny+1,1:nx+1)
    deallocate(uf_integr_tmp)

    do i = 1,nx
!    do i = 0,nx+1
       do j = 1,ny+1

          k = 1 ! lower level

             v(k,j,i) = v(k,j,i) - ( &
                 Ary(k,j,i)*vf_integr(j,i)/sv_integr(j,i) &
               + qrt *( &
                 zydx(k,j  ,i)* two * dzw(k  ,j  ,i)*w(k  ,j  ,i) & 
               + zydx(k,j  ,i)*       dzw(k+1,j  ,i)*w(k+1,j  ,i) & 
               + zydx(k,j-1,i)* two * dzw(k  ,j-1,i)*w(k  ,j-1,i) & 
               + zydx(k,j-1,i)*       dzw(k+1,j-1,i)*w(k+1,j-1,i) ) & 
                                   ) / Ary(k,j,i)

          do k = 2,nz-1 ! interior levels

             v(k,j,i) = v(k,j,i) - ( &
                 Ary(k,j,i)*vf_integr(j,i)/sv_integr(j,i) &
               + qrt *( &
                 zydx(k,j  ,i)* dzw(k  ,j  ,i)*w(k  ,j  ,i) & 
               + zydx(k,j  ,i)* dzw(k+1,j  ,i)*w(k+1,j  ,i) &  
               + zydx(k,j-1,i)* dzw(k  ,j-1,i)*w(k  ,j-1,i) & 
               + zydx(k,j-1,i)* dzw(k+1,j-1,i)*w(k+1,j-1,i) ) &  
                                   ) / Ary(k,j,i)
          
          enddo

          k = nz !upper level

             v(k,j,i) = v(k,j,i) - ( &
                 Ary(k,j,i)*vf_integr(j,i)/sv_integr(j,i) &
               + qrt *( &
                 zydx(k,j  ,i)*       dzw(k  ,j  ,i)*w(k  ,j  ,i) & 
               + zydx(k,j  ,i)* two * dzw(k+1,j  ,i)*w(k+1,j  ,i) & 
               + zydx(k,j-1,i)*       dzw(k  ,j-1,i)*w(k  ,j-1,i) & 
               + zydx(k,j-1,i)* two * dzw(k+1,j-1,i)*w(k+1,j-1,i) ) &  
                                   ) / Ary(k,j,i)

       enddo
    enddo

    allocate(vf_integr_tmp(1:nz,0:ny+1,0:nx+1))
    vf_integr_tmp(:,:,:) = zero
    vf_integr_tmp(:,1:ny+1,1:nx)=v(:,1:ny+1,1:nx)
    call fill_halo(1,vf_integr_tmp,lbc_null='v')
    v(:,1:ny+1,0:nx+1)=vf_integr_tmp(:,1:ny+1,0:nx+1)
    deallocate(vf_integr_tmp)

    do i = 1,nx
       do j = 1,ny
          
          do k = 1,nz+1 ! all levels

             w(k,j,i) = w(k,j,i) - wc(k,j,i) 

          enddo

       enddo
    enddo

    call fill_halo(1,w)

    deallocate(su_integr)
    deallocate(sv_integr)
    deallocate(uf_integr)
    deallocate(vf_integr)
    deallocate(wc)

    if ((present(uf)).and.(present(vf))) then
    !----------------------------------------------
    ! compute horizontal transport at each depth

    do i = 1,nx+1  
       do j = 1,ny 
!       do j = 0,ny+1

          k = 1 ! lower level
          uf(k,j,i) = Arx(k,j,i) * u(k,j,i) &
               - qrt * ( &
               + zxdy(k,j,i  )* two * dzw(k  ,j,i  )*w(k  ,j,i  ) & 
               + zxdy(k,j,i  )*       dzw(k+1,j,i  )*w(k+1,j,i  ) & 
               + zxdy(k,j,i-1)* two * dzw(k  ,j,i-1)*w(k  ,j,i-1) &
               + zxdy(k,j,i-1)*       dzw(k+1,j,i-1)*w(k+1,j,i-1) & 
               )  ! umask

!!$               - qrt * ( &
!!$               + zxdy(k,j,i  )*dzw(k+1,j,i  )*w(k+1-1,j,i)   & !note the index trick for wc
!!$               + zxdy(k,j,i-1)*dzw(k+1,j,i-1)*w(k+1-1,j,i-1) ) & !note the index trick for wc
!!$               -( &
!!$               + zxdy(k,j,i  )*zxdy(k,j,i  )/(cw(k,j,i  )+cw(k+1,j,i  )) &
!!$               + zxdy(k,j,i-1)*zxdy(k,j,i-1)/(cw(k,j,i-1)+cw(k+1,j,i-1)) ) * &
!!$               (hlf * ( dx(j,i) + dx(j,i-1) )) * u(k,j,i) &
!!$               -( & 
!!$               + zxdy(k,j,i) * zydx(k,j,i)/(cw(k,j,i  )+cw(k+1,j,i  ))   &
!!$               * hlf * ( &
!!$               hlf * ( dy(j  ,i) + dy(j-1,i) )*v(k,j,i ) + &
!!$               hlf * ( dy(j+1,i) + dy(j  ,i) )*v(k,j+1,i ))   & 
!!$               + zxdy(k,j,i-1) * zydx(k,j,i-1)/(cw(k,j,i-1)+cw(k+1,j,i-1))   &
!!$               * hlf * ( &
!!$               hlf * ( dy(j  ,i-1) + dy(j-1,i-1) )*v(k,j,i-1) + &
!!$               hlf * ( dy(j+1,i-1) + dy(j  ,i-1) )*v(k,j+1,i-1)) )

          do k = 2,nz-1 ! interior levels            
             uf(k,j,i) = Arx(k,j,i) * u(k,j,i) &
                  - qrt * ( &
                  + zxdy(k,j,i  ) * dzw(k  ,j,i  )*w(k  ,j,i  ) & 
                  + zxdy(k,j,i  ) * dzw(k+1,j,i  )*w(k+1,j,i  ) & 
                  + zxdy(k,j,i-1) * dzw(k  ,j,i-1)*w(k  ,j,i-1) &
                  + zxdy(k,j,i-1) * dzw(k+1,j,i-1)*w(k+1,j,i-1) &
                  )  ! umask
          enddo

          k = nz ! upper level
          uf(k,j,i) = Arx(k,j,i) * u(k,j,i) &
               - qrt * ( &
               + zxdy(k,j,i  )*       dzw(k  ,j,i  )*w(k  ,j,i  ) & 
               + zxdy(k,j,i  )* two * dzw(k+1,j,i  )*w(k+1,j,i  ) & 
               + zxdy(k,j,i-1)*       dzw(k  ,j,i-1)*w(k  ,j,i-1) &
               + zxdy(k,j,i-1)* two * dzw(k+1,j,i-1)*w(k+1,j,i-1) & 
               )  ! umask

       enddo
    enddo

    allocate(uf_tmp(1:nz,0:ny+1,0:nx+1))
    uf_tmp(:,1:ny,1:nx+1)=uf(:,1:ny,1:nx+1)
    call fill_halo(1,uf_tmp,lbc_null='u')
    uf(:,0:ny+1,1:nx+1)=uf_tmp(:,0:ny+1,1:nx+1)
    deallocate(uf_tmp)

    do i = 1,nx
!    do i = 0,nx+1
       do j = 1,ny+1

          k = 1 ! lower level
          vf(k,j,i) = Ary(k,j,i) * v(k,j,i) &
               - qrt * ( &
               + zydx(k,j  ,i)* two * dzw(k  ,j  ,i)*w(k  ,j  ,i) & 
               + zydx(k,j  ,i)*       dzw(k+1,j  ,i)*w(k+1,j  ,i) &  
               + zydx(k,j-1,i)* two * dzw(k  ,j-1,i)*w(k  ,j-1,i) & 
               + zydx(k,j-1,i)*       dzw(k+1,j-1,i)*w(k+1,j-1,i) & 
               ) !* vmask(j,i)

!!$               - qrt * ( &
!!$               + zydx(k,j  ,i) * dzw(k+1,j  ,i)*w(k+1-1,j,i) & !note the index trick for wc
!!$               + zydx(k,j-1,i) * dzw(k+1,j-1,i)*w(k+1-1,j-1,i) ) & !note the index trick for wc        
!!$               -( &
!!$               + zydx(k,j  ,i) * zydx(k,j  ,i)/(cw(k,j  ,i)+cw(k+1,j  ,i)) &
!!$               + zydx(k,j-1,i) * zydx(k,j-1,i)/(cw(k,j-1,i)+cw(k+1,j-1,i)) ) * &
!!$               hlf * ( dy(j,i) + dy(j-1,i) ) * v(k,j,i) &
!!$               - ( &
!!$               + zxdy(k,j  ,i) * zydx(k,j  ,i)/(cw(k,j  ,i)+cw(k+1,j  ,i))  &
!!$               * hlf * ( &
!!$               hlf * ( dx(j,i  ) + dx(j,i-1) ) * u(k,j,i) + &
!!$               hlf * ( dx(j,i+1) + dx(j,i  ) ) * u(k,j,i+1)) &
!!$               + zxdy(k,j-1,i) * zydx(k,j-1,i)/(cw(k,j-1,i)+cw(k+1,j-1,i))   &
!!$               * hlf * ( &
!!$               hlf * ( dx(j-1,i  ) + dx(j-1,i-1) ) * u(k,j-1,i) + &
!!$               hlf * ( dx(j-1,i+1) + dx(j-1,i  ) ) * u(k,j-1,i+1)) &
!!$               ) 

          do k = 2,nz-1 ! interior levels
             vf(k,j,i) = Ary(k,j,i) * v(k,j,i) &
                  - qrt * ( &
                  + zydx(k,j  ,i) * dzw(k  ,j  ,i)*w(k  ,j  ,i) &
                  + zydx(k,j  ,i) * dzw(k+1,j  ,i)*w(k+1,j  ,i) &
                  + zydx(k,j-1,i) * dzw(k  ,j-1,i)*w(k  ,j-1,i) &
                  + zydx(k,j-1,i) * dzw(k+1,j-1,i)*w(k+1,j-1,i) &
                  )  !* vmask(j,i)
          enddo

          k = nz ! upper level
          vf(k,j,i) = Ary(k,j,i) * v(k,j,i) &
               - qrt * ( &
               + zydx(k,j  ,i)*       dzw(k  ,j  ,i)*w(k  ,j  ,i) & 
               + zydx(k,j  ,i)* two * dzw(k+1,j  ,i)*w(k+1,j  ,i) & 
               + zydx(k,j-1,i)*       dzw(k  ,j-1,i)*w(k  ,j-1,i) & 
               + zydx(k,j-1,i)* two * dzw(k+1,j-1,i)*w(k+1,j-1,i) & 
               ) !* vmask(j,i)

       enddo
    enddo

    allocate(vf_tmp(1:nz,0:ny+1,0:nx+1))
    vf_tmp(:,1:ny+1,1:nx)=vf(:,1:ny+1,1:nx)
    call fill_halo(1,vf_tmp,lbc_null='v')
    vf(:,1:ny+1,0:nx+1)=vf_tmp(:,1:ny+1,0:nx+1)
    deallocate(vf_tmp)

    !! check
    allocate(uf_integr(0:ny+1,1:nx+1))
    allocate(vf_integr(1:ny+1,0:nx+1))
    uf_integr(:,:) = zero
    vf_integr(:,:) = zero
    do i = 1,nx+1  
       do j = 1,ny 
!       do j = 0,ny+1
          do k = 1,nz
             uf_integr(j,i) = uf_integr(j,i) + uf(k,j,i)
          enddo
       enddo
    enddo
    do i = 1,nx
!    do i = 0,nx+1
       do j = 1,ny+1
          do k = 1,nz
             vf_integr(j,i) = vf_integr(j,i) + vf(k,j,i)
          enddo
       enddo
    enddo
    if (check_output) then
       call write_netcdf(uf_integr,vname='ufout',netcdf_file_name='coin.nc',rank=myrank,iter=iter_coupling_in)
       call write_netcdf(vf_integr,vname='vfout',netcdf_file_name='coin.nc',rank=myrank,iter=iter_coupling_in)
       call write_netcdf(uf_integr-uf_bar,vname='diff_uf',netcdf_file_name='coin.nc',rank=myrank,iter=iter_coupling_in)
       call write_netcdf(vf_integr-vf_bar,vname='diff_vf',netcdf_file_name='coin.nc',rank=myrank,iter=iter_coupling_in)
    endif
    deallocate(uf_integr)
    deallocate(vf_integr)
    !! check

    endif

  end subroutine bt2bc_coupling

end module mg_coupling
