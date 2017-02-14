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
  subroutine bc2bt_coupling(ru,rv,rw,Hz,Hzh,dt,rufrc,rvfrc)

    real(kind=rp), dimension(:,:,:), pointer, intent(in)  :: ru,rv,rw
    real(kind=rp), dimension(:,:,:), pointer, intent(in)  :: Hz,Hzh
    real(kind=rp),                            intent(in)  :: dt
    real(kind=rp), dimension(:,:)  , pointer, intent(out) :: rufrc,rvfrc

    integer(kind=ip):: k, j, i
    integer(kind=ip):: nx, ny, nz

    real(kind=rp) :: gamma
    real(kind=rp), dimension(:,:)  , pointer :: dx,dy,Arz
    real(kind=rp), dimension(:,:,:), pointer :: dz
    real(kind=rp), dimension(:,:,:), pointer :: Arx,Ary
    real(kind=rp), dimension(:,:,:), pointer :: zxdy,zydx
    real(kind=rp), dimension(:,:,:), pointer :: alpha
    real(kind=rp), dimension(:,:)  , pointer :: beta
    real(kind=rp), dimension(:,:,:), pointer :: p

    real(kind=rp), dimension(:,:,:), pointer :: rutmp,rvtmp,rwtmp

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz

    dx    => grid(1)%dx
    dy    => grid(1)%dy
    dz    => grid(1)%dz
    alpha => grid(1)%alpha
    beta  => grid(1)%beta
    Arx   => grid(1)%Arx
    Ary   => grid(1)%Ary
    Arz   => grid(1)%Arz
    zxdy  => grid(1)%zxdy
    zydx  => grid(1)%zydx
    p     => grid(1)%p

    !------------------------------------------
    ! add the volume integrated grad of nh pressure

    allocate(rutmp(1:nz,0:ny+1,1:nx+1))
    allocate(rvtmp(1:nz,1:ny+1,0:nx+1))
    allocate(rwtmp(1:nz+1,0:ny+1,0:nx+1))

    do i = 1,nx+1
       do j = 0,ny+1
          do k = 1,nz
             rutmp(k,j,i) = ru(k,j,i) - Arx(k,j,i) * (p(k,j,i  )-p(k,j,i-1)) / dt !&
!                                      *(Hz(k,j,i)+Hz(k,j,i-1)) &
!                                      /(Hzh(k,j,i)+Hzh(k,j,i-1))
          enddo
       enddo
    enddo

    do i = 0,nx+1
       do j = 1,ny+1
          do k = 1,nz
             rvtmp(k,j,i) = rv(k,j,i) - Ary(k,j,i) * (p(k,j  ,i)-p(k,j-1,i)) / dt !&
!                                      *(Hz(k,j,i)+Hz(k,j-1,i)) &
!                                      /(Hzh(k,j,i)+Hzh(k,j-1,i))
          enddo
       enddo
    enddo

    do i = 0,nx+1
       do j = 0,ny+1
         k = 1
         !do nothing
         do k = 2,nz
             rwtmp(k,j,i) = rw(k,j,i) - Arz(j,i) * (p(k  ,j,i)-p(k-1,j,i)) / dt !&
!                                      *(Hz(k,j,i)+Hz(k-1,j,i)) &
!                                      /(Hzh(k,j,i)+Hzh(k-1,j,i))
         enddo
         k = nz+1
         rwtmp(k,j,i) = rw(k,j,i) - Arz(j,i) * (        -p(k-1,j,i)) / dt !&
!                                     *Hz(k,j,i) &
!                                     /Hzh(k,j,i)
       enddo
    enddo

    !------------------------------------------
    ! fluxes and integrate

    do i = 1,nx+1
       do j = 1,ny

          k = 1
          gamma = one &
               - qrt *( &
               (zxdy(k,j,i  )/dy(j,i  ))**2/alpha(k,j,i  ) + &
               (zxdy(k,j,i-1)/dy(j,i-1))**2/alpha(k,j,i-1) )
          rufrc(j,i) = gamma * rutmp(k,j,i) &
               - qrt * ( &
               + zxdy(k,j,i  ) /dy(j,i  )   * rwtmp(k+1,j,i  ) &
               + zxdy(k,j,i-1) /dy(j,i-1)   * rwtmp(k+1,j,i-1) ) &
               - beta(j,i-1)   /dz(k,j,i-1) * rvtmp(k,j  ,i-1) &
               - beta(j,i-1)   /dz(k,j,i-1) * rvtmp(k,j+1,i-1) &
               - beta(j,i  )   /dz(k,j,i  ) * rvtmp(k,j  ,i  ) &
               - beta(j,i  )   /dz(k,j,i  ) * rvtmp(k,j+1,i  )

          do k = 2,nz-1
             rufrc(j,i) = rufrc(j,i) &
                  + rutmp(k,j,i) &
                  - qrt * ( &
                  + zxdy(k,j,i  ) /dy(j,i  ) * rwtmp(k  ,j,i  ) &
                  + zxdy(k,j,i  ) /dy(j,i  ) * rwtmp(k+1,j,i  ) &
                  + zxdy(k,j,i-1) /dy(j,i-1) * rwtmp(k  ,j,i-1) & 
                  + zxdy(k,j,i-1) /dy(j,i-1) * rwtmp(k+1,j,i-1) )
          enddo

          k = nz
          rufrc(j,i) = rufrc(j,i) &
               + rutmp(k,j,i) &
               - qrt * ( &
               + zxdy(k,j,i  )     /dy(j,i  ) * rwtmp(k  ,j,i  ) & 
               + zxdy(k,j,i  ) *two/dy(j,i  ) * rwtmp(k+1,j,i  ) & 
               + zxdy(k,j,i-1)     /dy(j,i-1) * rwtmp(k  ,j,i-1) &
               + zxdy(k,j,i-1) *two/dy(j,i-1) * rwtmp(k+1,j,i-1) )
          
       enddo
    enddo

    do i = 1,nx
       do j = 1,ny+1 

          k = 1
          gamma = one &
               - qrt *( &
               (zydx(k,j  ,i)/dx(j  ,i))**2/alpha(k,j,i  ) + &
               (zydx(k,j-1,i)/dx(j-1,i))**2/alpha(k,j-1,i) )
          rvfrc(j,i) = gamma * rvtmp(k,j,i) &
               - qrt * ( &
               + zydx(k,j  ,i) /dx(j  ,i)   * rwtmp(k+1,j  ,i  ) &
               + zydx(k,j-1,i) /dx(j-1,i)   * rwtmp(k+1,j-1,i  ) ) &
               - beta(j-1,i)   /dz(k,j-1,i) * rutmp(k  ,j-1,i  ) &
               - beta(j-1,i)   /dz(k,j-1,i) * rutmp(k  ,j-1,i+1) &
               - beta(j  ,i)   /dz(k,j  ,i) * rutmp(k  ,j  ,i  ) &
               - beta(j  ,i)   /dz(k,j  ,i) * rutmp(k  ,j  ,i+1)  

          do k = 2,nz-1
             rvfrc(j,i) = rvfrc(j,i) &
                  + rvtmp(k,j,i) &
                  - qrt * ( &
                  + zydx(k,j  ,i) /dx(j  ,i) * rwtmp(k  ,j  ,i) & 
                  + zydx(k,j  ,i) /dx(j  ,i) * rwtmp(k+1,j  ,i) &
                  + zydx(k,j-1,i) /dx(j-1,i) * rwtmp(k  ,j-1,i) & 
                  + zydx(k,j-1,i) /dx(j-1,i) * rwtmp(k+1,j-1,i) )
          enddo

          k = nz
          rvfrc(j,i) = rvfrc(j,i) &
               + rvtmp(k,j,i) &
               - qrt * ( &
               + zydx(k,j  ,i)     /dx(j  ,i) * rwtmp(k  ,j  ,i) & 
               + zydx(k,j  ,i) *two/dx(j  ,i) * rwtmp(k+1,j  ,i) &
               + zydx(k,j-1,i)     /dx(j-1,i) * rwtmp(k  ,j-1,i) & 
               + zydx(k,j-1,i) *two/dx(j-1,i) * rwtmp(k+1,j-1,i) )

       enddo
    enddo

    deallocate(rutmp)
    deallocate(rvtmp)
    deallocate(rwtmp)

  end subroutine bc2bt_coupling

  !-------------------------------------------------------------------------     
  subroutine bt2bc_coupling(uf_bar,vf_bar,u,v,w,uf,vf)

    real(kind=rp), dimension(:,:),   pointer, intent(in)     :: uf_bar,vf_bar
    real(kind=rp), dimension(:,:,:), pointer, intent(inout)  :: u,v,w
    real(kind=rp), dimension(:,:,:), pointer, optional, intent(out):: uf,vf

    integer(kind=ip):: k,j,i
    integer(kind=ip):: nx,ny,nz

    real(kind=rp) :: Hci

    real(kind=rp), dimension(:,:)  , pointer :: dx,dy
    real(kind=rp), dimension(:,:,:), pointer :: dz
    real(kind=rp), dimension(:,:,:), pointer :: dzw
    real(kind=rp), dimension(:,:,:), pointer :: Arx,Ary
    real(kind=rp), dimension(:,:,:), pointer :: zxdy,zydx
    real(kind=rp), dimension(:,:,:), pointer :: alpha

    real(kind=rp), dimension(:,:)  , pointer :: su_integr,sv_integr
    real(kind=rp), dimension(:,:)  , pointer :: uf_integr,vf_integr
    real(kind=rp), dimension(:,:,:), pointer :: wc

    integer(kind=ip), save :: iter_coupling_in=0
    iter_coupling_in = iter_coupling_in + 1

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz

    dx    => grid(1)%dx
    dy    => grid(1)%dy
    dz    => grid(1)%dz
    dzw   => grid(1)%dzw
    Arx   => grid(1)%Arx
    Ary   => grid(1)%Ary
    alpha => grid(1)%alpha
    zxdy  => grid(1)%zxdy
    zydx  => grid(1)%zydx

    !------------------------------------------
    ! compute integrated transport anomalies

    allocate(su_integr(0:ny+1,1:nx+1))
    allocate(sv_integr(1:ny+1,0:nx+1))
    allocate(uf_integr(0:ny+1,1:nx+1))
    allocate(vf_integr(1:ny+1,0:nx+1))

    do i = 1,nx+1  
       do j = 0,ny+1

          k = 1
          su_integr(j,i) =  Arx(k,j,i)
          uf_integr(j,i) =  Arx(k,j,i) * u(k,j,i) &
               - qrt * ( &
               + zxdy(k,j,i  ) * two * dzw(k  ,j,i  )*w(k  ,j,i  ) & 
               + zxdy(k,j,i  ) *       dzw(k+1,j,i  )*w(k+1,j,i  ) & 
               + zxdy(k,j,i-1) * two * dzw(k  ,j,i-1)*w(k  ,j,i-1) &
               + zxdy(k,j,i-1) *       dzw(k+1,j,i-1)*w(k+1,j,i-1) )

          do k = 2,nz-1
             su_integr(j,i) = su_integr(j,i) + Arx(k,j,i)
             uf_integr(j,i) = uf_integr(j,i) + Arx(k,j,i) * u(k,j,i) &
                  - qrt * ( &
                  + zxdy(k,j,i  ) * dzw(k  ,j,i  ) * w(k  ,j,i  ) &
                  + zxdy(k,j,i  ) * dzw(k+1,j,i  ) * w(k+1,j,i  ) &
                  + zxdy(k,j,i-1) * dzw(k  ,j,i-1) * w(k  ,j,i-1) & 
                  + zxdy(k,j,i-1) * dzw(k+1,j,i-1) * w(k+1,j,i-1) )
          enddo

          k = nz
          su_integr(j,i) = su_integr(j,i) + Arx(k,j,i)
          uf_integr(j,i) = uf_integr(j,i) + Arx(k,j,i) * u(k,j,i) &
               - qrt * ( &
               + zxdy(k,j,i  ) *       dzw(k  ,j,i  )*w(k  ,j,i  ) & 
               + zxdy(k,j,i  ) * two * dzw(k+1,j,i  )*w(k+1,j,i  ) & 
               + zxdy(k,j,i-1) *       dzw(k  ,j,i-1)*w(k  ,j,i-1) &
               + zxdy(k,j,i-1) * two * dzw(k+1,j,i-1)*w(k+1,j,i-1) ) 

       enddo
    enddo

    do i = 0,nx+1
       do j = 1,ny+1

          k = 1
          sv_integr(j,i) = Ary(k,j,i)
          vf_integr(j,i) = Ary(k,j,i) * v(k,j,i) &
               - qrt * ( &
               + zydx(k,j  ,i) * two * dzw(k  ,j  ,i) * w(k  ,j  ,i) & 
               + zydx(k,j  ,i) *       dzw(k+1,j  ,i) * w(k+1,j  ,i) & 
               + zydx(k,j-1,i) * two * dzw(k  ,j-1,i) * w(k  ,j-1,i) & 
               + zydx(k,j-1,i) *       dzw(k+1,j-1,i) * w(k+1,j-1,i) )

          do k = 2,nz-1
             sv_integr(j,i) = sv_integr(j,i) + Ary(k,j,i)
             vf_integr(j,i) = vf_integr(j,i) + Ary(k,j,i) * v(k,j,i) &
                  - qrt * ( &
                  + zydx(k,j  ,i) * dzw(k  ,j  ,i) * w(k  ,j  ,i) & 
                  + zydx(k,j  ,i) * dzw(k+1,j  ,i) * w(k+1,j  ,i) & 
                  + zydx(k,j-1,i) * dzw(k  ,j-1,i) * w(k  ,j-1,i) & 
                  + zydx(k,j-1,i) * dzw(k+1,j-1,i) * w(k+1,j-1,i) )
          enddo

          k = nz
          sv_integr(j,i) = sv_integr(j,i) + Ary(k,j,i)
          vf_integr(j,i) = vf_integr(j,i) + Ary(k,j,i) * v(k,j,i) &
               - qrt * ( &
               + zydx(k,j  ,i) *       dzw(k  ,j  ,i) * w(k  ,j  ,i) & 
               + zydx(k,j  ,i) * two * dzw(k+1,j  ,i) * w(k+1,j  ,i) &
               + zydx(k,j-1,i) *       dzw(k  ,j-1,i) * w(k  ,j-1,i) & 
               + zydx(k,j-1,i) * two * dzw(k+1,j-1,i) * w(k+1,j-1,i) )

       enddo
    enddo

    !if (check_output) then
    !   if ((iter_coupling_in .EQ. 1) .OR. (iter_coupling_in .EQ. 2)) then
    !      call write_netcdf(uf_integr,vname='ufin',netcdf_file_name='coin.nc',rank=myrank,iter=iter_coupling_in)
    !      call write_netcdf(vf_integr,vname='vfin',netcdf_file_name='coin.nc',rank=myrank,iter=iter_coupling_in)
    !   endif
    !endif

    do i = 1,nx+1  
       do j = 0,ny+1
          uf_integr(j,i) = uf_integr(j,i) - uf_bar(j,i)
       enddo
    enddo

    do i = 0,nx+1
       do j = 1,ny+1
          vf_integr(j,i) = vf_integr(j,i) - vf_bar(j,i)
       enddo
    enddo

    !if (check_output) then
    !   if ((iter_coupling_in .EQ. 1) .OR. (iter_coupling_in .EQ. 2)) then
    !      !if ((iter_coupling_in .EQ. 199) .OR. (iter_coupling_in .EQ. 200)) then
    !      call write_netcdf(uf_integr,vname='diff_ufin',netcdf_file_name='coin.nc',rank=myrank,iter=iter_coupling_in)
    !      call write_netcdf(vf_integr,vname='diff_vfin',netcdf_file_name='coin.nc',rank=myrank,iter=iter_coupling_in)
    !   endif
    !endif

    !-------------------------------
    ! evaluate correction for w

!ND 11/01 : should rather be computed over 0,nx+1 0:ny+1, without the final fill_halo

    allocate(wc(1:nz+1,0:ny+1,0:nx+1))

    do i = 1,nx
       do j = 1,ny

          k = 1 
          wc(k,j,i) = ( &
               + hlf *( &
               (zxdy(k,j,i)/dy(j,i)) * &
               (uf_integr(j,i  )/su_integr(j,i  ) & 
               +uf_integr(j,i+1)/su_integr(j,i+1) ) ) &
               + hlf *( &
               (zydx(k,j,i)/dx(j,i)) * &
               (vf_integr(j  ,i)/sv_integr(j  ,i) & 
               +vf_integr(j+1,i)/sv_integr(j+1,i) ) ) &
               )

          Hci= one / sum(dz(:,j,i))

          do k = 2,nz 
             wc(k,j,i) = ( &
                  - sum(dz(1:k-1,j,i)) * Hci * &
                  ( uf_integr(j,i+1)-uf_integr(j,i)+vf_integr(j+1,i)-vf_integr(j,i) ) / (dx(j,i)*dy(j,i)) &
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

          k = nz+1
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
!ND 11/01

    !if (check_output) then
    !   if ((iter_coupling_in .EQ. 1) .OR. (iter_coupling_in .EQ. 2)) then
    !      call write_netcdf(wc,vname='wc',netcdf_file_name='coin.nc',rank=myrank,iter=iter_coupling_in)
    !   endif
    !endif

    !-------------------------------
    ! correct u,v,w at each depth

    do i = 1,nx+1  
       do j = 0,ny+1

          k = 1
          u(k,j,i) = u(k,j,i) - ( &
               + uf_integr(j,i)/su_integr(j,i) &
               + qrt / Arx(k,j,i) *( &
               + zxdy(k,j,i  )* two *dzw(k  ,j,i  )*wc(k  ,j,i  ) & 
               + zxdy(k,j,i  )      *dzw(k+1,j,i  )*wc(k+1,j,i  ) &
               + zxdy(k,j,i-1)* two *dzw(k  ,j,i-1)*wc(k  ,j,i-1) & 
               + zxdy(k,j,i-1)      *dzw(k+1,j,i-1)*wc(k+1,j,i-1) ) &
               )                  

          do k = 2,nz-1
             u(k,j,i) = u(k,j,i) - ( &
                  + uf_integr(j,i)/su_integr(j,i) &
                  + qrt / Arx(k,j,i) *( &
                  + zxdy(k,j,i  ) *dzw(k  ,j,i  )*wc(k  ,j,i  ) & 
                  + zxdy(k,j,i  ) *dzw(k+1,j,i  )*wc(k+1,j,i  ) &
                  + zxdy(k,j,i-1) *dzw(k  ,j,i-1)*wc(k  ,j,i-1) & 
                  + zxdy(k,j,i-1) *dzw(k+1,j,i-1)*wc(k+1,j,i-1) ) &
                  )
          enddo

          k = nz
          u(k,j,i) = u(k,j,i) - ( &
               + uf_integr(j,i)/su_integr(j,i) &
               + qrt / Arx(k,j,i) *( &
               + zxdy(k,j,i  )      *dzw(k  ,j,i  )*wc(k  ,j,i  ) & 
               + zxdy(k,j,i  )* two *dzw(k+1,j,i  )*wc(k+1,j,i  ) &
               + zxdy(k,j,i-1)      *dzw(k  ,j,i-1)*wc(k  ,j,i-1) & 
               + zxdy(k,j,i-1)* two *dzw(k+1,j,i-1)*wc(k+1,j,i-1) ) &
               )

       enddo
    enddo

    do i = 0,nx+1
       do j = 1,ny+1

          k = 1
          v(k,j,i) = v(k,j,i) - ( &
               + vf_integr(j,i)/sv_integr(j,i) &
               + qrt / Ary(k,j,i) *( &
               + zydx(k,j  ,i)* two * dzw(k  ,j  ,i)*wc(k  ,j  ,i) & 
               + zydx(k,j  ,i)*       dzw(k+1,j  ,i)*wc(k+1,j  ,i) & 
               + zydx(k,j-1,i)* two * dzw(k  ,j-1,i)*wc(k  ,j-1,i) & 
               + zydx(k,j-1,i)*       dzw(k+1,j-1,i)*wc(k+1,j-1,i) ) &
               )

          do k = 2,nz-1
             v(k,j,i) = v(k,j,i) - ( &
                  + vf_integr(j,i)/sv_integr(j,i) &
                  + qrt / Ary(k,j,i) *( &
                  + zydx(k,j  ,i)* dzw(k  ,j  ,i)*wc(k  ,j  ,i) & 
                  + zydx(k,j  ,i)* dzw(k+1,j  ,i)*wc(k+1,j  ,i) &  
                  + zydx(k,j-1,i)* dzw(k  ,j-1,i)*wc(k  ,j-1,i) & 
                  + zydx(k,j-1,i)* dzw(k+1,j-1,i)*wc(k+1,j-1,i) ) &
                  )
          enddo

          k = nz
          v(k,j,i) = v(k,j,i) - ( &
               + vf_integr(j,i)/sv_integr(j,i) &
               + qrt / Ary(k,j,i) *( &
               + zydx(k,j  ,i)*       dzw(k  ,j  ,i)*wc(k  ,j  ,i) & 
               + zydx(k,j  ,i)* two * dzw(k+1,j  ,i)*wc(k+1,j  ,i) & 
               + zydx(k,j-1,i)*       dzw(k  ,j-1,i)*wc(k  ,j-1,i) & 
               + zydx(k,j-1,i)* two * dzw(k+1,j-1,i)*wc(k+1,j-1,i) ) &
               )

       enddo
    enddo

    do i = 0,nx+1
       do j = 0,ny+1
          do k = 1,nz+1
             w(k,j,i) = w(k,j,i) - wc(k,j,i)
          enddo
       enddo
    enddo

!    call fill_halo(1,w)

    deallocate(su_integr)
    deallocate(sv_integr)
    deallocate(uf_integr)
    deallocate(vf_integr)
    deallocate(wc)

    !----------------------------------------------
    ! compute horizontal transport at each depth
    if ((present(uf)).and.(present(vf))) then

       do i = 1,nx+1  
          do j = 0,ny+1

             k = 1
             uf(k,j,i) = Arx(k,j,i) * u(k,j,i) &
                  - qrt * ( &
                  + zxdy(k,j,i  ) * two * dzw(k  ,j,i  )*w(k  ,j,i  ) & 
                  + zxdy(k,j,i  ) *       dzw(k+1,j,i  )*w(k+1,j,i  ) & 
                  + zxdy(k,j,i-1) * two * dzw(k  ,j,i-1)*w(k  ,j,i-1) &
                  + zxdy(k,j,i-1) *       dzw(k+1,j,i-1)*w(k+1,j,i-1) ) 

             do k = 2,nz-1
                uf(k,j,i) = Arx(k,j,i) * u(k,j,i) &
                     - qrt * ( &
                     + zxdy(k,j,i  ) * dzw(k  ,j,i  )*w(k  ,j,i  ) & 
                     + zxdy(k,j,i  ) * dzw(k+1,j,i  )*w(k+1,j,i  ) & 
                     + zxdy(k,j,i-1) * dzw(k  ,j,i-1)*w(k  ,j,i-1) &
                     + zxdy(k,j,i-1) * dzw(k+1,j,i-1)*w(k+1,j,i-1) )
             enddo

             k = nz
             uf(k,j,i) = Arx(k,j,i) * u(k,j,i) &
                  - qrt * ( &
                  + zxdy(k,j,i  ) *       dzw(k  ,j,i  )*w(k  ,j,i  ) & 
                  + zxdy(k,j,i  ) * two * dzw(k+1,j,i  )*w(k+1,j,i  ) & 
                  + zxdy(k,j,i-1) *       dzw(k  ,j,i-1)*w(k  ,j,i-1) &
                  + zxdy(k,j,i-1) * two * dzw(k+1,j,i-1)*w(k+1,j,i-1) )

          enddo
       enddo

       do i = 0,nx+1
          do j = 1,ny+1

             k = 1
             vf(k,j,i) = Ary(k,j,i) * v(k,j,i) &
                  - qrt * ( &
                  + zydx(k,j  ,i) * two * dzw(k  ,j  ,i)*w(k  ,j  ,i) & 
                  + zydx(k,j  ,i) *       dzw(k+1,j  ,i)*w(k+1,j  ,i) &  
                  + zydx(k,j-1,i) * two * dzw(k  ,j-1,i)*w(k  ,j-1,i) & 
                  + zydx(k,j-1,i) *       dzw(k+1,j-1,i)*w(k+1,j-1,i) ) 

             do k = 2,nz-1
                vf(k,j,i) = Ary(k,j,i) * v(k,j,i) &
                     - qrt * ( &
                     + zydx(k,j  ,i) * dzw(k  ,j  ,i)*w(k  ,j  ,i) &
                     + zydx(k,j  ,i) * dzw(k+1,j  ,i)*w(k+1,j  ,i) &
                     + zydx(k,j-1,i) * dzw(k  ,j-1,i)*w(k  ,j-1,i) &
                     + zydx(k,j-1,i) * dzw(k+1,j-1,i)*w(k+1,j-1,i) ) 
             enddo

             k = nz
             vf(k,j,i) = Ary(k,j,i) * v(k,j,i) &
                  - qrt * ( &
                  + zydx(k,j  ,i) *       dzw(k  ,j  ,i)*w(k  ,j  ,i) & 
                  + zydx(k,j  ,i) * two * dzw(k+1,j  ,i)*w(k+1,j  ,i) & 
                  + zydx(k,j-1,i) *       dzw(k  ,j-1,i)*w(k  ,j-1,i) & 
                  + zydx(k,j-1,i) * two * dzw(k+1,j-1,i)*w(k+1,j-1,i) ) 

          enddo
       enddo

       !!! check
       !allocate(uf_integr(0:ny+1,0:nx+1))
       !allocate(vf_integr(0:ny+1,0:nx+1))
       !uf_integr(:,:) = zero
       !vf_integr(:,:) = zero
       !do i = 1,nx+1  
       !   do j = 1,ny 
       !      !       do j = 0,ny+1
       !      do k = 1,nz
       !         uf_integr(j,i) = uf_integr(j,i) + uf(k,j,i)
       !      enddo
       !   enddo
       !enddo
       !do i = 1,nx
       !   !    do i = 0,nx+1
       !   do j = 1,ny+1
       !      do k = 1,nz
       !         vf_integr(j,i) = vf_integr(j,i) + vf(k,j,i)
       !      enddo
       !   enddo
       !enddo
       !if (check_output) then
       !   if ((iter_coupling_in .EQ. 1) .OR. (iter_coupling_in .EQ. 2)) then
       !      !if ((iter_coupling_in .EQ. 199) .OR. (iter_coupling_in .EQ. 200)) then
       !      call write_netcdf(uf_integr,vname='ufout',netcdf_file_name='coin.nc',rank=myrank,iter=iter_coupling_in)
       !      call write_netcdf(vf_integr,vname='vfout',netcdf_file_name='coin.nc',rank=myrank,iter=iter_coupling_in)
       !      call write_netcdf(uf_integr-uf_bar,vname='diff_ufout',netcdf_file_name='coin.nc',rank=myrank,iter=iter_coupling_in)
       !      call write_netcdf(vf_integr-vf_bar,vname='diff_vfout',netcdf_file_name='coin.nc',rank=myrank,iter=iter_coupling_in)
       !   endif
       !endif
       !deallocate(uf_integr)
       !deallocate(vf_integr)
       !! check

    endif

  end subroutine bt2bc_coupling

end module mg_coupling