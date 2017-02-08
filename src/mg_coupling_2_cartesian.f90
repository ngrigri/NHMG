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
  subroutine bc2bt_coupling(ru,rv,rw,dt,rufrc,rvfrc)

    real(kind=rp), dimension(:,:,:), pointer, intent(in)  :: ru,rv,rw
    real(kind=rp),                            intent(in)  :: dt
    real(kind=rp), dimension(:,:)  , pointer, intent(out) :: rufrc,rvfrc

    integer(kind=ip):: k, j, i
    integer(kind=ip):: nx, ny, nz

    real(kind=rp), dimension(:,:)  , pointer :: Arz
    real(kind=rp), dimension(:,:,:), pointer :: Arx,Ary
    real(kind=rp), dimension(:,:,:), pointer :: p

    real(kind=rp), dimension(:,:,:), pointer :: rutmp,rvtmp

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz

    Arx   => grid(1)%Arx
    Ary   => grid(1)%Ary
    Arz   => grid(1)%Arz
    p     => grid(1)%p

    !------------------------------------------
    ! add the volume integrated grad of nh pressure

    allocate(rutmp(1:nz,0:ny+1,1:nx+1))
    allocate(rvtmp(1:nz,1:ny+1,0:nx+1))

    do i = 1,nx+1
       do j = 0,ny+1
          do k = 1,nz
             rutmp(k,j,i) = ru(k,j,i) !- Arx(k,j,i) * (p(k,j,i  )-p(k,j,i-1)) / dt
          enddo
       enddo
    enddo

    do i = 0,nx+1
       do j = 1,ny+1
          do k = 1,nz
             rvtmp(k,j,i) = rv(k,j,i) !- Ary(k,j,i) * (p(k,j  ,i)-p(k,j-1,i)) / dt
          enddo
       enddo
    enddo

    !------------------------------------------
    ! fluxes and integrate

    do i = 1,nx+1
       do j = 1,ny

          k = 1
          rufrc(j,i) = rutmp(k,j,i)

          do k = 2,nz-1
             rufrc(j,i) = rufrc(j,i) + rutmp(k,j,i)
          enddo

          k = nz
          rufrc(j,i) = rufrc(j,i) + rutmp(k,j,i)

       enddo
    enddo

    do i = 1,nx
       do j = 1,ny+1 

          k = 1
          rvfrc(j,i) = rvtmp(k,j,i)

          do k = 2,nz-1
             rvfrc(j,i) = rvfrc(j,i) + rvtmp(k,j,i)
          enddo
          
          k = nz
          rvfrc(j,i) = rvfrc(j,i) + rvtmp(k,j,i)

       enddo
    enddo

    deallocate(rutmp)
    deallocate(rvtmp)

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
          su_integr(j,i) = Arx(k,j,i)
          uf_integr(j,i) = Arx(k,j,i) * u(k,j,i)

          do k = 2,nz-1
             su_integr(j,i) = su_integr(j,i) + Arx(k,j,i)
             uf_integr(j,i) = uf_integr(j,i) + Arx(k,j,i) * u(k,j,i)
          enddo

          k = nz
          su_integr(j,i) = su_integr(j,i) + Arx(k,j,i)
          uf_integr(j,i) = uf_integr(j,i) + Arx(k,j,i) * u(k,j,i)

       enddo
    enddo

    do i = 0,nx+1
       do j = 1,ny+1

          k = 1
          sv_integr(j,i) = Ary(k,j,i)
          vf_integr(j,i) = Ary(k,j,i) * v(k,j,i)

          do k = 2,nz-1
             sv_integr(j,i) = sv_integr(j,i) + Ary(k,j,i)
             vf_integr(j,i) = vf_integr(j,i) + Ary(k,j,i) * v(k,j,i)
          enddo

          k = nz
          sv_integr(j,i) = sv_integr(j,i) + Ary(k,j,i)
          vf_integr(j,i) = vf_integr(j,i) + Ary(k,j,i) * v(k,j,i)

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
!ND 07/02 : which definition for cartesian case?

    allocate(wc(1:nz+1,0:ny+1,0:nx+1))

    do i = 1,nx
       do j = 1,ny

          do k = 1,nz+1 
             wc(k,j,i) = zero
          enddo

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

          do k = 1,nz
             u(k,j,i) = u(k,j,i) - (uf_integr(j,i)/su_integr(j,i))
          enddo

       enddo
    enddo

    do i = 0,nx+1
       do j = 1,ny+1

          do k = 1,nz
             v(k,j,i) = v(k,j,i) - (vf_integr(j,i)/sv_integr(j,i))
          enddo

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

             do k = 1,nz
                uf(k,j,i) = Arx(k,j,i) * u(k,j,i)
             enddo

          enddo
       enddo

       do i = 0,nx+1
          do j = 1,ny+1

             do k = 1,nz
                vf(k,j,i) = Ary(k,j,i) * v(k,j,i) 
             enddo

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
