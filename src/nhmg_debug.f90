module nhmg_debug

  use mg_mpi
  use mg_netcdf_out

  implicit none

contains

  !--------------------------------------------------------------
  subroutine nhmg_write_pred_in(nx,ny,nz,Hzba,Hza,Hzha,uba,vba,wba,ua,va,wa,rufrca,rvfrca)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: Hzba
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: Hza
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: Hzha
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: uba
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: vba
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz+1), target, intent(in) :: wba
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: ua
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: va
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz+1), target, intent(in) :: wa
    real(kind=rp), dimension(-1:nx+2,-1:ny+2),        target, intent(in) :: rufrca
    real(kind=rp), dimension(-1:nx+2,-1:ny+2),        target, intent(in) :: rvfrca

    real(kind=rp), dimension(:,:,:), pointer :: Hzb,Hz,Hzh,ub,vb,wb,u,v,w
    real(kind=rp), dimension(:,:),   pointer :: rufrc,rvfrc

    integer(kind=ip), save :: iter_write_pred_in=0
    iter_write_pred_in = iter_write_pred_in + 1

    Hzb => Hzba
    Hz => Hza
    Hzh => Hzha
    ub => uba
    vb => vba
    wb => wba
    u => ua
    v => va
    w => wa
    rufrc => rufrca
    rvfrc => rvfrca

    call write_netcdf(Hzb,vname='Hzb',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)
    call write_netcdf(Hz,vname='Hz',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)
    call write_netcdf(Hzh,vname='Hzh',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)
    call write_netcdf(ub,vname='ub',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)
    call write_netcdf(vb,vname='vb',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)
    call write_netcdf(wb,vname='wb',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)
    call write_netcdf(u,vname='u',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)
    call write_netcdf(v,vname='v',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)
    call write_netcdf(w,vname='w',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)
    call write_netcdf(rufrc,vname='rufrc',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)
    call write_netcdf(rvfrc,vname='rvfrc',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)

  end subroutine nhmg_write_pred_in

  !--------------------------------------------------------------
  subroutine nhmg_write_pred_inter(nx,ny,nz,ua,va,wa,rufrca,rvfrca)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: ua
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: va
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz+1), target, intent(in) :: wa
    real(kind=rp), dimension(-1:nx+2,-1:ny+2),        target, intent(in) :: rufrca
    real(kind=rp), dimension(-1:nx+2,-1:ny+2),        target, intent(in) :: rvfrca

    real(kind=rp), dimension(:,:,:), pointer :: u,v,w
    real(kind=rp), dimension(:,:),   pointer :: rufrc,rvfrc

    integer(kind=ip), save :: iter_write_pred_inter=0
    iter_write_pred_inter = iter_write_pred_inter + 1

    u => ua
    v => va
    w => wa
    rufrc => rufrca
    rvfrc => rvfrca

    call write_netcdf(u,vname='u',netcdf_file_name='wrt_pred_inter.nc',rank=myrank,iter=iter_write_pred_inter)
    call write_netcdf(v,vname='v',netcdf_file_name='wrt_pred_inter.nc',rank=myrank,iter=iter_write_pred_inter)
    call write_netcdf(w,vname='w',netcdf_file_name='wrt_pred_inter.nc',rank=myrank,iter=iter_write_pred_inter)
    call write_netcdf(rufrc,vname='rufrc',netcdf_file_name='wrt_pred_inter.nc',rank=myrank,iter=iter_write_pred_inter)
    call write_netcdf(rvfrc,vname='rvfrc',netcdf_file_name='wrt_pred_inter.nc',rank=myrank,iter=iter_write_pred_inter)

  end subroutine nhmg_write_pred_inter

  !--------------------------------------------------------------
  subroutine nhmg_write_pred_out(nx,ny,nz,ua,va,wa)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: ua
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: va
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz+1), target, intent(in) :: wa

    real(kind=rp), dimension(:,:,:), pointer :: u,v,w

    integer(kind=ip), save :: iter_write_pred_out=0
    iter_write_pred_out = iter_write_pred_out + 1

    u => ua
    v => va
    w => wa

    call write_netcdf(u,vname='u',netcdf_file_name='wrt_pred_out.nc',rank=myrank,iter=iter_write_pred_out)
    call write_netcdf(v,vname='v',netcdf_file_name='wrt_pred_out.nc',rank=myrank,iter=iter_write_pred_out)
    call write_netcdf(w,vname='w',netcdf_file_name='wrt_pred_out.nc',rank=myrank,iter=iter_write_pred_out)

  end subroutine nhmg_write_pred_out

  !--------------------------------------------------------------
  subroutine nhmg_write_corr_in(nx,ny,nz,ua,va,wa)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: ua
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: va
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz+1), target, intent(in) :: wa

    real(kind=rp), dimension(:,:,:), pointer :: u,v,w

    integer(kind=ip), save :: iter_write_corr_in=0
    iter_write_corr_in = iter_write_corr_in + 1

    u => ua
    v => va
    w => wa

    call write_netcdf(u,vname='u',netcdf_file_name='wrt_corr_in.nc',rank=myrank,iter=iter_write_corr_in)
    call write_netcdf(v,vname='v',netcdf_file_name='wrt_corr_in.nc',rank=myrank,iter=iter_write_corr_in)
    call write_netcdf(w,vname='w',netcdf_file_name='wrt_corr_in.nc',rank=myrank,iter=iter_write_corr_in)

  end subroutine nhmg_write_corr_in

  !--------------------------------------------------------------
  subroutine nhmg_write_corr_inter(nx,ny,nz,ua,va,wa)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: ua
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: va
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz+1), target, intent(in) :: wa

    real(kind=rp), dimension(:,:,:), pointer :: u,v,w

    integer(kind=ip), save :: iter_write_corr_inter=0
    iter_write_corr_inter = iter_write_corr_inter + 1

    u => ua
    v => va
    w => wa

    call write_netcdf(u,vname='u',netcdf_file_name='wrt_corr_inter.nc',rank=myrank,iter=iter_write_corr_inter)
    call write_netcdf(v,vname='v',netcdf_file_name='wrt_corr_inter.nc',rank=myrank,iter=iter_write_corr_inter)
    call write_netcdf(w,vname='w',netcdf_file_name='wrt_corr_inter.nc',rank=myrank,iter=iter_write_corr_inter)

  end subroutine nhmg_write_corr_inter

  !--------------------------------------------------------------
  subroutine nhmg_write_corr_out(nx,ny,nz,ua,va,wa)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: ua
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: va
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz+1), target, intent(in) :: wa

    real(kind=rp), dimension(:,:,:), pointer :: u,v,w

    integer(kind=ip), save :: iter_write_corr_out=0
    iter_write_corr_out = iter_write_corr_out + 1

    u => ua
    v => va
    w => wa

    call write_netcdf(u,vname='u',netcdf_file_name='wrt_corr_out.nc',rank=myrank,iter=iter_write_corr_out)
    call write_netcdf(v,vname='v',netcdf_file_name='wrt_corr_out.nc',rank=myrank,iter=iter_write_corr_out)
    call write_netcdf(w,vname='w',netcdf_file_name='wrt_corr_out.nc',rank=myrank,iter=iter_write_corr_out)

  end subroutine nhmg_write_corr_out

  !--------------------------------------------------------------
  subroutine h_write_pred_in(nx,ny,nz,Hzba,Hza,Hzha,uba,vba,ua,va)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: Hzba
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: Hza
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: Hzha
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: uba
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: vba
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: ua
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: va

    real(kind=rp), dimension(:,:,:), pointer :: Hzb,Hz,Hzh,ub,vb,u,v

    integer(kind=ip), save :: iter_write_pred_in=0
    iter_write_pred_in = iter_write_pred_in + 1

    Hzb => Hzba
    Hz => Hza
    Hzh => Hzha
    ub => uba
    vb => vba
    u => ua
    v => va

    call write_netcdf(Hzb,vname='Hzb',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)
    call write_netcdf(Hz,vname='Hz',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)
    call write_netcdf(Hzh,vname='Hzh',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)
    call write_netcdf(ub,vname='ub',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)
    call write_netcdf(vb,vname='vb',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)
    call write_netcdf(u,vname='u',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)
    call write_netcdf(v,vname='v',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)

  end subroutine h_write_pred_in

  !--------------------------------------------------------------
  subroutine h_write_pred_inter(nx,ny,nz,ua,va,rua,rva)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: ua
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: va
    real(kind=rp), dimension(-1:nx+2,-1:ny+2),        target, intent(in) :: rua
    real(kind=rp), dimension(-1:nx+2,-1:ny+2),        target, intent(in) :: rva

    real(kind=rp), dimension(:,:,:), pointer :: u,v
    real(kind=rp), dimension(:,:),   pointer :: ru,rv

    integer(kind=ip), save :: iter_write_pred_inter=0
    iter_write_pred_inter = iter_write_pred_inter + 1

    u => ua
    v => va
    ru => rua
    rv => rva

    call write_netcdf(u,vname='u',netcdf_file_name='wrt_pred_inter.nc',rank=myrank,iter=iter_write_pred_inter)
    call write_netcdf(v,vname='v',netcdf_file_name='wrt_pred_inter.nc',rank=myrank,iter=iter_write_pred_inter)
    call write_netcdf(ru,vname='ru',netcdf_file_name='wrt_pred_inter.nc',rank=myrank,iter=iter_write_pred_inter)
    call write_netcdf(rv,vname='rv',netcdf_file_name='wrt_pred_inter.nc',rank=myrank,iter=iter_write_pred_inter)

  end subroutine h_write_pred_inter

  !--------------------------------------------------------------
  subroutine h_write_pred_out(nx,ny,nz,ua,va)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: ua
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: va

    real(kind=rp), dimension(:,:,:), pointer :: u,v

    integer(kind=ip), save :: iter_write_pred_out=0
    iter_write_pred_out = iter_write_pred_out + 1

    u => ua
    v => va

    call write_netcdf(u,vname='u',netcdf_file_name='wrt_pred_out.nc',rank=myrank,iter=iter_write_pred_out)
    call write_netcdf(v,vname='v',netcdf_file_name='wrt_pred_out.nc',rank=myrank,iter=iter_write_pred_out)

  end subroutine h_write_pred_out

  !--------------------------------------------------------------
  subroutine h_write_corr_in(nx,ny,nz,ua,va)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: ua
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: va

    real(kind=rp), dimension(:,:,:), pointer :: u,v

    integer(kind=ip), save :: iter_write_corr_in=0
    iter_write_corr_in = iter_write_corr_in + 1

    u => ua
    v => va

    call write_netcdf(u,vname='u',netcdf_file_name='wrt_corr_in.nc',rank=myrank,iter=iter_write_corr_in)
    call write_netcdf(v,vname='v',netcdf_file_name='wrt_corr_in.nc',rank=myrank,iter=iter_write_corr_in)

  end subroutine h_write_corr_in

  !--------------------------------------------------------------
  subroutine h_write_corr_inter(nx,ny,nz,ua,va)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: ua
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: va

    real(kind=rp), dimension(:,:,:), pointer :: u,v

    integer(kind=ip), save :: iter_write_corr_inter=0
    iter_write_corr_inter = iter_write_corr_inter + 1

    u => ua
    v => va

    call write_netcdf(u,vname='u',netcdf_file_name='wrt_corr_inter.nc',rank=myrank,iter=iter_write_corr_inter)
    call write_netcdf(v,vname='v',netcdf_file_name='wrt_corr_inter.nc',rank=myrank,iter=iter_write_corr_inter)

  end subroutine h_write_corr_inter

  !--------------------------------------------------------------
  subroutine h_write_corr_out(nx,ny,nz,ua,va)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: ua
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: va

    real(kind=rp), dimension(:,:,:), pointer :: u,v

    integer(kind=ip), save :: iter_write_corr_out=0
    iter_write_corr_out = iter_write_corr_out + 1

    u => ua
    v => va

    call write_netcdf(u,vname='u',netcdf_file_name='wrt_corr_out.nc',rank=myrank,iter=iter_write_corr_out)
    call write_netcdf(v,vname='v',netcdf_file_name='wrt_corr_out.nc',rank=myrank,iter=iter_write_corr_out)

  end subroutine h_write_corr_out

end module nhmg_debug

