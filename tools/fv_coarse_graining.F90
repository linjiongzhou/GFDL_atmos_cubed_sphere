!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the FV3 dynamical core.
!*
!* The FV3 dynamical core is free software: you can redistribute it
!* and/or modify it under the terms of the
!* GNU Lesser General Public License as published by the
!* Free Software Foundation, either version 3 of the License, or
!* (at your option) any later version.
!*
!* The FV3 dynamical core is distributed in the hope that it will be
!* useful, but WITHOUT ANYWARRANTY; without even the implied warranty
!* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!* See the GNU General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with the FV3 dynamical core.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
module fv_coarse_graining_mod

  use constants_mod,   only: grav, rdgas, pi=>pi_8
  use diag_manager_mod, only: diag_axis_init, register_diag_field, register_static_field, send_data
  use fms_io_mod,      only: register_restart_field, save_restart
  use fv_arrays_mod,   only: fv_atmos_type
  use fv_mp_mod,       only: domain_decomp
  use mpp_mod,         only: mpp_error, mpp_npes, FATAL, input_nml_file
  use mpp_domains_mod, only: mpp_get_compute_domain, mpp_define_mosaic, domain2d, mpp_define_io_domain
  use mpp_domains_mod, only: EAST, NORTH
  use time_manager_mod, only: time_type
  use tracer_manager_mod, only: get_tracer_index, get_tracer_names, set_tracer_profile
  use field_manager_mod,  only: MODEL_ATMOS

  implicit none
  private

  public :: fv_coarse_graining_init, fv_io_write_restart_coarse
  public :: fv_coarse_grained_diagnostics_init, fv_coarse_grained_diagnostics
  public :: block_average, block_edge_average_x, block_edge_average_y

  interface block_average
     module procedure block_average_2d
     module procedure block_average_3d_variable_2d_weights
     module procedure block_average_3d_variable_3d_weights
  end interface block_average

  interface block_edge_average_x
     module procedure block_edge_average_x_2d
     module procedure block_edge_average_x_3d
  end interface block_edge_average_x

  interface block_edge_average_y
     module procedure block_edge_average_y_2d
     module procedure block_edge_average_y_3d
  end interface block_edge_average_y

  type(domain2d) :: coarse_domain
  real, parameter:: missing_value = -1.e10
  real, parameter :: rad2deg = 180. / pi
  integer :: is, ie, js, je, is_coarse, ie_coarse, js_coarse, je_coarse, npz, nt_prog, nt_diag, ntracers
  integer :: target_resolution

  ! Diagnostic IDs
  integer :: id_coarse_lon, id_coarse_lat, id_coarse_lont, id_coarse_latt
  integer :: id_pfnh_coarse
  integer :: id_qv_dt_phys_coarse, id_ql_dt_phys_coarse, id_qi_dt_phys_coarse
  integer :: id_qr_dt_phys_coarse, id_qg_dt_phys_coarse, id_qs_dt_phys_coarse
  integer :: id_t_dt_phys_coarse, id_u_dt_phys_coarse, id_v_dt_phys_coarse
  integer :: id_t_dt_nudge_coarse, id_ps_dt_nudge_coarse, id_delp_dt_nudge_coarse
  integer :: id_u_dt_nudge_coarse, id_v_dt_nudge_coarse
  integer :: id_qv_dt_gfdlmp_coarse, id_ql_dt_gfdlmp_coarse, id_qr_dt_gfdlmp_coarse
  integer :: id_qg_dt_gfdlmp_coarse, id_qi_dt_gfdlmp_coarse, id_qs_dt_gfdlmp_coarse
  integer :: id_t_dt_gfdlmp_coarse, id_u_dt_gfdlmp_coarse, id_v_dt_gfdlmp_coarse
  integer :: id_liq_wat_dt_phys_coarse, id_ice_wat_dt_phys_coarse
  integer :: id_liq_wat_dt_gfdlmp_coarse, id_ice_wat_dt_gfdlmp_coarse

  ! Namelist parameters (with default values)
  integer :: coarsening_factor = 8
  logical :: write_coarse_grained_restart_files = .false.
  logical :: write_coarse_grained_diagnostics = .false.

  namelist /fv_coarse_graining_nml/ coarsening_factor,&
       write_coarse_grained_restart_files, write_coarse_grained_diagnostics

contains

  subroutine fv_coarse_graining_init(Atm, n)
    type(fv_atmos_type), intent(inout) :: Atm
    integer, intent(in) :: n

    character(len=256) :: errmsg
    logical :: exists
    integer :: nlunit, ios

    integer :: native_resolution, np_target

#ifdef INTERNAL_FILE_NML
    read(input_nml_file, nml=fv_coarse_graining_nml,iostat=ios)
#else
    inquire (file=trim(Atm%nml_filename), exist=exists)
    if (.not. exists) then
      write(errmsg,*) 'fv_coarse_graining_nml: namelist file ',trim(Atm%nml_filename),' does not exist'
      call mpp_error(FATAL, errmsg)
    else
      open (unit=nlunit, file=Atm%nml_filename, READONLY, status='OLD', iostat=ios)
    endif
    rewind(nlunit)
    read (nlunit, nml=fv_diag_column_nml, iostat=ios)
    close (nlunit)
#endif

    if (write_coarse_grained_restart_files .or. write_coarse_grained_diagnostics) then

       npz = Atm%npz
       nt_prog = Atm%flagstruct%nt_prog
       nt_diag = Atm%flagstruct%nt_phys
       ntracers = nt_prog + nt_diag

       native_resolution = Atm%flagstruct%npx - 1
       target_resolution = native_resolution / coarsening_factor
       np_target = target_resolution + 1

       ! Check that the coarsening factor evenly divides the native domain
       if (mod(native_resolution, coarsening_factor) .ne. 0) then
          write(errmsg, *) 'fv_coarse_graining_init: coarsening factor ',&
               coarsening_factor, ' does not evenly divide native resolution',&
               native_resolution, '.'
          call mpp_error(FATAL, errmsg)
       endif

!!$       ! Check that compute domains will be square (for simplicity)
!!$       if (Atm%layout(1) .ne. Atm%layout(2)) then
!!$          write(errmsg, *) 'fv_coarse_graining_init: domain decomposition layout ',&
!!$               Atm%layout, ' does not produce square subdomains.'
!!$          call mpp_error(FATAL, errmsg)
!!$       endif

       ! Check that the coarsening factor is appropriate for the domain
       ! decomposition
       if (mod(native_resolution, Atm%layout(1)) > 0 .or. mod(target_resolution, Atm%layout(1)) > 0) then
          write(errmsg, *) 'fv_coarse_graining_init: domain decomposition layout ',&
               Atm%layout(1), ' does not evenly divide both the native and coarse grid.'
          call mpp_error(FATAL, errmsg)
       endif

       call define_cubic_mosaic(coarse_domain, target_resolution,&
            target_resolution, Atm%layout)
       call mpp_define_io_domain(coarse_domain, Atm%io_layout)
       call mpp_get_compute_domain(coarse_domain, is_coarse, ie_coarse, &
            js_coarse, je_coarse)
       call mpp_get_compute_domain(Atm%domain, is, ie, js, je)

    endif

    if (write_coarse_grained_restart_files) then
       call allocate_coarse_restart_type(Atm)
       call register_coarse_grained_restart_files(Atm, n)
    endif

  end subroutine fv_coarse_graining_init

  subroutine allocate_coarse_restart_type(Atm)
    type(fv_atmos_type), intent(INOUT) :: Atm

    allocate (Atm%coarse_restart%u(is_coarse:ie_coarse,js_coarse:je_coarse+1,npz))
    allocate (Atm%coarse_restart%v(is_coarse:ie_coarse+1,js_coarse:je_coarse,npz))
    allocate (Atm%coarse_restart%ua(is_coarse:ie_coarse,js_coarse:je_coarse,npz))
    allocate (Atm%coarse_restart%va(is_coarse:ie_coarse,js_coarse:je_coarse,npz))
    allocate (Atm%coarse_restart%u_srf(is_coarse:ie_coarse,js_coarse:je_coarse))
    allocate (Atm%coarse_restart%v_srf(is_coarse:ie_coarse,js_coarse:je_coarse))
    allocate (Atm%coarse_restart%w(is_coarse:ie_coarse,js_coarse:je_coarse,npz))
    allocate (Atm%coarse_restart%delp(is_coarse:ie_coarse,js_coarse:je_coarse,npz))
    allocate (Atm%coarse_restart%delz(is_coarse:ie_coarse,js_coarse:je_coarse,npz))
    allocate (Atm%coarse_restart%pt(is_coarse:ie_coarse,js_coarse:je_coarse,npz))
    allocate (Atm%coarse_restart%q(is_coarse:ie_coarse,js_coarse:je_coarse,npz,nt_prog))
    allocate (Atm%coarse_restart%qdiag(is_coarse:ie_coarse,js_coarse:je_coarse,npz,nt_prog+1:nt_diag + nt_prog))
    allocate (Atm%coarse_restart%sgh(is_coarse:ie_coarse,js_coarse:je_coarse))
    allocate (Atm%coarse_restart%oro(is_coarse:ie_coarse,js_coarse:je_coarse))
    allocate (Atm%coarse_restart%phis(is_coarse:ie_coarse,js_coarse:je_coarse))
    allocate (Atm%coarse_restart%ze0(is_coarse:ie_coarse,js_coarse:je_coarse,npz))

    Atm%coarse_restart%allocated = .true.

  end subroutine allocate_coarse_restart_type

  subroutine register_coarse_grained_restart_files(Atm, n)
    type(fv_atmos_type), intent(inout) :: Atm
    integer, intent(in) :: n

    character(len=64) :: fname, tracer_name
    integer           :: id_restart
    integer           :: nt, ntracers

    fname = 'fv_core_coarse.res.nc'

    id_restart = register_restart_field(Atm%coarse_restart%fv_tile_restart_coarse, &
         fname, 'u', Atm%coarse_restart%u, &
         domain=coarse_domain, &
         position=NORTH, tile_count=n)
    id_restart = register_restart_field(Atm%coarse_restart%fv_tile_restart_coarse, &
         fname, 'v', Atm%coarse_restart%v, &
         domain=coarse_domain, &
         position=EAST,tile_count=n)

    if (.not. Atm%flagstruct%hydrostatic) then
       id_restart = register_restart_field(Atm%coarse_restart%fv_tile_restart_coarse, &
            fname, 'W', Atm%coarse_restart%w, &
            domain=coarse_domain, &
            mandatory=.false., tile_count=n)
       id_restart = register_restart_field(Atm%coarse_restart%fv_tile_restart_coarse, &
            fname, 'DZ', Atm%coarse_restart%delz, &
            domain=coarse_domain, &
            mandatory=.false., tile_count=n)

       if ( Atm%flagstruct%hybrid_z ) then
          id_restart = register_restart_field(Atm%coarse_restart%fv_tile_restart_coarse, &
               fname, 'ZE0', Atm%coarse_restart%ze0, &
               domain=coarse_domain, mandatory=.false., tile_count=n)
       endif
    endif

    id_restart = register_restart_field(Atm%coarse_restart%fv_tile_restart_coarse, &
         fname, 'T', Atm%coarse_restart%pt, &
         domain=coarse_domain, tile_count=n)
    id_restart = register_restart_field(Atm%coarse_restart%fv_tile_restart_coarse, &
         fname, 'delp', Atm%coarse_restart%delp, &
         domain=coarse_domain, tile_count=n)
    id_restart = register_restart_field(Atm%coarse_restart%fv_tile_restart_coarse, &
         fname, 'phis', Atm%coarse_restart%phis, &
         domain=coarse_domain, tile_count=n)

    if (Atm%flagstruct%agrid_vel_rst) then
       id_restart = register_restart_field(Atm%coarse_restart%fv_tile_restart_coarse, &
            fname, 'ua', Atm%coarse_restart%ua, &
            domain=coarse_domain, tile_count=n, mandatory=.false.)
       id_restart = register_restart_field(Atm%coarse_restart%fv_tile_restart_coarse, &
            fname, 'va', Atm%coarse_restart%va, &
            domain=coarse_domain, tile_count=n, mandatory=.false.)
    endif

    fname = 'fv_srf_wnd_coarse.res.nc'
    id_restart = register_restart_field(Atm%coarse_restart%rsf_restart_coarse, &
         fname, 'u_srf', Atm%coarse_restart%u_srf, &
         domain=coarse_domain, tile_count=n)
    id_restart = register_restart_field(Atm%coarse_restart%rsf_restart_coarse, &
         fname, 'v_srf', Atm%coarse_restart%v_srf, &
         domain=coarse_domain, tile_count=n)

    if ( Atm%flagstruct%fv_land ) then
       fname = 'mg_drag_coarse.res.nc'
       id_restart = register_restart_field(Atm%coarse_restart%mg_restart_coarse, &
            fname, 'ghprime', Atm%coarse_restart%sgh, &
            domain=coarse_domain, tile_count=n)

       fname = 'fv_land_coarse.res.nc'
       id_restart = register_restart_field(Atm%coarse_restart%lnd_restart_coarse, &
            fname, 'oro', Atm%coarse_restart%oro, &
            domain=coarse_domain, tile_count=n)
    endif

    fname = 'fv_tracer_coarse.res.nc'
    do nt = 1, nt_prog
       call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
       call set_tracer_profile(MODEL_ATMOS, nt, Atm%coarse_restart%q(:,:,:,nt))
       id_restart = register_restart_field(Atm%coarse_restart%tra_restart_coarse, &
            fname,tracer_name, Atm%coarse_restart%q(:,:,:,nt), &
            domain=coarse_domain, mandatory=.false., tile_count=n)
    enddo

    do nt = nt_prog+1, ntracers
       call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
       call set_tracer_profile (MODEL_ATMOS, nt, Atm%coarse_restart%qdiag(:,:,:,nt))
       id_restart = register_restart_field(Atm%coarse_restart%tra_restart_coarse, &
            fname, tracer_name, Atm%coarse_restart%qdiag(:,:,:,nt), &
            domain=coarse_domain, mandatory=.false., tile_count=n)
    enddo

  end subroutine register_coarse_grained_restart_files

  subroutine fv_io_write_restart_coarse(Atm, timestamp)
    type(fv_atmos_type),        intent(inout) :: Atm
    character(len=*), optional, intent(in) :: timestamp

    if (write_coarse_grained_restart_files) then
       if (present(timestamp)) then
          call write_restart_coarse(Atm, timestamp)
       else
          call write_restart_coarse(Atm)
       endif
    endif

  end subroutine fv_io_write_restart_coarse

  subroutine write_restart_coarse(Atm, timestamp)
    type(fv_atmos_type),        intent(inout) :: Atm
    character(len=*), optional, intent(in) :: timestamp

    call coarse_grain_restart_data(Atm)

    ! TODO: do we need to support this?
    ! if ( (use_ncep_sst .or. Atm%flagstruct%nudge) .and. .not. Atm%gridstruct%nested ) then
    !    call save_restart(Atm%SST_restart, timestamp)
    ! endif

    call save_restart(Atm%coarse_restart%fv_tile_restart_coarse, timestamp)
    call save_restart(Atm%coarse_restart%rsf_restart_coarse, timestamp)

    if ( Atm%flagstruct%fv_land ) then
       call save_restart(Atm%coarse_restart%mg_restart_coarse, timestamp)
       call save_restart(Atm%coarse_restart%lnd_restart_coarse, timestamp)
    endif

    call save_restart(Atm%coarse_restart%tra_restart_coarse, timestamp)

  end subroutine write_restart_coarse

  subroutine coarse_grain_restart_data(Atm)
    type(fv_atmos_type), intent(inout) :: Atm
    character(len=64) :: tracer_name
    integer :: k, nt
    real, allocatable :: mass(:,:,:)

    allocate(mass(is:ie,js:je,1:npz))
    do k = 1, npz
       mass(is:ie,js:je,k) = Atm%gridstruct%area(is:ie,js:je) * Atm%delp(is:ie,js:je,k)
    enddo

    call block_edge_average_x(Atm%gridstruct%dx(is:ie,js:je+1), &
          Atm%u(is:ie,js:je+1,1:npz), Atm%coarse_restart%u)
    call block_edge_average_y(Atm%gridstruct%dy(is:ie+1,js:je), &
          Atm%v(is:ie+1,js:je,1:npz), Atm%coarse_restart%v)
    call block_average(mass(is:ie,js:je,1:npz), &
         Atm%pt(is:ie,js:je,1:npz), Atm%coarse_restart%pt)

    if (.not. Atm%flagstruct%hydrostatic) then
       call block_average(mass(is:ie,js:je,1:npz), &
            Atm%w(is:ie,js:je,1:npz), Atm%coarse_restart%w)
       call block_average(Atm%gridstruct%area(is:ie,js:je), &
            Atm%delz(is:ie,js:je,1:npz), Atm%coarse_restart%delz)
       if ( Atm%flagstruct%hybrid_z ) then
          call block_average(Atm%gridstruct%area(is:ie,js:je), &
            Atm%ze0(is:ie,js:je,1:npz), Atm%coarse_restart%ze0)
       endif
    endif

    call block_average(mass(is:ie,js:je,1:npz), &
         Atm%pt(is:ie,js:je,1:npz), Atm%coarse_restart%pt)
    call block_average(Atm%gridstruct%area(is:ie,js:je), &
         Atm%delp(is:ie,js:je,1:npz), Atm%coarse_restart%delp)
    call block_average(Atm%gridstruct%area(is:ie,js:je), &
         Atm%phis(is:ie,js:je), Atm%coarse_restart%phis)

    if (Atm%flagstruct%agrid_vel_rst) then
       call block_average(Atm%gridstruct%area(is:ie,js:je), &
            Atm%ua(is:ie,js:je,1:npz), Atm%coarse_restart%ua)
       call block_average(Atm%gridstruct%area(is:ie,js:je), &
            Atm%va(is:ie,js:je,1:npz), Atm%coarse_restart%va)
    endif

    call block_average(Atm%gridstruct%area(is:ie,js:je), &
         Atm%u_srf(is:ie,js:je), Atm%coarse_restart%u_srf)
    call block_average(Atm%gridstruct%area(is:ie,js:je), &
         Atm%v_srf(is:ie,js:je), Atm%coarse_restart%v_srf)

    if ( Atm%flagstruct%fv_land ) then
       call block_average(Atm%gridstruct%area(is:ie,js:je), &
            Atm%sgh(is:ie,js:je), Atm%coarse_restart%sgh)
       call block_average(Atm%gridstruct%area(is:ie,js:je), &
            Atm%oro(is:ie,js:je), Atm%coarse_restart%oro)
    endif

    do nt = 1, nt_prog
       call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
       if (tracer_name .eq. 'cld_amt') then
          ! cld_amt is cloud fraction (and so is not mass-weighted)
          call block_average(Atm%gridstruct%area(is:ie,js:je), &
               Atm%q(is:ie,js:je,1:npz,nt), &
               Atm%coarse_restart%q(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz,nt))
       else
          call block_average(mass(is:ie,js:je,1:npz), &
               Atm%q(is:ie,js:je,1:npz,nt), &
               Atm%coarse_restart%q(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz,nt))
       endif
    enddo

    do nt = nt_prog+1, ntracers
       call block_average(mass(is:ie,js:je,1:npz), &
            Atm%qdiag(is:ie,js:je,1:npz,nt), &
            Atm%coarse_restart%qdiag(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz,nt))
    enddo

    deallocate(mass)
  end subroutine coarse_grain_restart_data

  subroutine fv_coarse_grained_diagnostics_init(Atm, Time, id_pfull)
    type(fv_atmos_type), intent(inout) :: Atm(:)
    type(time_type), intent(in) :: Time
    integer, intent(in) :: id_pfull

    if (write_coarse_grained_diagnostics) then
       call register_coarse_grained_diagnostics(Atm, Time, id_pfull)
    endif

  end subroutine fv_coarse_grained_diagnostics_init

  subroutine register_coarse_grained_diagnostics(Atm, Time, id_pfull)
    type(fv_atmos_type), intent(inout) :: Atm(:)
    type(time_type), intent(in) :: Time
    integer, intent(in) :: id_pfull

    integer :: coarse_axes(3), coarse_axes_t(3)
    integer :: i, j, n, offset
    integer :: id_coarse_x, id_coarse_y, id_coarse_xt, id_coarse_yt
    logical :: used
    real, allocatable :: grid_x_coarse(:), grid_y_coarse(:), grid_xt_coarse(:), grid_yt_coarse(:)
    real, allocatable :: agrid_coarse(:,:,:), grid_coarse(:,:,:)

    n = 1

    allocate(grid_x_coarse(target_resolution + 1), grid_y_coarse(target_resolution + 1))
    grid_x_coarse = (/ (i, i=1, target_resolution + 1) /)
    grid_y_coarse = (/ (j, j=1, target_resolution + 1) /)

    allocate(grid_xt_coarse(target_resolution), grid_yt_coarse(target_resolution))
    grid_xt_coarse = (/ (i, i=1, target_resolution) /)
    grid_yt_coarse = (/ (j, j=1, target_resolution) /)

    id_coarse_x = diag_axis_init('grid_x_coarse', grid_x_coarse, &
         'degrees_E', 'x', 'Corner longitude', set_name='coarse_grid', &
         Domain2=coarse_domain, tile_count=n)
    id_coarse_y = diag_axis_init('grid_y_coarse', grid_y_coarse, &
         'degrees_N', 'y', 'Corner latitude', set_name='coarse_grid', &
         Domain2=coarse_domain, tile_count=n)

    id_coarse_xt = diag_axis_init('grid_xt_coarse', grid_xt_coarse, &
         'degrees_E', 'x', 'T-cell longitude', set_name='coarse_grid', &
         Domain2=coarse_domain, tile_count=n)
    id_coarse_yt = diag_axis_init('grid_yt_coarse', grid_yt_coarse, &
         'degrees_N', 'y', 'T-cell latitude', set_name='coarse_grid', &
         Domain2=coarse_domain, tile_count=n)

    deallocate(grid_x_coarse, grid_y_coarse, grid_xt_coarse, grid_yt_coarse)

    coarse_axes(1) = id_coarse_x
    coarse_axes(2) = id_coarse_y
    coarse_axes(3) = id_pfull

    coarse_axes_t(1) = id_coarse_xt
    coarse_axes_t(2) = id_coarse_yt
    coarse_axes_t(3) = id_pfull

    if( .not. Atm(n)%flagstruct%hydrostatic ) then
       id_pfnh_coarse = register_diag_field('dynamics', 'pfnh_coarse', &
            coarse_axes_t(1:3), Time, &
            'non-hydrostatic pressure', &
            'pa', missing_value=missing_value)
    endif

    id_qv_dt_phys_coarse = register_diag_field('dynamics', &
         'qv_dt_phys_coarse', coarse_axes_t(1:3), Time, &
         'water vapor specific humidity tendency from physics', &
         'kg/kg/s', missing_value=missing_value)
    if ((id_qv_dt_phys_coarse > 0) .and. (.not. allocated(Atm(n)%phys_diag%phys_qv_dt))) then
       allocate (Atm(n)%phys_diag%phys_qv_dt(is:ie,js:je,npz))
    endif

    id_ql_dt_phys_coarse = register_diag_field('dynamics', &
         'ql_dt_phys_coarse', coarse_axes_t(1:3), Time, &
         'total liquid water tendency from physics', &
         'kg/kg/s', missing_value=missing_value)
    if ((id_ql_dt_phys_coarse > 0) .and. (.not. allocated(Atm(n)%phys_diag%phys_ql_dt))) then
       allocate (Atm(n)%phys_diag%phys_ql_dt(is:ie,js:je,npz))
    endif

    id_qi_dt_phys_coarse = register_diag_field('dynamics', &
         'qi_dt_phys_coarse', coarse_axes_t(1:3), Time, &
         'total ice water tendency from physics', &
         'kg/kg/s', missing_value=missing_value)
    if ((id_qi_dt_phys_coarse > 0) .and. (.not. allocated(Atm(n)%phys_diag%phys_qi_dt))) then
       allocate (Atm(n)%phys_diag%phys_qi_dt(is:ie,js:je,npz))
    endif

    id_liq_wat_dt_phys_coarse = register_diag_field('dynamics', &
         'liq_wat_dt_phys_coarse', coarse_axes_t(1:3), Time, &
         'liquid water tracer tendency from physics', &
         'kg/kg/s', missing_value=missing_value)
    if ((id_liq_wat_dt_phys_coarse > 0) .and. (.not. allocated(Atm(n)%phys_diag%phys_liq_wat_dt))) then
       allocate (Atm(n)%phys_diag%phys_liq_wat_dt(is:ie,js:je,npz))
    endif

    id_ice_wat_dt_phys_coarse = register_diag_field('dynamics', &
         'ice_wat_dt_phys_coarse', coarse_axes_t(1:3), Time, &
         'ice water tendency from physics', &
         'kg/kg/s', missing_value=missing_value)
    if ((id_ice_wat_dt_phys_coarse > 0) .and. (.not. allocated(Atm(n)%phys_diag%phys_ice_wat_dt))) then
       allocate (Atm(n)%phys_diag%phys_ice_wat_dt(is:ie,js:je,npz))
    endif

    id_qr_dt_phys_coarse = register_diag_field('dynamics', &
         'qr_dt_phys_coarse', coarse_axes_t(1:3), Time, &
         'rain water tendency from physics', &
         'kg/kg/s', missing_value=missing_value)
    if ((id_qr_dt_phys_coarse > 0) .and. (.not. allocated(Atm(n)%phys_diag%phys_qr_dt))) then
       allocate (Atm(n)%phys_diag%phys_qr_dt(is:ie,js:je,npz))
    endif

    id_qg_dt_phys_coarse = register_diag_field('dynamics', &
         'qg_dt_phys_coarse', coarse_axes_t(1:3), Time, &
         'graupel tendency from physics', &
         'kg/kg/s', missing_value=missing_value)
    if ((id_qg_dt_phys_coarse > 0) .and. (.not. allocated(Atm(n)%phys_diag%phys_qg_dt))) then
       allocate (Atm(n)%phys_diag%phys_qg_dt(is:ie,js:je,npz))
    endif

    id_qs_dt_phys_coarse = register_diag_field('dynamics', &
         'qs_dt_phys_coarse', coarse_axes_t(1:3), Time, &
         'snow water tendency from physics', &
         'kg/kg/s', missing_value=missing_value)
    if ((id_qs_dt_phys_coarse > 0) .and. (.not. allocated(Atm(n)%phys_diag%phys_qs_dt))) then
       allocate (Atm(n)%phys_diag%phys_qs_dt(is:ie,js:je,npz))
    endif

    id_t_dt_phys_coarse = register_diag_field('dynamics', &
         't_dt_phys_coarse', coarse_axes_t(1:3), Time, &
         'temperature tendency from physics', &
         'K/s', missing_value=missing_value)
    if ((id_t_dt_phys_coarse > 0) .and. (.not. allocated(Atm(n)%phys_diag%phys_t_dt))) then
       allocate (Atm(n)%phys_diag%phys_t_dt(is:ie,js:je,npz))
    endif

    id_u_dt_phys_coarse = register_diag_field('dynamics', &
         'u_dt_phys_coarse', coarse_axes_t(1:3), Time, &
         'zonal wind tendency from physics', &
         'm/s/s', missing_value=missing_value)
    if ((id_u_dt_phys_coarse > 0) .and. (.not. allocated(Atm(n)%phys_diag%phys_u_dt))) then
       allocate (Atm(n)%phys_diag%phys_u_dt(is:ie,js:je,npz))
    endif

    id_v_dt_phys_coarse = register_diag_field('dynamics', &
         'v_dt_phys_coarse', coarse_axes_t(1:3), Time, &
         'meridional wind tendency from physics', &
         'm/s/s', missing_value=missing_value)
    if ((id_v_dt_phys_coarse > 0) .and. (.not. allocated(Atm(n)%phys_diag%phys_v_dt))) then
       allocate (Atm(n)%phys_diag%phys_v_dt(is:ie,js:je,npz))
    endif

    ! Add coarse-grain tendencies from nudging
    id_t_dt_nudge_coarse = register_diag_field('dynamics', &
         't_dt_nudge_coarse', coarse_axes_t(1:3), Time, &
         'temperature tendency from nudging', &
         'K/s', missing_value=missing_value)
    if ((id_t_dt_nudge_coarse > 0) .and. (.not. allocated(Atm(n)%nudge_diag%nudge_t_dt))) then
       allocate (Atm(n)%nudge_diag%nudge_t_dt(is:ie,js:je,npz))
    endif

    id_ps_dt_nudge_coarse = register_diag_field('dynamics', &
         'ps_dt_nudge_coarse', coarse_axes_t(1:2), Time, &
         'surface pressure tendency from nudging', &
         'Pa/s', missing_value=missing_value)
    if ((id_ps_dt_nudge_coarse > 0) .and. (.not. allocated(Atm(n)%nudge_diag%nudge_ps_dt))) then
       allocate (Atm(n)%nudge_diag%nudge_ps_dt(is:ie,js:je))
    endif

    id_delp_dt_nudge_coarse = register_diag_field('dynamics', &
         'delp_dt_nudge_coarse', coarse_axes_t(1:3), Time, &
         'pressure thickness tendency from nudging', &
         'Pa/s', missing_value=missing_value)
    if ((id_delp_dt_nudge_coarse > 0) .and. (.not. allocated(Atm(n)%nudge_diag%nudge_delp_dt))) then
       allocate (Atm(n)%nudge_diag%nudge_delp_dt(is:ie,js:je,npz))
    endif

    id_u_dt_nudge_coarse = register_diag_field('dynamics', &
         'u_dt_nudge_coarse', coarse_axes_t(1:3), Time, &
         'zonal wind tendency from nudging', &
         'm/s/s', missing_value=missing_value)
    if ((id_u_dt_nudge_coarse > 0) .and. (.not. allocated(Atm(n)%nudge_diag%nudge_u_dt))) then
       allocate (Atm(n)%nudge_diag%nudge_u_dt(is:ie,js:je,npz))
    endif

    id_v_dt_nudge_coarse = register_diag_field('dynamics', &
         'v_dt_nudge_coarse', coarse_axes_t(1:3), Time, &
         'meridional wind tendency from nudging', &
         'm/s/s', missing_value=missing_value)
    if ((id_v_dt_nudge_coarse > 0) .and. (.not. allocated(Atm(n)%nudge_diag%nudge_v_dt))) then
       allocate (Atm(n)%nudge_diag%nudge_v_dt(is:ie,js:je,npz))
    endif

    ! Add coarse-grained GFDL MP tendencies
    id_qv_dt_gfdlmp_coarse = register_diag_field('dynamics', &
         'qv_dt_gfdlmp_coarse', coarse_axes_t(1:3), Time, &
         'water vapor tendency from GFDL MP', &
         'kg/kg/s', missing_value=missing_value)
    if ((id_qv_dt_gfdlmp_coarse > 0) .and. (.not. allocated(Atm(n)%inline_mp%qv_dt))) then
       allocate (Atm(n)%inline_mp%qv_dt(is:ie,js:je,npz))
    endif

    id_ql_dt_gfdlmp_coarse = register_diag_field('dynamics', &
         'ql_dt_gfdlmp_coarse', coarse_axes_t(1:3), Time, &
         'total liquid water tendency from GFDL MP', &
         'kg/kg/s', missing_value=missing_value)
    if ((id_ql_dt_gfdlmp_coarse > 0) .and. (.not. allocated(Atm(n)%inline_mp%ql_dt))) then
       allocate (Atm(n)%inline_mp%ql_dt(is:ie,js:je,npz))
    endif

    id_qi_dt_gfdlmp_coarse = register_diag_field('dynamics', &
         'qi_dt_gfdlmp_coarse', coarse_axes_t(1:3), Time, &
         'total ice water tendency from GFDL MP', &
         'kg/kg/s', missing_value=missing_value)
    if ((id_qi_dt_gfdlmp_coarse > 0) .and. (.not. allocated(Atm(n)%inline_mp%qi_dt))) then
       allocate (Atm(n)%inline_mp%qi_dt(is:ie,js:je,npz))
    endif

    id_liq_wat_dt_gfdlmp_coarse = register_diag_field('dynamics', &
         'liq_wat_dt_gfdlmp_coarse', coarse_axes_t(1:3), Time, &
         'liquid water tracer tendency from GFDL MP', &
         'kg/kg/s', missing_value=missing_value)
    if ((id_liq_wat_dt_gfdlmp_coarse > 0) .and. (.not. allocated(Atm(n)%inline_mp%liq_wat_dt))) then
       allocate (Atm(n)%inline_mp%liq_wat_dt(is:ie,js:je,npz))
    endif

    id_ice_wat_dt_gfdlmp_coarse = register_diag_field('dynamics', &
         'ice_wat_dt_gfdlmp_coarse', coarse_axes_t(1:3), Time, &
         'ice water tracer tendency from GFDL MP', &
         'kg/kg/s', missing_value=missing_value)
    if ((id_ice_wat_dt_gfdlmp_coarse > 0) .and. (.not. allocated(Atm(n)%inline_mp%ice_wat_dt))) then
       allocate (Atm(n)%inline_mp%ice_wat_dt(is:ie,js:je,npz))
    endif

    id_qr_dt_gfdlmp_coarse = register_diag_field('dynamics', &
         'qr_dt_gfdlmp_coarse', coarse_axes_t(1:3), Time, &
         'rain water tendency from GFDL MP', &
         'kg/kg/s', missing_value=missing_value)
    if ((id_qr_dt_gfdlmp_coarse > 0) .and. (.not. allocated(Atm(n)%inline_mp%qr_dt))) then
       allocate (Atm(n)%inline_mp%qr_dt(is:ie,js:je,npz))
    endif

    id_qg_dt_gfdlmp_coarse = register_diag_field('dynamics', &
         'qg_dt_gfdlmp_coarse', coarse_axes_t(1:3), Time, &
         'graupel tendency from GFDL MP', &
         'kg/kg/s', missing_value=missing_value)
    if ((id_qg_dt_gfdlmp_coarse > 0) .and. (.not. allocated(Atm(n)%inline_mp%qg_dt))) then
       allocate (Atm(n)%inline_mp%qg_dt(is:ie,js:je,npz))
    endif

    id_qs_dt_gfdlmp_coarse = register_diag_field('dynamics', &
         'qs_dt_gfdlmp_coarse', coarse_axes_t(1:3), Time, &
         'snow water tendency from GFDL MP', &
         'kg/kg/s', missing_value=missing_value)
    if ((id_qs_dt_gfdlmp_coarse > 0) .and. (.not. allocated(Atm(n)%inline_mp%qs_dt))) then
       allocate (Atm(n)%inline_mp%qs_dt(is:ie,js:je,npz))
    endif

    id_t_dt_gfdlmp_coarse = register_diag_field('dynamics', &
         't_dt_gfdlmp_coarse', coarse_axes_t(1:3), Time, &
         'temperature tendency from GFDL MP', &
         'K/s', missing_value=missing_value)
    if ((id_t_dt_gfdlmp_coarse > 0) .and. (.not. allocated(Atm(n)%inline_mp%T_dt))) then
       allocate (Atm(n)%inline_mp%T_dt(is:ie,js:je,npz))
    endif

    id_u_dt_gfdlmp_coarse = register_diag_field('dynamics', &
         'u_dt_gfdlmp_coarse', coarse_axes_t(1:3), Time, &
         'zonal wind tendency from GFDL MP', &
         'm/s/s', missing_value=missing_value)
    if ((id_u_dt_gfdlmp_coarse > 0) .and. (.not. allocated(Atm(n)%inline_mp%u_dt))) then
       allocate (Atm(n)%inline_mp%u_dt(is:ie,js:je,npz))
    endif

    id_v_dt_gfdlmp_coarse = register_diag_field('dynamics', &
         'v_dt_gfdlmp_coarse', coarse_axes_t(1:3), Time, &
         'zonal wind tendency from GFDL MP', &
         'm/s/s', missing_value=missing_value)
    if ((id_v_dt_gfdlmp_coarse > 0) .and. (.not. allocated(Atm(n)%inline_mp%v_dt))) then
       allocate (Atm(n)%inline_mp%v_dt(is:ie,js:je,npz))
    endif

    ! Construct coarse latitude and longitude coordinates
    allocate(grid_coarse(is_coarse:ie_coarse+1,js_coarse:je_coarse+1, 1:2))
    allocate(agrid_coarse(is_coarse:ie_coarse,js_coarse:je_coarse, 1:2))

    grid_coarse = Atm(n)%gridstruct%grid(is:ie+1:coarsening_factor,js:je+1:coarsening_factor,:)

    ! agrid_coarse is constructed differently depending on whether coarsening
    ! factor is even or odd (in the even case the values are taken from the
    ! native grid corners; in the odd case the values are taken from the native
    ! grid centers).
    offset = coarsening_factor / 2
    if (mod(coarsening_factor, 2) .eq. 0) then
       agrid_coarse = Atm(n)%gridstruct%grid(is+offset:ie+1:coarsening_factor,js+offset:je+1:coarsening_factor,:)
    else
       agrid_coarse = Atm(n)%gridstruct%agrid(is+offset:ie:coarsening_factor,js+offset:je:coarsening_factor,:)
    endif

    id_coarse_lon = register_static_field('dynamics', &
         'grid_lon_coarse', coarse_axes(1:2), &
         'longitude', 'degrees_E')
    id_coarse_lat = register_static_field('dynamics', &
         'grid_lat_coarse', coarse_axes(1:2), &
         'latitude', 'degrees_N')

    id_coarse_lont = register_static_field('dynamics', &
         'grid_lont_coarse', coarse_axes_t(1:2), &
         'longitude', 'degrees_E')
    id_coarse_latt = register_static_field('dynamics', &
         'grid_latt_coarse', coarse_axes_t(1:2), &
         'latitude', 'degrees_N')

    if (id_coarse_lon  > 0) used = send_data(id_coarse_lon, &
         rad2deg * grid_coarse(is_coarse:ie_coarse+1,js_coarse:je_coarse+1,1), Time)
    if (id_coarse_lat  > 0) used = send_data(id_coarse_lat, &
         rad2deg * grid_coarse(is_coarse:ie_coarse+1,js_coarse:je_coarse+1,2), Time)
    if (id_coarse_lont > 0) used = send_data(id_coarse_lont, &
         rad2deg * agrid_coarse(is_coarse:ie_coarse,js_coarse:je_coarse,1), Time)
    if (id_coarse_latt > 0) used = send_data(id_coarse_latt, &
         rad2deg * agrid_coarse(is_coarse:ie_coarse,js_coarse:je_coarse,2), Time)

    deallocate(grid_coarse, agrid_coarse)

  end subroutine register_coarse_grained_diagnostics

  subroutine fv_coarse_grained_diagnostics(Atm, Time, zvir)
    type(fv_atmos_type), intent(in), target :: Atm(:)
    type(time_type), intent(in) :: Time
    real, intent(in) :: zvir

    if (write_coarse_grained_diagnostics) then
       call compute_coarse_grained_diagnostics(Atm, Time, zvir)
    endif
  end subroutine fv_coarse_grained_diagnostics

  subroutine compute_coarse_grained_diagnostics(Atm, Time, zvir)
    type(fv_atmos_type), intent(in), target :: Atm(:)
    type(time_type), intent(in) :: Time
    real, intent(in) :: zvir

    logical :: used
    integer :: var_mass_weighted(19), var_2d(1), var_3d(27)
    integer :: k
    integer :: n = 1
    real, allocatable :: work_2d(:,:), work_3d(:,:,:), mass(:,:,:), nhpres(:,:,:)

    var_2d = (/ id_ps_dt_nudge_coarse /)

    var_3d = (/ &
         id_u_dt_phys_coarse, &
         id_v_dt_phys_coarse, &
         id_t_dt_phys_coarse, &
         id_qv_dt_phys_coarse, &
         id_ql_dt_phys_coarse, &
         id_qi_dt_phys_coarse, &
         id_qr_dt_phys_coarse, &
         id_qg_dt_phys_coarse, &
         id_qs_dt_phys_coarse, &
         id_t_dt_nudge_coarse, &
         id_u_dt_nudge_coarse, &
         id_v_dt_nudge_coarse, &
         id_delp_dt_nudge_coarse, &
         id_qv_dt_gfdlmp_coarse, &
         id_ql_dt_gfdlmp_coarse, &
         id_qr_dt_gfdlmp_coarse, &
         id_qi_dt_gfdlmp_coarse, &
         id_qg_dt_gfdlmp_coarse, &
         id_qs_dt_gfdlmp_coarse, &
         id_t_dt_gfdlmp_coarse, &
         id_u_dt_gfdlmp_coarse, &
         id_v_dt_gfdlmp_coarse, &
         id_liq_wat_dt_gfdlmp_coarse, &
         id_ice_wat_dt_gfdlmp_coarse, &
         id_liq_wat_dt_phys_coarse, &
         id_ice_wat_dt_phys_coarse, &
         id_pfnh_coarse &
         /)

    var_mass_weighted = (/ &
         id_t_dt_phys_coarse, &
         id_qv_dt_phys_coarse, &
         id_ql_dt_phys_coarse, &
         id_qi_dt_phys_coarse, &
         id_qr_dt_phys_coarse, &
         id_qg_dt_phys_coarse, &
         id_qs_dt_phys_coarse, &
         id_t_dt_nudge_coarse, &
         id_qv_dt_gfdlmp_coarse, &
         id_ql_dt_gfdlmp_coarse, &
         id_qr_dt_gfdlmp_coarse, &
         id_qi_dt_gfdlmp_coarse, &
         id_qg_dt_gfdlmp_coarse, &
         id_qs_dt_gfdlmp_coarse, &
         id_t_dt_gfdlmp_coarse, &
         id_liq_wat_dt_gfdlmp_coarse, &
         id_ice_wat_dt_gfdlmp_coarse, &
         id_liq_wat_dt_phys_coarse, &
         id_ice_wat_dt_phys_coarse &
         /)

    if (any(var_2d > 0)) then
       allocate(work_2d(is_coarse:ie_coarse,js_coarse:je_coarse))
    endif

    if (any(var_3d > 0)) then
       allocate(work_3d(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz))
    endif

    if (any(var_mass_weighted > 0)) then
       allocate(mass(is:ie,js:je,1:npz))
       do k = 1, npz
          mass(is:ie,js:je,k) = Atm(n)%gridstruct%area(is:ie,js:je) * Atm(n)%delp(is:ie,js:je,k)
       enddo
    endif

    if (id_pfnh_coarse > 0) then
       allocate(nhpres(is:ie,js:je,1:npz))
       call compute_nonhydrostatic_pressure(Atm(n), zvir, nhpres)
       call block_average(Atm(n)%gridstruct%area(is:ie,js:je), &
            nhpres(is:ie,js:je,1:npz), work_3d)
       used = send_data(id_pfnh_coarse, work_3d, Time)
       deallocate(nhpres)
    endif

    if (id_u_dt_phys_coarse > 0) then
       call block_average(Atm(n)%gridstruct%area(is:ie,js:je), &
            Atm(n)%phys_diag%phys_u_dt(is:ie,js:je,1:npz), work_3d)
       used = send_data(id_u_dt_phys_coarse, work_3d, Time)
    endif

    if (id_v_dt_phys_coarse > 0) then
       call block_average(Atm(n)%gridstruct%area(is:ie,js:je), &
            Atm(n)%phys_diag%phys_v_dt(is:ie,js:je,1:npz), work_3d)
       used = send_data(id_v_dt_phys_coarse, work_3d, Time)
    endif

    if (id_t_dt_phys_coarse > 0) then
       call block_average(mass, &
            Atm(n)%phys_diag%phys_t_dt(is:ie,js:je,1:npz), work_3d)
       used = send_data(id_t_dt_phys_coarse, work_3d, Time)
    endif

    if (id_qv_dt_phys_coarse > 0) then
       call block_average(mass, &
            Atm(n)%phys_diag%phys_qv_dt(is:ie,js:je,1:npz), work_3d)
       used = send_data(id_qv_dt_phys_coarse, work_3d, Time)
    endif

    if (id_ql_dt_phys_coarse > 0) then
       call block_average(mass, &
            Atm(n)%phys_diag%phys_ql_dt(is:ie,js:je,1:npz), work_3d)
       used = send_data(id_ql_dt_phys_coarse, work_3d, Time)
    endif

    if (id_qi_dt_phys_coarse > 0) then
       call block_average(mass, &
            Atm(n)%phys_diag%phys_qi_dt(is:ie,js:je,1:npz), work_3d)
       used = send_data(id_qi_dt_phys_coarse, work_3d, Time)
    endif

    if (id_liq_wat_dt_phys_coarse > 0) then
       call block_average(mass, &
            Atm(n)%phys_diag%phys_liq_wat_dt(is:ie,js:je,1:npz), work_3d)
       used = send_data(id_liq_wat_dt_phys_coarse, work_3d, Time)
    endif

    if (id_ice_wat_dt_phys_coarse > 0) then
       call block_average(mass, &
            Atm(n)%phys_diag%phys_ice_wat_dt(is:ie,js:je,1:npz), work_3d)
       used = send_data(id_ice_wat_dt_phys_coarse, work_3d, Time)
    endif

    if (id_qr_dt_phys_coarse > 0) then
       call block_average(mass, &
            Atm(n)%phys_diag%phys_qr_dt(is:ie,js:je,1:npz), work_3d)
       used = send_data(id_qr_dt_phys_coarse, work_3d, Time)
    endif

    if (id_qg_dt_phys_coarse > 0) then
       call block_average(mass, &
            Atm(n)%phys_diag%phys_qg_dt(is:ie,js:je,1:npz), work_3d)
       used = send_data(id_qg_dt_phys_coarse, work_3d, Time)
    endif

    if (id_qs_dt_phys_coarse > 0) then
       call block_average(mass, &
            Atm(n)%phys_diag%phys_qs_dt(is:ie,js:je,1:npz), work_3d)
       used = send_data(id_qs_dt_phys_coarse, work_3d, Time)
    endif

    if (id_u_dt_nudge_coarse > 0) then
       call block_average(Atm(n)%gridstruct%area(is:ie,js:je), &
            Atm(n)%nudge_diag%nudge_u_dt(is:ie,js:je,1:npz), work_3d)
       used = send_data(id_u_dt_nudge_coarse, work_3d, Time)
    endif

    if (id_v_dt_nudge_coarse > 0) then
       call block_average(Atm(n)%gridstruct%area(is:ie,js:je), &
            Atm(n)%nudge_diag%nudge_v_dt(is:ie,js:je,1:npz), work_3d)
       used = send_data(id_v_dt_nudge_coarse, work_3d, Time)
    endif

    if (id_t_dt_nudge_coarse > 0) then
       call block_average(mass, &
            Atm(n)%nudge_diag%nudge_t_dt(is:ie,js:je,1:npz), work_3d)
       used = send_data(id_t_dt_nudge_coarse, work_3d, Time)
    endif

    if (id_delp_dt_nudge_coarse > 0) then
       call block_average(Atm(n)%gridstruct%area(is:ie,js:je), &
            Atm(n)%nudge_diag%nudge_delp_dt(is:ie,js:je,1:npz), work_3d)
       used = send_data(id_delp_dt_nudge_coarse, work_3d, Time)
    endif

    if (id_ps_dt_nudge_coarse > 0) then
       call block_average(Atm(n)%gridstruct%area(is:ie,js:je), &
            Atm(n)%nudge_diag%nudge_ps_dt(is:ie,js:je), work_2d)
       used = send_data(id_ps_dt_nudge_coarse, work_2d, Time)
    endif

    if (id_t_dt_gfdlmp_coarse > 0) then
       call block_average(mass, &
            Atm(n)%inline_mp%T_dt(is:ie,js:je,1:npz), work_3d)
       used = send_data(id_t_dt_gfdlmp_coarse, work_3d, Time)
    endif

    if (id_qv_dt_gfdlmp_coarse > 0) then
       call block_average(mass, &
            Atm(n)%inline_mp%qv_dt(is:ie,js:je,1:npz), work_3d)
       used = send_data(id_qv_dt_gfdlmp_coarse, work_3d, Time)
    endif

    if (id_ql_dt_gfdlmp_coarse > 0) then
       call block_average(mass, &
            Atm(n)%inline_mp%ql_dt(is:ie,js:je,1:npz), work_3d)
       used = send_data(id_ql_dt_gfdlmp_coarse, work_3d, Time)
    endif

    if (id_qr_dt_gfdlmp_coarse > 0) then
       call block_average(mass, &
            Atm(n)%inline_mp%qr_dt(is:ie,js:je,1:npz), work_3d)
       used = send_data(id_qr_dt_gfdlmp_coarse, work_3d, Time)
    endif

    if (id_qi_dt_gfdlmp_coarse > 0) then
       call block_average(mass, &
            Atm(n)%inline_mp%qi_dt(is:ie,js:je,1:npz), work_3d)
       used = send_data(id_qi_dt_gfdlmp_coarse, work_3d, Time)
    endif

    if (id_ice_wat_dt_gfdlmp_coarse > 0) then
       call block_average(mass, &
            Atm(n)%inline_mp%ice_wat_dt(is:ie,js:je,1:npz), work_3d)
       used = send_data(id_ice_wat_dt_gfdlmp_coarse, work_3d, Time)
    endif

    if (id_liq_wat_dt_gfdlmp_coarse > 0) then
       call block_average(mass, &
            Atm(n)%inline_mp%liq_wat_dt(is:ie,js:je,1:npz), work_3d)
       used = send_data(id_liq_wat_dt_gfdlmp_coarse, work_3d, Time)
    endif

    if (id_qg_dt_gfdlmp_coarse > 0) then
       call block_average(mass, &
            Atm(n)%inline_mp%qg_dt(is:ie,js:je,1:npz), work_3d)
       used = send_data(id_qg_dt_gfdlmp_coarse, work_3d, Time)
    endif

    if (id_qs_dt_gfdlmp_coarse > 0) then
       call block_average(mass, &
            Atm(n)%inline_mp%qs_dt(is:ie,js:je,1:npz), work_3d)
       used = send_data(id_qs_dt_gfdlmp_coarse, work_3d, Time)
    endif

    if (id_u_dt_gfdlmp_coarse > 0) then
       call block_average(Atm(n)%gridstruct%area(is:ie,js:je), &
            Atm(n)%inline_mp%u_dt(is:ie,js:je,1:npz), work_3d)
       used = send_data(id_u_dt_gfdlmp_coarse, work_3d, Time)
    endif

    if (id_v_dt_gfdlmp_coarse > 0) then
       call block_average(Atm(n)%gridstruct%area(is:ie,js:je), &
            Atm(n)%inline_mp%v_dt(is:ie,js:je,1:npz), work_3d)
       used = send_data(id_v_dt_gfdlmp_coarse, work_3d, Time)
    endif

    if (allocated(work_2d)) deallocate(work_2d)
    if (allocated(work_3d)) deallocate(work_3d)
    if (allocated(mass)) deallocate(mass)

  end subroutine compute_coarse_grained_diagnostics

  subroutine block_average_2d(weights, fine, coarse)
    real, intent(in) :: weights(is:ie,js:je)
    real, intent(in) :: fine(is:ie,js:je)
    real, intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse)

    real, allocatable :: weighted_fine(:,:)
    integer :: i, j, i_coarse, j_coarse, a

    allocate(weighted_fine(is:ie,js:je))
    weighted_fine = weights * fine

    a = coarsening_factor - 1
    do i = is, ie, coarsening_factor
       i_coarse = (i - 1) / coarsening_factor + 1
       do j = js, je, coarsening_factor
          j_coarse = (j - 1) / coarsening_factor + 1
          coarse(i_coarse,j_coarse) = sum(weighted_fine(i:i + a,j:j + a)) / sum(weights(i:i + a,j:j + a))
       enddo
    enddo

    deallocate(weighted_fine)

  end subroutine block_average_2d

  subroutine block_average_3d_variable_2d_weights(weights, fine, coarse)
    real, intent(in) :: weights(is:ie,js:je)
    real, intent(in) :: fine(is:ie,js:je,1:npz)
    real, intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz)

    integer :: k

    do k = 1, npz
       call block_average_2d(weights, fine(is:ie,js:je,k), coarse(is_coarse:ie_coarse,js_coarse:je_coarse,k))
    enddo
  end subroutine block_average_3d_variable_2d_weights

  subroutine block_average_3d_variable_3d_weights(weights, fine, coarse)
    real, intent(in) :: weights(is:ie,js:je,1:npz)
    real, intent(in) :: fine(is:ie,js:je,1:npz)
    real, intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse,1:npz)

    integer :: k

    do k = 1, npz
       call block_average_2d(weights(is:ie,js:je,k), fine(is:ie,js:je,k), &
            coarse(is_coarse:ie_coarse,js_coarse:je_coarse,k))
    enddo
  end subroutine block_average_3d_variable_3d_weights

  ! For now we always assume that the variables passed to this function are
  ! staggered in the j dimension, but unstaggered in the i dimension.  This
  ! assumption holds for winds either on the C or D grid, which is what I
  ! anticipate this subroutine will exclusively be used for.
  subroutine block_edge_average_x_2d(dx, fine, coarse)
    real, intent(in) :: dx(is:ie,js:je+1)
    real, intent(in) :: fine(is:ie,js:je+1)
    real, intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse+1)

    integer :: i, j, i_coarse, j_coarse, a

    a = coarsening_factor - 1
    do i = is, ie, coarsening_factor
       i_coarse = (i - 1) / coarsening_factor + 1
       do j = js, je + 1, coarsening_factor
          j_coarse = (j - 1) / coarsening_factor + 1
          coarse(i_coarse,j_coarse) = sum(dx(i:i+a,j) * fine(i:i+a,j)) / sum(dx(i:i+a,j))
       enddo
    enddo
  end subroutine block_edge_average_x_2d

  ! For now we always assume that the variables passed to this function are
  ! staggered in the j dimension, but unstaggered in the i dimension.  This
  ! assumption holds for winds either on the C or D grid, which is what I
  ! anticipate this subroutine will exclusively be used for.
  subroutine block_edge_average_x_3d(dx, fine, coarse)
    real, intent(in) :: dx(is:ie,js:je+1)
    real, intent(in) :: fine(is:ie,js:je+1,1:npz)
    real, intent(out) :: coarse(is_coarse:ie_coarse,js_coarse:je_coarse+1,1:npz)

    integer :: k

    do k = 1, npz
       call block_edge_average_x_2d(dx, fine(is:ie,js:je+1,k), &
            coarse(is_coarse:ie_coarse,js_coarse:je_coarse+1,k))
    enddo
  end subroutine block_edge_average_x_3d

  ! For now we always assume that the variables passed to this function are
  ! staggered in the i dimension, but unstaggered in the j dimension.  This
  ! assumption holds for winds either on the C or D grid, which is what I
  ! anticipate this subroutine will exclusively be used for.
  subroutine block_edge_average_y_2d(dy, fine, coarse)
    real, intent(in) :: dy(is:ie+1,js:je)
    real, intent(in) :: fine(is:ie+1,js:je)
    real, intent(out) :: coarse(is_coarse:ie_coarse+1,js_coarse:je_coarse)

    integer :: i, j, i_coarse, j_coarse, a

    a = coarsening_factor - 1
    do i = is, ie + 1, coarsening_factor
       i_coarse = (i - 1) / coarsening_factor + 1
       do j = js, je, coarsening_factor
          j_coarse = (j - 1) / coarsening_factor + 1
          coarse(i_coarse,j_coarse) = sum(dy(i,j:j+a) * fine(i,j:j+a)) / sum(dy(i,j:j+a))
       enddo
    enddo
  end subroutine block_edge_average_y_2d

  ! For now we always assume that the variables passed to this function are
  ! staggered in the i dimension, but unstaggered in the j dimension.  This
  ! assumption holds for winds either on the C or D grid, which is what I
  ! anticipate this subroutine will exclusively be used for.
  subroutine block_edge_average_y_3d(dy, fine, coarse)
    real, intent(in) :: dy(is:ie+1,js:je)
    real, intent(in) :: fine(is:ie+1,js:je,1:npz)
    real, intent(out) :: coarse(is_coarse:ie_coarse+1,js_coarse:je_coarse,1:npz)

    integer :: k

    do k = 1, npz
       call block_edge_average_y_2d(dy, fine(is:ie+1,js:je,k), &
            coarse(is_coarse:ie_coarse+1,js_coarse:je_coarse,k))
    enddo
  end subroutine block_edge_average_y_3d

  ! This subroutine is copied from fms/horiz_interp/horiz_interp_test.F90;
  ! domain_decomp in fv_mp_mod.F90 does something similar, but it does a
  ! few other unnecessary things (and requires more arguments).
  subroutine define_cubic_mosaic(domain, ni, nj, layout)
    type(domain2d), intent(inout) :: domain
    integer,        intent(in)    :: layout(:)
    integer,        intent(in)    :: ni, nj
    integer   :: pe_start(6), pe_end(6)
    integer   :: global_indices(4,6), layout2d(2,6)
    integer, dimension(12)        :: istart1, iend1, jstart1, jend1, tile1
    integer, dimension(12)        :: istart2, iend2, jstart2, jend2, tile2
    integer                       :: ntiles, num_contact
    integer :: p, npes_per_tile, i

    ntiles = 6
    num_contact = 12
    p = 0
    npes_per_tile = mpp_npes()/ntiles
    do i = 1, 6
       layout2d(:,i) = layout(:)
       global_indices(1,i) = 1
       global_indices(2,i) = ni
       global_indices(3,i) = 1
       global_indices(4,i) = nj
       pe_start(i) = p
       p = p + npes_per_tile
       pe_end(i) = p-1
    enddo

    !--- Contact line 1, between tile 1 (EAST) and tile 2 (WEST)
    tile1(1) = 1;     tile2(1) = 2
    istart1(1) = ni;  iend1(1) = ni;  jstart1(1) = 1;      jend1(1) = nj
    istart2(1) = 1;   iend2(1) = 1;   jstart2(1) = 1;      jend2(1) = nj
    !--- Contact line 2, between tile 1 (NORTH) and tile 3 (WEST)
    tile1(2) = 1; tile2(2) = 3
    istart1(2) = 1;      iend1(2) = ni;  jstart1(2) = nj;  jend1(2) = nj
    istart2(2) = 1;      iend2(2) = 1;   jstart2(2) = nj;  jend2(2) = 1
    !--- Contact line 3, between tile 1 (WEST) and tile 5 (NORTH)
    tile1(3) = 1;     tile2(3) = 5
    istart1(3) = 1;   iend1(3) = 1;      jstart1(3) = 1;   jend1(3) = nj
    istart2(3) = ni;  iend2(3) = 1;      jstart2(3) = nj;  jend2(3) = nj
    !--- Contact line 4, between tile 1 (SOUTH) and tile 6 (NORTH)
    tile1(4) = 1; tile2(4) = 6
    istart1(4) = 1;      iend1(4) = ni;  jstart1(4) = 1;   jend1(4) = 1
    istart2(4) = 1;      iend2(4) = ni;  jstart2(4) = nj;  jend2(4) = nj
    !--- Contact line 5, between tile 2 (NORTH) and tile 3 (SOUTH)
    tile1(5) = 2;        tile2(5) = 3
    istart1(5) = 1;      iend1(5) = ni;  jstart1(5) = nj;  jend1(5) = nj
    istart2(5) = 1;      iend2(5) = ni;  jstart2(5) = 1;   jend2(5) = 1
    !--- Contact line 6, between tile 2 (EAST) and tile 4 (SOUTH)
    tile1(6) = 2; tile2(6) = 4
    istart1(6) = ni;  iend1(6) = ni;  jstart1(6) = 1;      jend1(6) = nj
    istart2(6) = ni;  iend2(6) = 1;   jstart2(6) = 1;      jend2(6) = 1
    !--- Contact line 7, between tile 2 (SOUTH) and tile 6 (EAST)
    tile1(7) = 2; tile2(7) = 6
    istart1(7) = 1;   iend1(7) = ni;  jstart1(7) = 1;   jend1(7) = 1
    istart2(7) = ni;  iend2(7) = ni;  jstart2(7) = nj;  jend2(7) = 1
    !--- Contact line 8, between tile 3 (EAST) and tile 4 (WEST)
    tile1(8) = 3; tile2(8) = 4
    istart1(8) = ni;  iend1(8) = ni;  jstart1(8) = 1;      jend1(8) = nj
    istart2(8) = 1;   iend2(8) = 1;   jstart2(8) = 1;      jend2(8) = nj
    !--- Contact line 9, between tile 3 (NORTH) and tile 5 (WEST)
    tile1(9) = 3; tile2(9) = 5
    istart1(9) = 1;      iend1(9) = ni;  jstart1(9) = nj;  jend1(9) = nj
    istart2(9) = 1;      iend2(9) = 1;   jstart2(9) = nj;  jend2(9) = 1
    !--- Contact line 10, between tile 4 (NORTH) and tile 5 (SOUTH)
    tile1(10) = 4; tile2(10) = 5
    istart1(10) = 1;     iend1(10) = ni; jstart1(10) = nj; jend1(10) = nj
    istart2(10) = 1;     iend2(10) = ni; jstart2(10) = 1;  jend2(10) = 1
    !--- Contact line 11, between tile 4 (EAST) and tile 6 (SOUTH)
    tile1(11) = 4; tile2(11) = 6
    istart1(11) = ni; iend1(11) = ni; jstart1(11) = 1;     jend1(11) = nj
    istart2(11) = ni; iend2(11) = 1;  jstart2(11) = 1;     jend2(11) = 1
    !--- Contact line 12, between tile 5 (EAST) and tile 6 (WEST)
    tile1(12) = 5; tile2(12) = 6
    istart1(12) = ni; iend1(12) = ni; jstart1(12) = 1;     jend1(12) = nj
    istart2(12) = 1;  iend2(12) = 1;  jstart2(12) = 1;     jend2(12) = nj
    call mpp_define_mosaic(global_indices, layout2d, domain, ntiles, num_contact, tile1, tile2, &
         istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
         pe_start, pe_end, symmetry = .true., name = 'coarse cubic mosaic'  )

    return

  end subroutine define_cubic_mosaic

  subroutine compute_nonhydrostatic_pressure(Atm, zvir, nhpres)
    type(fv_atmos_type), intent(in) :: Atm
    real, intent(in) :: zvir
    real, intent(out) :: nhpres(is:ie,js:je,1:npz)

    integer ::  i, j, k
    integer :: sphum

    sphum = get_tracer_index(MODEL_ATMOS, 'sphum')

#ifdef GFS_PHYS
    do k= 1, npz
       do j=js, je
          do i=is, ie
             nhpres(i,j,k) = Atm%delp(i,j,k)*(1.-sum(Atm%q(i,j,k,2:Atm%flagstruct%nwat)))
          enddo
       enddo
    enddo

    do k = 1, npz
       do j = js, je
          do i = is, ie
             nhpres(i,j,k) = -nhpres(i,j,k)/(Atm%delz(i,j,k)*grav)*rdgas*    &
                  Atm%pt(i,j,k)*(1.+zvir*Atm%q(i,j,k,sphum))
          enddo
       enddo
    enddo
#else
      do k = 1, npz
         do j = js, je
            do i = is, ie
               nhpres(i,j,k) = -Atm%delp(i,j,k)/(Atm%delz(i,j,k)*grav)*rdgas*  &
                    Atm%pt(i,j,k)*(1.+zvir*Atm%q(i,j,k,sphum))
            enddo
         enddo
      enddo
#endif
  end subroutine compute_nonhydrostatic_pressure

end module fv_coarse_graining_mod
