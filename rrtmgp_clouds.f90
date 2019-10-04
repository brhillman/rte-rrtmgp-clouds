
subroutine stop_on_err(error_msg)
  use iso_fortran_env, only : error_unit
  character(len=*), intent(in) :: error_msg

  if(error_msg /= "") then
    write (error_unit,*) trim(error_msg)
    write (error_unit,*) "rrtmgp_garand_atmos stopping"
    stop
  end if
end subroutine stop_on_err

program rte_rrtmgp_clouds
  use mo_rte_kind,           only: wp
  use mo_optical_props,      only: ty_optical_props, &
                                   ty_optical_props_arry, ty_optical_props_1scl, ty_optical_props_2str
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
  use mo_cloud_optics,       only: ty_cloud_optics
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_source_functions,   only: ty_source_func_lw
  use mo_fluxes,             only: ty_fluxes_broadband
  use mo_rte_lw,             only: rte_lw
  use mo_rte_sw,             only: rte_sw
  use mo_load_coefficients,  only: load_and_init
  use mo_load_cloud_coefficients, &
                             only: load_cld_lutcoeff, load_cld_padecoeff
  use mo_garand_atmos_io,    only: read_atmos, write_lw_fluxes, write_sw_fluxes
  implicit none
  ! ----------------------------------------------------------------------------------
  ! Variables
  ! ----------------------------------------------------------------------------------
  ! Arrays: dimensions (col, lay)
  real(wp), dimension(:,:),   allocatable :: p_lay, t_lay, p_lev
  real(wp), dimension(:,:),   allocatable :: col_dry
  real(wp), dimension(:,:),   allocatable, target :: temp_array
  real(wp), dimension(:),     pointer     :: temp_vec
  !
  ! Longwave only
  !
  real(wp), dimension(:,:),   allocatable :: t_lev
  real(wp), dimension(:),     allocatable :: t_sfc
  real(wp), dimension(:,:),   allocatable :: emis_sfc ! First dimension is band
  !
  ! Shortwave only
  !
  real(wp), dimension(:),     allocatable :: mu0
  real(wp), dimension(:,:),   allocatable :: sfc_alb_dir, sfc_alb_dif ! First dimension is band
  !
  ! Source functions
  !
  !   Longwave
  type(ty_source_func_lw)               :: lw_sources
  !   Shortwave
  real(wp), dimension(:,:), allocatable :: toa_flux
  !
  ! Clouds
  !
  real(wp), allocatable, dimension(:,:) :: lwp, iwp, rel, rei
  logical,  allocatable, dimension(:,:) :: cloud_mask
  !
  ! Output variables
  !
  real(wp), dimension(:,:), target, &
                            allocatable :: flux_up, flux_dn, flux_dir
  !
  ! Derived types from the RTE and RRTMGP libraries
  !
  type(ty_gas_optics_rrtmgp) :: k_dist
  type(ty_cloud_optics)      :: cloud_optics
  type(ty_gas_concs)         :: gas_concs, gas_concs_garand, gas_concs_1col
  class(ty_optical_props_arry), &
                 allocatable :: atmos, clouds
  type(ty_fluxes_broadband)  :: fluxes

  !
  ! Inputs to RRTMGP
  !
  logical :: top_at_1, is_sw, is_lw

  integer :: ncol, nlay, nbnd, ngpt, nUserArgs=0
  character(len=6) :: ncol_char

  character(len=256) :: input_file, k_dist_file, cloud_optics_file
  ! ----------------------------------------------------------------------------------
  ! Code
  ! ----------------------------------------------------------------------------------
  !
  ! Parse command line for any file names, block size
  !
  ! rrtmgp_clouds rrtmgp-clouds.nc $RRTMGP_ROOT/rrtmgp/data/rrtmgp-data-lw-g256-2018-12-04.nc $RRTMGP_ROOT/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-lw.nc  128
  ! rrtmgp_clouds rrtmgp-clouds.nc $RRTMGP_ROOT/rrtmgp/data/rrtmgp-data-sw-g224-2018-12-04.nc $RRTMGP_ROOT/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-sw.nc  128
  nUserArgs = command_argument_count()
  if (nUserArgs <  4) call stop_on_err("Need to supply input_file k_distribution_file ncol.")
  if (nUserArgs >= 1) call get_command_argument(1,input_file)
  if (nUserArgs >= 2) call get_command_argument(2,k_dist_file)
  if (nUserArgs >= 3) call get_command_argument(3,cloud_optics_file)
  if (nUserArgs >= 4) then
    call get_command_argument(4, ncol_char)
    read(ncol_char, '(i6)') ncol
    if(ncol <= 0) call stop_on_err("Specify positive ncol.")
  end if
  if (nUserArgs >  5) print *, "Ignoring command line arguments beyond the first three..."
  if(trim(input_file) == '-h' .or. trim(input_file) == "--help") then
    call stop_on_err("rrtmgp_clouds input_file absorption_coefficients_file cloud_optics_file ncol")
  end if
  !
  ! Read temperature, pressure, gas concentrations.
  !   Arrays are allocated as they are read
  !
  call read_atmos(input_file,                 &
                  p_lay, t_lay, p_lev, t_lev, &
                  gas_concs_garand, col_dry)
  deallocate(col_dry)
  nlay = size(p_lay, 2)
  ! For clouds we'll use the first column, repeated over and over
  !   These shenanigans are so we supply the gases as vectors, considered constant over the column dimension
  allocate(temp_array(size(p_lay, 1), nlay))
  temp_vec => temp_array(1,:)
  call stop_on_err(gas_concs_garand%get_vmr('h2o', temp_array))
  call stop_on_err(gas_concs%set_vmr       ('h2o', temp_vec))
  call stop_on_err(gas_concs_garand%get_vmr('co2', temp_array))
  call stop_on_err(gas_concs%set_vmr       ('co2', temp_vec))
  call stop_on_err(gas_concs_garand%get_vmr('o3' , temp_array))
  call stop_on_err(gas_concs%set_vmr       ('o3' , temp_vec))
  call stop_on_err(gas_concs_garand%get_vmr('n2o', temp_array))
  call stop_on_err(gas_concs%set_vmr       ('n2o', temp_vec))
  call stop_on_err(gas_concs_garand%get_vmr('co' , temp_array))
  call stop_on_err(gas_concs%set_vmr       ('co' , temp_vec))
  call stop_on_err(gas_concs_garand%get_vmr('ch4', temp_array))
  call stop_on_err(gas_concs%set_vmr       ('ch4', temp_vec))
  call stop_on_err(gas_concs_garand%get_vmr('o2' , temp_array))
  call stop_on_err(gas_concs%set_vmr       ('o2' , temp_vec))
  call stop_on_err(gas_concs_garand%get_vmr('n2' , temp_array))
  call stop_on_err(gas_concs%set_vmr       ('n2' , temp_vec))
  deallocate(temp_array); nullify(temp_vec)
  !  If we trusted in Fortran allocate-on-assign we could skip the temp_array here
  allocate(temp_array(ncol, nlay))
  temp_array = spread(p_lay(1,:), dim = 1, ncopies=ncol)
  call move_alloc(temp_array, p_lay)
  allocate(temp_array(ncol, nlay))
  temp_array = spread(t_lay(1,:), dim = 1, ncopies=ncol)
  call move_alloc(temp_array, t_lay)
  allocate(temp_array(ncol, nlay+1))
  temp_array = spread(p_lev(1,:), dim = 1, ncopies=ncol)
  call move_alloc(temp_array, p_lev)
  allocate(temp_array(ncol, nlay+1))
  temp_array = spread(t_lev(1,:), dim = 1, ncopies=ncol)
  call move_alloc(temp_array, t_lev)

  ! ----------------------------------------------------------------------------
  ! load data into classes
  call load_and_init(k_dist, k_dist_file, gas_concs)
  is_sw = k_dist%source_is_external()
  is_lw = .not. is_sw
  !
  ! Should also try with Pade calculations
  !  call load_cld_padecoeff(cloud_optics, cloud_optics_file)
  !
  call load_cld_lutcoeff(cloud_optics, cloud_optics_file)
  ! integer division to find a nominal ice roughness
  call stop_on_err(cloud_optics%set_ice_roughness(2))
  ! ----------------------------------------------------------------------------
  !
  ! Problem sizes
  !
  nbnd = k_dist%get_nband()
  ngpt = k_dist%get_ngpt()
  top_at_1 = p_lay(1, 1) < p_lay(1, nlay)

  ! ----------------------------------------------------------------------------
  ! LW calculations neglect scattering; SW calculations use the 2-stream approximation
  !   Here we choose the right variant of optical_props.
  !
  if(is_sw) then
    allocate(ty_optical_props_2str::atmos)
    allocate(ty_optical_props_2str::clouds)
  else
    allocate(ty_optical_props_1scl::atmos)
    allocate(ty_optical_props_1scl::clouds)
  end if
  ! Clouds optical props are defined by band
  call stop_on_err(clouds%init(k_dist%get_band_lims_wavenumber()))
  !
  ! Allocate arrays for the optical properties themselves.
  !
  select type(atmos)
    class is (ty_optical_props_1scl)
      call stop_on_err(atmos%alloc_1scl(ncol, nlay, k_dist))
    class is (ty_optical_props_2str)
      call stop_on_err(atmos%alloc_2str( ncol, nlay, k_dist))
    class default
      call stop_on_err("rrtmgp_garand_atmos: Don't recognize the kind of optical properties ")
  end select
  select type(clouds)
    class is (ty_optical_props_1scl)
      call stop_on_err(clouds%alloc_1scl(ncol, nlay))
    class is (ty_optical_props_2str)
      call stop_on_err(clouds%alloc_2str(ncol, nlay))
    class default
      call stop_on_err("rrtmgp_garand_atmos: Don't recognize the kind of optical properties ")
  end select
  ! ----------------------------------------------------------------------------
  !  Boundary conditions depending on whether the k-distribution being supplied
  !   is LW or SW
  if(is_sw) then
    allocate(toa_flux(ncol, ngpt), sfc_alb_dir(nbnd, ncol), sfc_alb_dif(nbnd, ncol), mu0(ncol))
    ! Ocean-ish values for no particular reason
    sfc_alb_dir = 0.06_wp
    sfc_alb_dif = 0.06_wp
    mu0 = .86_wp
  else
    call stop_on_err(lw_sources%alloc(ncol, nlay, k_dist))
    allocate(t_sfc(ncol), emis_sfc(nbnd, ncol))
    ! Surface temperature
    t_sfc = t_lev(1, merge(nlay+1, 1, top_at_1))
    emis_sfc = 0.98_wp
  end if
  ! ----------------------------------------------------------------------------
  !
  ! Fluxes
  !
  allocate(flux_up(ncol,nlay+1), flux_dn(ncol,nlay+1))
  fluxes%flux_up => flux_up
  fluxes%flux_dn => flux_dn
  if(is_sw) then
    allocate(flux_dir(ncol,nlay+1))
    fluxes%flux_dn_dir => flux_dir
  end if
  !
  ! Clouds
  !
  allocate(lwp(ncol,nlay), iwp(ncol,nlay), &
           rel(ncol,nlay), rei(ncol,nlay), cloud_mask(ncol,nlay))
  ! Restrict clouds to troposphere (< 100 hPa = 100*100 Pa)
  !   and not very close to the ground
  cloud_mask = p_lay < 100._wp * 100._wp .and. p_lay > 900._wp
  !
  ! Ice and liquid will overlap in a few layers
  !
  lwp = merge(10._wp, 0._wp, cloud_mask .and. t_lay > 263._wp)
  iwp = merge(10._wp, 0._wp, cloud_mask .and. t_lay < 273._wp)
  rel = merge(0.5 * (cloud_optics%get_min_radius_liq() + cloud_optics%get_max_radius_liq()), &
              0._wp, lwp > 0._wp)
  rei = merge(0.5 * (cloud_optics%get_min_radius_ice() + cloud_optics%get_max_radius_ice()), &
              0._wp, iwp > 0._wp)
  ! ----------------------------------------------------------------------------
  !
  ! All work from here to the writing of flxues should happen on the GPU
  !
  call stop_on_err(                                      &
    cloud_optics%cloud_optics(lwp, iwp, rel, rei, clouds))
  !
  ! Solvers
  !
  if(is_lw) then
    call stop_on_err(k_dist%gas_optics(p_lay, p_lev, &
                                       t_lay, t_sfc, &
                                       gas_concs,    &
                                       atmos,        &
                                       lw_sources,   &
                                       tlev = t_lev))
    call stop_on_err(clouds%increment(atmos))
    call stop_on_err(rte_lw(atmos, top_at_1, &
                            lw_sources,      &
                            emis_sfc,        &
                            fluxes))
    call write_lw_fluxes(input_file, flux_up, flux_dn)
  else
    call stop_on_err(k_dist%gas_optics(p_lay, p_lev, &
                                       t_lay,        &
                                       gas_concs,    &
                                       atmos,        &
                                       toa_flux))
    call stop_on_err(clouds%delta_scale())
    call stop_on_err(clouds%increment(atmos))
    call stop_on_err(rte_sw(atmos, top_at_1, &
                            mu0,   toa_flux, &
                            sfc_alb_dir, sfc_alb_dif, &
                            fluxes))
  call write_sw_fluxes(input_file, flux_up, flux_dn, flux_dir)
  end if

end program rte_rrtmgp_clouds
