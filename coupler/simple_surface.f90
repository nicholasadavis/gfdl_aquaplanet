module simple_surface_mod

use atmos_model_mod, only: atmos_data_type

use surface_flux_mod, only: surface_flux

use diag_integral_mod, only: diag_integral_field_init, &
                             sum_diag_integral_field

use           fms_mod, only: file_exist, open_namelist_file, check_nml_error, &
                             error_mesg, FATAL, mpp_pe, mpp_root_pe, &
                             open_file, close_file, read_data, write_data, &
                             write_version_number, stdlog, set_domain, read_data, NOTE, &
                             mpp_npes

use  diag_manager_mod, only: register_diag_field,  &
                             register_static_field, send_data

use  time_manager_mod, only: time_type

use     constants_mod, only: rdgas, rvgas, cp_air, hlv, hlf, cp_ocean, radius, omega

use ocean_rough_mod, only: compute_ocean_roughness

use transforms_mod, only:  get_grid_domain, grid_domain, get_deg_lat, get_sin_lat, &
                           get_cos_lat, get_lat_max

use spherical_mod, only: compute_lat_deriv_cos

use mpp_domains_reduce_mod, only: mpp_global_field

use mpp_util_mod, only: mpp_sync

implicit none
private

public :: simple_surface_init,   &
          compute_flux,          &
          update_simple_surface, &
          simple_surface_end

!-----------------------------------------------------------------------
character(len=128) :: version = '$Id: simple_surface.f90,v 1.1.2.3 2005/05/21 02:02:04 pjp Exp $'
character(len=128) :: tagname = '$Name:  $'

!-----------------------------------------------------------------------
!-------- namelist (for diagnostics) ------

character(len=14), parameter :: mod_name = 'simple_surface'

integer :: id_drag_moist,  id_drag_heat,  id_drag_mom,              &
           id_rough_moist, id_rough_heat, id_rough_mom,             &
           id_u_star, id_b_star, id_u_flux, id_v_flux, id_t_surf,   &
           id_t_flux, id_q_flux, id_o_flux, id_r_flux,              &
           id_t_atm,  id_u_atm,  id_v_atm,  id_wind,                &
           id_t_ref,  id_rh_ref, id_u_ref,  id_v_ref,               &
           id_del_h,  id_del_m,  id_del_q, id_albedo,               &
           id_entrop_evap, id_entrop_shflx, id_entrop_lwflx,        &
           id_eddy_u, id_eddy_v, id_eddy_t, id_eddy_q,              &
           id_o_flux_mer, is, ie, js, je, lat_max

logical :: first_static = .true.
logical :: do_init = .true.

!-----------------------------------------------------------------------

real ::  z_ref_heat      = 2.,       &
         z_ref_mom       = 10.,      &
          heat_capacity   = 1.e07,    &
          const_roughness = 3.21e-05, &
          const_albedo    = 0.12,     &
    max_of          = 25.,      &
    lonmax_of       = 180.,     &
    latmax_of       = 0.,       &
          latwidth_of     = 15.,      &
    lonwidth_of     = 90.,      &
    higher_albedo    = 0.38,    &
    lat_glacier      = 45.,     &
    maxofmerid       = .5,      &
    latmaxofmerid    = 25.,     &
    Tm               = 305.,    &
    deltaT           = 40.
          

integer :: surface_choice   = 1
integer :: roughness_choice = 1
integer :: albedo_choice    = 1
logical :: do_oflx          = .false.
logical :: do_oflxmerid     = .false.
logical :: apply_fluxes     = .false.
logical :: do_sst_clamp     = .false.
character(len=128) :: input_file = 'fluxes.nc'
character(len=128) :: sst_file = 'sst.nc'

namelist /simple_surface_nml/ z_ref_heat, z_ref_mom,             &
                              surface_choice,  heat_capacity,    &
                              roughness_choice, const_roughness, &
                              albedo_choice, const_albedo, do_oflx, &
                              max_of, lonmax_of, latmax_of, latwidth_of, &
                              lonwidth_of, higher_albedo, lat_glacier, &
                              do_oflxmerid, maxofmerid, latmaxofmerid, Tm, &
                              deltaT, do_sst_clamp, apply_fluxes, input_file, sst_file

!-----------------------------------------------------------------------

!---- allocatable module storage ------

  real, allocatable, dimension(:,:) :: e_t_n, f_t_delt_n, &
                                       e_q_n, f_q_delt_n   

  real, allocatable, dimension(:,:) :: dhdt_surf, dedt_surf, dedq_surf, &
                                       drdt_surf, dhdt_atm, dedq_atm, &
                                       flux_t, flux_q, flux_lw

  real, allocatable, dimension(:,:) :: sst, flux_u, flux_v, flux_o_global, &
                                       flux_meridional_global, f_global, cos_lat_global, &
                                       delta_t_surf_zm_global, delta_lat_global, lat_global, &
                                       delta_flux_global, cos_deriv_global, tau_x_zm_global, &
                                       t_surf_zm_global, enthalpy_global, flux_u_zm_global, &
                                       sin_lat_global, sst_zm_global, p

  real, allocatable, dimension(:)   :: flux_o_global_vector

contains

!#######################################################################

 subroutine compute_flux     ( dt, Time, Atm, land_frac,            &
                               t_surf_atm, albedo, rough_mom,       &
                               flux_u_atm, flux_v_atm, dtaudv_atm,  &
                               u_star, b_star                       )

 real,                   intent(in)  :: dt
 type       (time_type), intent(in)  :: Time
 type (atmos_data_type), intent(in)  :: Atm
 real, dimension(:,:),   intent(out) :: albedo,    rough_mom,       &
                                        land_frac, dtaudv_atm,      &
                                        flux_u_atm, flux_v_atm,     &
                                        u_star, b_star


real, dimension(:,:),   intent(out) :: t_surf_atm
       
real, dimension(size(Atm%t_bot,1), size(Atm%t_bot,2)) :: &
       u_surf, v_surf, rough_heat, rough_moist,          &
       stomatal, snow, water, max_water,                 &
       q_star, q_surf, cd_q, cd_t, cd_m, wind, dtaudu_atm, &
       u_stress_eddy, v_stress_eddy,delta_t,work,sst_clamp,  &
       zm_speed_spread, drag_m, rho_drag, tv_atm, rho, &
       e_sat, work_diff, rho_t, rho_q, drag_t, drag_q, p_ratio, th_atm, q_sat, t_eddy, q_eddy, work_2

real, dimension(size(Atm%t_bot,1))  :: zm_u_vect, zm_v_vect, eddy_u, eddy_v, work_zm

real, dimension(size(Atm%t_bot,2))  :: interp_field_v, interp_field_u, deg_lat, &
                                       t_bot_zm, q_bot_zm, u_bot_zm, v_bot_zm, &
                                       p_bot_zm, z_bot_zm, p_surf_zm, t_surf_atm_zm, &
                                       q_surf_zm, u_surf_zm, v_surf_zm, rough_mom_zm, &
                                       rough_heat_zm, rough_moist_zm, gust_zm, &
                                       flux_t_zm, flux_q_zm, flux_lw_zm, flux_u_zm, &
                                       flux_v_zm, cd_m_zm, cd_t_zm, cd_q_zm, wind_zm, &
                                       u_star_zm, b_star_zm, q_star_zm, dhdt_surf_zm, &
                                       dedt_surf_zm, dedq_surf_zm, drdt_surf_zm, &
                                       dhdt_atm_zm, dedq_atm_zm, dtaudu_atm_zm, &
                                       dtaudv_atm_zm, t_eddy_zm, &
                                       v_eddy_zm, u_eddy_zm, q_eddy_zm

real                                :: tau_applied = 5.
       
logical, dimension(size(Atm%t_bot,1), size(Atm%t_bot,2)) :: &
       mask, glacier, seawater

logical, dimension(size(Atm%t_bot,2)) :: mask_zm, seawater_zm

logical :: used

integer :: j, dummy, i
real :: lat, pi, rand_num_lat, zm_v, zm_u, work_mean, zm_speed

real, parameter :: d622   = rdgas/rvgas
real, parameter :: d378   = 1.-d622
real            :: d608   = d378/d622
real            :: kappa  = rdgas/cp_air
character(len=64) :: ch1
call set_domain(grid_domain)
pi = 4.0*atan(1.)



!-----------------------------------------------------------------------

   if (do_init) call error_mesg ('compute_flux',  &
                 'must call simple_surface_init first', FATAL)

!-----------------------------------------------------------------------
!------ allocate storage also needed in flux_up_to_atmos -----

allocate (e_t_n      (size(Atm%t_bot,1), size(Atm%t_bot,2)), &
         e_q_n      (size(Atm%t_bot,1), size(Atm%t_bot,2)), &
         f_t_delt_n (size(Atm%t_bot,1), size(Atm%t_bot,2)), &
         f_q_delt_n (size(Atm%t_bot,1), size(Atm%t_bot,2)), &
         dhdt_surf  (size(Atm%t_bot,1), size(Atm%t_bot,2)), &
         dedt_surf  (size(Atm%t_bot,1), size(Atm%t_bot,2)), &
         dedq_surf  (size(Atm%t_bot,1), size(Atm%t_bot,2)), &
         drdt_surf  (size(Atm%t_bot,1), size(Atm%t_bot,2)), &
         dhdt_atm   (size(Atm%t_bot,1), size(Atm%t_bot,2)), &
         dedq_atm   (size(Atm%t_bot,1), size(Atm%t_bot,2)), &
         flux_t     (size(Atm%t_bot,1), size(Atm%t_bot,2)), &
         flux_q     (size(Atm%t_bot,1), size(Atm%t_bot,2)), &
        ! flux_u     (size(Atm%t_bot,1), size(Atm%t_bot,2)), &
         flux_lw    (size(Atm%t_bot,1), size(Atm%t_bot,2))  )

u_surf     = 0.0
v_surf     = 0.0
stomatal   = 0.0
snow       = 0.0
water      = 1.0
max_water  = 1.0
dummy      = 1
mask    = .true.
glacier = .false.
seawater= .false. 
   
if(roughness_choice == 1) then
 rough_mom   = const_roughness
 rough_heat  = const_roughness
 rough_moist = const_roughness
elseif(roughness_choice == 2) then
 call compute_ocean_roughness (mask, u_star, &           ! Inchon and beyond. Changes answers.
                     rough_mom, rough_heat, rough_moist)
endif

if(albedo_choice == 1) then
 albedo = const_albedo
elseif(albedo_choice == 2) then
 do j = 1, size(Atm%t_bot,2)
   lat = 0.5*(Atm%lat_bnd(j+1) + Atm%lat_bnd(j))*180/pi

   if ( lat > lat_glacier ) then

     albedo(:,j) = higher_albedo

   else

     albedo(:,j) = const_albedo

   endif

 enddo

endif
   
delta_t = Atm%Surf_Diff%delta_t
cd_t = 0.0
cd_m = 0.0
cd_q = 0.0

if (apply_fluxes) then
  call set_domain(grid_domain)
  call read_data(input_file, 'forcing_taux',   u_stress_eddy(:,:), grid_domain)
  call read_data(input_file, 'forcing_tauy',   v_stress_eddy(:,:), grid_domain)
  call read_data(input_file, 'forcing_shf',   t_eddy(:,:), grid_domain)
  call read_data(input_file, 'forcing_evap',   q_eddy(:,:), grid_domain)
else
  u_stress_eddy = 0
  v_stress_eddy = 0
  t_eddy        = 0
  q_eddy        = 0
endif

if (do_sst_clamp) then 
  call set_domain(grid_domain)
  call read_data(sst_file, 'ts_field',   sst_clamp(:,:), grid_domain)
  t_surf_atm = sst_clamp
else
  t_surf_atm = sst
endif

call surface_flux (Atm%t_bot, Atm%q_bot, Atm%u_bot, Atm%v_bot,        & ! Lima
                  Atm%p_bot, Atm%z_bot, Atm%p_surf, t_surf_atm,      &
                  t_surf_atm,                                        & ! Required argument, intent(in). t_surf_atm instead of Land%t_ca
                  q_surf, u_surf, v_surf,                            &
                  rough_mom, rough_heat, rough_moist,                &
                  rough_mom,                                         & ! Required argument, intent(in). rough_mom instead of Land%rough_scale
                  Atm%gust, flux_t, flux_q, flux_lw, flux_u, flux_v, &
                  cd_m, cd_t, cd_q, wind, u_star, b_star, q_star,    &
                  dhdt_surf, dedt_surf,                              &
                  dedq_surf,                                         & ! Required argument, intent(out), but not needed by this model.
                  drdt_surf, dhdt_atm, dedq_atm,                     &
                  dtaudu_atm,                                        & ! Required argument, intent(out), but not needed by this model.
                  dtaudv_atm, dt,                                    & ! Required argument, intent(in). Looks like it should be .false. everywhere.
                  .not.mask,                                         &
                  seawater,                                          & ! Required argument, intent(in). Looks like fudgefactor for salt water. Use .false.
                  mask, u_stress_eddy, v_stress_eddy, t_eddy, q_eddy)                                               ! Required argument, intent(in). Looks like it should be .true. everywhere.

if ( id_eddy_t > 0 .or. id_eddy_q > 0 .or. id_eddy_u > 0 .or. id_eddy_v > 0) then
  u_surf_zm = 0.0
  v_surf_zm = 0.0
  cd_t_zm = 0.0
  cd_m_zm = 0.0
  cd_q_zm = 0.0
  t_eddy_zm = 0.0
  q_eddy_zm = 0.0
  u_eddy_zm = 0.0
  v_eddy_zm = 0.0

  call zonal_mean(Atm%t_bot,t_bot_zm)
  call zonal_mean(Atm%q_bot,q_bot_zm)
  call zonal_mean(Atm%u_bot,u_bot_zm)
  call zonal_mean(Atm%v_bot,v_bot_zm)
  call zonal_mean(Atm%p_bot,p_bot_zm)
  call zonal_mean(Atm%p_surf,p_surf_zm)
  call zonal_mean(Atm%z_bot,z_bot_zm)
  call zonal_mean(t_surf_atm,t_surf_atm_zm)
  !no real way to do these for mixed land/sea stuff
  rough_mom_zm   = const_roughness
  rough_heat_zm  = const_roughness
  rough_moist_zm = const_roughness
  call zonal_mean(Atm%gust,gust_zm)
  mask_zm = .true.
  seawater_zm = .false.

  call surface_flux (t_bot_zm, q_bot_zm, u_bot_zm, v_bot_zm, &
                    p_bot_zm, z_bot_zm, p_surf_zm, t_surf_atm_zm, &
                    t_surf_atm_zm, q_surf_zm, u_surf_zm, v_surf_zm, &
                    rough_mom_zm, rough_heat_zm, rough_moist_zm, &
                    rough_mom_zm, gust_zm, flux_t_zm, flux_q_zm, flux_lw_zm, &
                    flux_u_zm, flux_v_zm, cd_m_zm, cd_t_zm, cd_q_zm, wind_zm, u_star_zm, &
                    b_star_zm, q_star_zm, dhdt_surf_zm, dedt_surf_zm, &
                    dedq_surf_zm,drdt_surf_zm, dhdt_atm_zm, dedq_atm_zm, &
                    dtaudu_atm_zm,dtaudv_atm_zm, dt, .not.mask_zm,seawater_zm,&
                    mask_zm, u_eddy_zm, v_eddy_zm, t_eddy_zm, q_eddy_zm)
endif

flux_u_atm = flux_u 
flux_v_atm = flux_v 

land_frac = 0.0

!=======================================================================
!-------------------- diagnostics section ------------------------------

if ( id_wind       > 0 ) used = send_data ( id_wind,       wind,               Time )
if ( id_drag_moist > 0 ) used = send_data ( id_drag_moist, cd_q,               Time )
if ( id_drag_heat  > 0 ) used = send_data ( id_drag_heat,  cd_t,               Time )
if ( id_drag_mom   > 0 ) used = send_data ( id_drag_mom,   cd_m,               Time )
if ( id_rough_heat > 0 ) used = send_data ( id_rough_heat, rough_heat,         Time )
if ( id_rough_mom  > 0 ) used = send_data ( id_rough_mom,  rough_mom,          Time )
if ( id_u_star     > 0 ) used = send_data ( id_u_star,     u_star,             Time )
if ( id_b_star     > 0 ) used = send_data ( id_b_star,     b_star,             Time )
if ( id_t_atm      > 0 ) used = send_data ( id_t_atm,      Atm%t_bot,          Time )
if ( id_u_atm      > 0 ) used = send_data ( id_u_atm,      Atm%u_bot,          Time )
if ( id_v_atm      > 0 ) used = send_data ( id_v_atm,      Atm%v_bot,          Time )
if ( id_albedo     > 0 ) used = send_data ( id_albedo,     albedo,             Time )
if ( id_u_flux     > 0 ) used = send_data ( id_u_flux,     flux_u,             Time )
if ( id_v_flux     > 0 ) used = send_data ( id_v_flux,     flux_v,             Time )

if ( id_eddy_t > 0) then
    work = spread(flux_t_zm,1,size(flux_t,1))
    call zonal_mean(flux_t,work_zm)
    work_2 = spread(work_zm,1,size(flux_t,1))
    work = work_2 - work
    used = send_data(id_eddy_t, work, Time)
endif
if ( id_eddy_q > 0) then
    work = spread(flux_q_zm,1,size(flux_q,1))
    call zonal_mean(flux_q,work_zm)
    work_2 = spread(work_zm,1,size(flux_t,1))
    work = work_2 - work
    used = send_data(id_eddy_q, work, Time)
endif
if ( id_eddy_u > 0) then
    work = spread(flux_u_zm,1,size(flux_u,1))
    call zonal_mean(flux_u,work_zm)
    work_2 = spread(work_zm,1,size(flux_t,1))
    work = work_2 - work
    used = send_data(id_eddy_u, work, Time)
endif
if ( id_eddy_v > 0) then
    work = spread(flux_v_zm,1,size(flux_v,1))
    call zonal_mean(flux_v,work_zm)
    work_2 = spread(work_zm,1,size(flux_t,1))
    work = work_2 - work
    used = send_data(id_eddy_v, work, Time)
endif


!=======================================================================

 end subroutine compute_flux

!#######################################################################

subroutine update_simple_surface (dt, Time, Atm, dt_t_atm, dt_q_atm  )

real, intent(in) :: dt
type       (time_type), intent(in)  :: Time
type (atmos_data_type), intent(in)  :: Atm
 
real, dimension(:,:),   intent(out) :: dt_t_atm, dt_q_atm
   
real, dimension(size(Atm%t_bot,1), size(Atm%t_bot,2)) :: &
                     gamma, dtmass, delta_t, delta_q, dflux_t, dflux_q, &
                     flux, deriv, dt_t_surf, &
                     entrop_evap, entrop_shflx, entrop_lwflx,sst_clamp

integer, dimension(2)                          :: bounds
real, dimension(is:ie,js:je)                   :: sin_lat, cos_lat, lat, f, cos_deriv, &
                                                  flux_u_zm, t_surf_zm, sst_zm, delta_lat, &
                                                  delta_t_surf, delta_flux
real, dimension(js:je)                         :: sin_lat_1d, cos_lat_1d, lat_1d
real    :: zt, dw, cp_inv, pi, coeff
integer :: j, half, i
logical :: used
character(len=64) :: ch1

allocate( f_global (is:ie,1:lat_max), &
          sin_lat_global (is:ie,1:lat_max), &
          cos_lat_global (is:ie,1:lat_max), &
          delta_t_surf_zm_global (is:ie,1:lat_max), &
          delta_lat_global (is:ie,1:lat_max), &
          lat_global (is:ie,1:lat_max), &
          delta_flux_global (is:ie,1:lat_max), &
          cos_deriv_global (is:ie,1:lat_max), &
          tau_x_zm_global (is:ie,1:lat_max), &
          t_surf_zm_global (is:ie,1:lat_max), & 
          enthalpy_global (is:ie,1:lat_max), &
          flux_o_global (is:ie,1:lat_max), & 
          flux_u_zm_global (is:ie,1:lat_max), & 
          sst_zm_global (is:ie,1:lat_max), & 
          flux_o_global_vector (1:lat_max), & 
          p (1:lat_max,1:lat_max), & 
          flux_meridional_global (is:ie,1:lat_max) )

half = lat_max/2
pi = 4.0*atan(1.)

call get_sin_lat(sin_lat_1d)
call get_cos_lat(cos_lat_1d)
call get_deg_lat(lat_1d)
lat_1d=lat_1d*pi/180

do i=is,ie
  lat(i,js:je) = lat_1d(js:je)
  cos_lat(i,js:je) = cos_lat_1d(js:je)
  sin_lat(i,js:je) = sin_lat_1d(js:je)
enddo
f(:,js:je) = (2*omega*sin_lat(:,js:je))
dw = 7*pi/180
cos_deriv(:,js:je) = 1 / (cos_lat(:,js:je) * radius)

call mpp_global_field(grid_domain, f, f_global)
call mpp_global_field(grid_domain, sin_lat, sin_lat_global)
call mpp_global_field(grid_domain, cos_lat, cos_lat_global)
call mpp_global_field(grid_domain, lat, lat_global)
call mpp_global_field(grid_domain, cos_deriv, cos_deriv_global)

flux_lw = Atm%flux_lw - flux_lw
   
dtmass  = Atm%Surf_Diff%dtmass
delta_t = Atm%Surf_Diff%delta_t
delta_q = Atm%Surf_Diff%delta_q
dflux_t = Atm%Surf_Diff%dflux_t
dflux_q = Atm%Surf_Diff%dflux_q

cp_inv = 1.0/cp_air

! temperature

gamma      =  1./ (1.0 - dtmass*(dflux_t + dhdt_atm*cp_inv))
e_t_n      =  dtmass*dhdt_surf*cp_inv*gamma
f_t_delt_n = (delta_t + dtmass * flux_t*cp_inv) * gamma    

flux_t     =  flux_t        + dhdt_atm * f_t_delt_n

dhdt_surf  =  dhdt_surf     + dhdt_atm * e_t_n   

! moisture

gamma      =  1./ (1.0 - dtmass*(dflux_q + dedq_atm))
e_q_n      =  dtmass*dedt_surf*gamma
f_q_delt_n = (delta_q  + dtmass * flux_q) * gamma    

flux_q     =  flux_q        + dedq_atm * f_q_delt_n

dedt_surf  =  dedt_surf     + dedq_atm * e_q_n

!################################################
!               OCEAN FLUX
!################################################
!Based on Levine and Schneider (2011)

call zonal_mean(flux_u,flux_u_zm(1,js:je))
call zonal_mean(sst,sst_zm(1,js:je))

do i=is,ie
  sst_zm(i,js:je) = sst_zm(1,js:je)
  flux_u_zm(i,js:je) = flux_u_zm(1,js:je)
enddo

call mpp_global_field(grid_domain,flux_u_zm,flux_u_zm_global)
call mpp_global_field(grid_domain,sst_zm,sst_zm_global)

do i=is,ie
  call centered_difference(sst_zm_global(i,:),delta_t_surf_zm_global(i,:),1,lat_max)
  call centered_difference(lat_global(i,:),delta_lat_global(i,:),1,lat_max)
enddo

!Find tropical belt latitudes
call easterly_bounds(flux_u_zm_global(1,:),bounds(:),lat_global(1,:))

delta_flux_global(:,:) = 0
flux_meridional_global(:,:) = 0
flux_o_global(:,:) = 0

!Calculate adiabatic ocean circulation
enthalpy_global = cp_ocean * flux_u_zm_global * delta_t_surf_zm_global / f_global

do j=bounds(1),half
  flux_meridional_global(:,j) = (cos_lat_global(:,j)**2) * sum(enthalpy_global(:,bounds(1):j),2)
enddo
do j=half+1,bounds(2)
  flux_meridional_global(:,j) = -(cos_lat_global(:,j)**2) * sum(enthalpy_global(:,j:bounds(2)),2)
enddo

!Calculate flux divergence
do i=is,ie
  call centered_difference(flux_meridional_global(i,:),delta_flux_global(i,:),bounds(1),bounds(2))
  do j=bounds(1),bounds(2)
    delta_flux_global(i,j) = -delta_flux_global(i,j) * cos_deriv_global(i,j) / delta_lat_global(i,j)
  enddo
enddo

flux_o_global_vector = 0
flux_o_global = 0
coeff = (1/(sqrt(2*pi)*dw))

!Gaussian smoother 
do j=1,lat_max
  do i=1,lat_max
    p(j,i) = exp(-((lat_global(is,j)-lat_global(is,i))**2)/(2*(dw**2)))*delta_lat_global(is,i)
  enddo
  zt = sum(p(j,:))
  p(j,:) = p(j,:) / zt
enddo

do j=2,lat_max-1
  flux_o_global_vector(j) = sum(delta_flux_global(is,2:lat_max-1)*p(j,2:lat_max-1))
enddo
do i=is,ie
  flux_o_global(i,1:lat_max) = flux_o_global_vector
enddo

call mpp_sync()

if(surface_choice == 1) then
   
  flux    = (flux_lw + Atm%flux_sw - hlf*Atm%fprec &
          - (flux_t + hlv*flux_q) + flux_o_global(is:ie,js:je))*dt/heat_capacity

  deriv   = - (dhdt_surf + hlv*dedt_surf + drdt_surf)*dt/heat_capacity 

  dt_t_surf = flux/(1.0 -deriv)

  if (do_sst_clamp) then 

    call set_domain(grid_domain)
    call read_data(sst_file, 'ts_field',   sst_clamp(:,:), grid_domain)
    sst = sst_clamp
  else
    sst = sst + dt_t_surf
    
  endif
  
elseif(surface_choice == 2) then

  dt_t_surf  = 0.0
  
endif

flux_t     = flux_t      + dt_t_surf*dhdt_surf
flux_q     = flux_q      + dt_t_surf*dedt_surf
flux_lw    = flux_lw     - dt_t_surf*drdt_surf
dt_t_atm   = f_t_delt_n  + dt_t_surf*e_t_n
dt_q_atm   = f_q_delt_n  + dt_t_surf*e_q_n

!=======================================================================
!-------------------- diagnostics section ------------------------------

if ( id_t_surf > 0 ) used = send_data ( id_t_surf, sst, Time )
if ( id_t_flux > 0 ) used = send_data ( id_t_flux, flux_t, Time )
if ( id_r_flux > 0 ) used = send_data ( id_r_flux, flux_lw, Time )
if ( id_q_flux > 0 ) used = send_data ( id_q_flux, flux_q, Time )
if ( id_o_flux > 0 ) used = send_data ( id_o_flux, flux_o_global(:,js:je), Time )
if ( id_o_flux_mer > 0 ) used = send_data ( id_o_flux_mer, flux_meridional_global(:,js:je), Time )
if ( id_entrop_evap > 0 ) then 
  entrop_evap = flux_q/sst
  used = send_data ( id_entrop_evap, entrop_evap, Time)
endif
if (id_entrop_shflx > 0 ) then
  entrop_shflx = flux_t/sst
  used = send_data ( id_entrop_shflx, entrop_shflx, Time)
endif
if (id_entrop_lwflx > 0 ) then
  entrop_lwflx = flux_lw/sst
  used = send_data ( id_entrop_lwflx, entrop_lwflx, Time)
endif

call sum_diag_integral_field ('evap', flux_q*86400.)

!=======================================================================
!---- deallocate module storage ----

deallocate (f_global, sin_lat_global, cos_lat_global, delta_t_surf_zm_global, &
           delta_lat_global, lat_global, delta_flux_global, cos_deriv_global, &
           tau_x_zm_global, enthalpy_global, flux_o_global, flux_u_zm_global, &
           sst_zm_global, flux_o_global_vector, p, flux_meridional_global)
deallocate (f_t_delt_n, f_q_delt_n, e_t_n, e_q_n)
deallocate(dhdt_surf, dedt_surf, dedq_surf, drdt_surf, dhdt_atm, dedq_atm, &
          flux_t, flux_q, flux_lw)  


!-----------------------------------------------------------------------

end subroutine update_simple_surface

!#######################################################################

subroutine simple_surface_init ( Time, Atm)

type       (time_type), intent(in)  :: Time
type (atmos_data_type), intent(in)  :: Atm

integer :: unit, ierr, io
integer :: j
integer :: layout(2)
real :: xx, xx2, lat

call get_grid_domain(is, ie, js, je)
call get_lat_max(lat_max)
!-----------------------------------------------------------------------
!------ read namelist ------

   if ( file_exist('input.nml')) then
      unit = open_namelist_file ( )
      ierr=1; do while (ierr /= 0)
         read  (unit, nml=simple_surface_nml, iostat=io, end=10)
         ierr = check_nml_error(io,'simple_surface_nml')
      enddo
 10   call close_file (unit)
   endif

!--------- write version number and namelist ------------------

 call write_version_number ( version, tagname )
 if ( mpp_pe() == mpp_root_pe() ) then
      write (stdlog(), nml=simple_surface_nml)
 endif

 call diag_integral_field_init ('evap', 'f6.3')
 call diag_field_init ( Time, Atm%axes(1:2) )
   
allocate(sst   (size(Atm%t_bot,1),size(Atm%t_bot,2)))
allocate(flux_u(size(Atm%t_bot,1),size(Atm%t_bot,2)))
allocate(flux_v(size(Atm%t_bot,1),size(Atm%t_bot,2)))
   
if(file_exist('INPUT/simple_surface.res.nc')) then
  call read_data('INPUT/simple_surface.res.nc', 'sst',    sst,    domain=Atm%domain)
  call read_data('INPUT/simple_surface.res.nc', 'flux_u', flux_u, domain=Atm%domain)
  call read_data('INPUT/simple_surface.res.nc', 'flux_v', flux_v, domain=Atm%domain)
else if(file_exist('INPUT/simple_surface.res')) then
  unit = open_file(file='INPUT/simple_surface.res',form='unformatted',&
                   action='read')
  call set_domain(Atm%domain)
  
  call read_data(unit, sst)
  call read_data(unit, flux_u)
  call read_data(unit, flux_v)
  
  call close_file(unit)
else 
  do j = 1, size(Atm%t_bot,2)
    lat = 0.5*(Atm%lat_bnd(j+1) + Atm%lat_bnd(j))
    xx = sin(lat)*sin(lat)
    xx2 = xx*xx
    sst(:,j) = 305.0 - 10.0*xx2 - 30.0*xx
  end do
  flux_u = 0.0
  flux_v = 0.0
endif
   

   do_init = .false.
!-----------------------------------------------------------------------

end subroutine simple_surface_init

!#######################################################################

subroutine diag_field_init ( Time, atmos_axes )

type(time_type), intent(in) :: Time
integer,         intent(in) :: atmos_axes(2)

integer :: iref
character(len=6) :: label_zm, label_zh
real, dimension(2) :: trange = (/  100., 400. /), &
                      vrange = (/ -400., 400. /), &
                      frange = (/ -0.01, 1.01 /)
!-----------------------------------------------------------------------
!  initializes diagnostic fields that may be output from this module
!  (the id numbers may be referenced anywhere in this module)
!-----------------------------------------------------------------------

!------ labels for diagnostics -------
!  (z_ref_mom, z_ref_heat are namelist variables)

iref = int(z_ref_mom+0.5)
if ( real(iref) == z_ref_mom ) then
                 write (label_zm,105) iref
  if (iref < 10) write (label_zm,100) iref
else
    write (label_zm,110) z_ref_mom
endif

iref = int(z_ref_heat+0.5)
if ( real(iref) == z_ref_heat ) then
                write (label_zh,105) iref
 if (iref < 10) write (label_zh,100) iref
else
    write (label_zh,110) z_ref_heat
endif

100 format (i1,' m',3x)
105 format (i2,' m',2x)
110 format (f4.1,' m')


id_wind = &
register_diag_field ( mod_name, 'wind', atmos_axes, Time, &
                    'wind speed for flux calculations', 'm/s', &
                     range=(/0.,vrange(2)/) )

id_drag_moist = &
register_diag_field ( mod_name, 'drag_moist', atmos_axes, Time, &
                    'drag coeff for moisture',    'none'     )

id_drag_heat  = &
register_diag_field ( mod_name, 'drag_heat', atmos_axes, Time, &
                    'drag coeff for heat',    'none'     )

id_drag_mom   = &
register_diag_field ( mod_name, 'drag_mom',  atmos_axes, Time, &
                    'drag coeff for momentum',     'none'     )

id_rough_moist = &
register_diag_field ( mod_name, 'rough_moist', atmos_axes, Time, &
                    'surface roughness for moisture',  'm'  )

id_rough_heat = &
register_diag_field ( mod_name, 'rough_heat', atmos_axes, Time, &
                    'surface roughness for heat',  'm'  )

id_rough_mom  = &
register_diag_field ( mod_name, 'rough_mom',  atmos_axes, Time, &
                    'surface roughness for momentum',  'm'  )

id_u_star     = &
register_diag_field ( mod_name, 'u_star',     atmos_axes, Time, &
                    'friction velocity',   'm/s'   )
  
id_b_star     = &
register_diag_field ( mod_name, 'b_star',     atmos_axes, Time, &
                    'buoyancy scale',      'm/s2'   )

id_u_flux     = &
register_diag_field ( mod_name, 'tau_x',      atmos_axes, Time, &
                    'zonal wind stress',     'pa'   )

id_v_flux     = &
register_diag_field ( mod_name, 'tau_y',      atmos_axes, Time, &
                    'meridional wind stress',     'pa'   )

id_eddy_u     = &
register_diag_field ( mod_name, 'eddy_u_stress',      atmos_axes, Time, &
                    'eddy zonal wind stress',     'm/s'   )

id_eddy_v     = &
register_diag_field ( mod_name, 'eddy_v_stress',      atmos_axes, Time, &
                    'eddy meridional wind stress',     'm/s'   )

id_eddy_t     = &
register_diag_field ( mod_name, 'eddy_t',      atmos_axes, Time, &
                    'eddy shflx',     'W/m**2/s'   )

id_eddy_q     = &
register_diag_field ( mod_name, 'eddy_q',      atmos_axes, Time, &
                    'eddy evap',     'kg/kg/s'   )

id_t_surf     = &
register_diag_field ( mod_name, 't_surf',     atmos_axes, Time, &
                    'surface temperature',    'deg_k', &
                    range=trange    )

id_t_flux     = &
register_diag_field ( mod_name, 'shflx',      atmos_axes, Time, &
                    'sensible heat flux',     'w/m2'    )

id_q_flux     = &
register_diag_field ( mod_name, 'evap',       atmos_axes, Time, &
                    'evaporation rate',        'kg/m2/s'  )

id_o_flux     = &
register_diag_field (mod_name, 'oflx',        atmos_axes, Time, &
                    'ocean heat divergence', 'w/m2' )

id_o_flux_mer    = &
register_diag_field (mod_name, 'oflx_meridional', atmos_axes, Time, &
                    'ocean heat flux', 'w' )

id_r_flux     = &
register_diag_field ( mod_name, 'lwflx',      atmos_axes, Time, &
                    'net (down-up) longwave flux',   'w/m2'    )

id_t_atm      = &
register_diag_field ( mod_name, 't_atm',      atmos_axes, Time, &
                    'temperature at btm level',    'deg_k', &
                    range=trange     )

id_u_atm      = &
register_diag_field ( mod_name, 'u_atm',      atmos_axes, Time, &
                    'u wind component at btm level',  'm/s', &
                    range=vrange    )

id_v_atm      = &
register_diag_field ( mod_name, 'v_atm',      atmos_axes, Time, &
                    'v wind component at btm level',  'm/s', &
                    range=vrange    )

id_t_ref      = &
register_diag_field ( mod_name, 't_ref',      atmos_axes, Time, &
                    'temperature at '//label_zh, 'deg_k' , &
                    range=trange      )

id_rh_ref     = &
register_diag_field ( mod_name, 'rh_ref',     atmos_axes, Time,   &
                    'relative humidity at '//label_zh, 'percent' )

id_u_ref      = &
register_diag_field ( mod_name, 'u_ref',      atmos_axes, Time, &
                    'zonal wind component at '//label_zm,  'm/s', &
                    range=vrange )

id_v_ref      = &
register_diag_field ( mod_name, 'v_ref',      atmos_axes, Time,     &
                  'meridional wind component at '//label_zm, 'm/s', &
                    range=vrange )

id_del_h      = &
register_diag_field ( mod_name, 'del_h',      atmos_axes, Time,  &
                    'ref height interp factor for heat', 'none' )
id_del_m      = &
register_diag_field ( mod_name, 'del_m',      atmos_axes, Time,     &
                    'ref height interp factor for momentum','none' )
id_del_q      = &
register_diag_field ( mod_name, 'del_q',      atmos_axes, Time,     &
                    'ref height interp factor for moisture','none' )
id_albedo      = &
register_diag_field ( mod_name, 'albedo',      atmos_axes, Time,     &
                    'surface albedo','none' )
id_entrop_evap      = &
register_diag_field ( mod_name, 'entrop_evap', atmos_axes, Time,     &
                    'entropy source from evap','kg/m2/s/K' )

id_entrop_shflx      = &
register_diag_field ( mod_name, 'entrop_shflx', atmos_axes, Time,     &
                    'entropy source from SH flux','w/m2/K' )

id_entrop_lwflx      = &
register_diag_field ( mod_name, 'entrop_lwflx', atmos_axes, Time,     &
                    'entropy source from LW flux','w/m2/K' )

!-----------------------------------------------------------------------

end subroutine diag_field_init


!########################################################################

subroutine simple_surface_end (Atm)

type (atmos_data_type), intent(in)  :: Atm
integer :: unit

call write_data('RESTART/simple_surface.res.nc', 'sst',    sst,    domain=Atm%domain)
call write_data('RESTART/simple_surface.res.nc', 'flux_u', flux_u, domain=Atm%domain)
call write_data('RESTART/simple_surface.res.nc', 'flux_v', flux_v, domain=Atm%domain)

end subroutine simple_surface_end

!#######################################################################

subroutine zonal_mean(input_data,zonal_data)
real, dimension(:,:), intent(in)    :: input_data
real, dimension(:), intent(out)     :: zonal_data
real                                :: length_data
integer                             :: j

length_data = size(input_data,1)
do j=1,size(input_data,2)
    zonal_data(j) = sum(input_data(:,j),1)/length_data
enddo

end subroutine zonal_mean

!#######################################################################
subroutine easterly_bounds(tau_x,bounds,lat)

!finds the first latitude poleward of the equator where the zonal surface 
!stress changes sign from easterly to westerly
!includes a catch for very high latitude edges, should be relaxed if studying substan

real, intent(in), dimension(1:lat_max)     :: tau_x
real, intent(in), dimension(1:lat_max)     :: lat
integer, intent(out), dimension(2) :: bounds
integer :: j, half, loc
integer :: sign_change
real    :: current, lat_abs

half = lat_max/2

sign_change = 1
j = half - 2
do while (sign_change .eq. 1) 
  if (tau_x(j) .lt. 0) then
    sign_change = 0
    bounds(1) = j+1
  endif
  j = j - 1
  if (lat(j) .le. -60) then
    sign_change = 0
    bounds(1) = j+1
  endif 
enddo 
sign_change = 1
j = half + 3
do while (sign_change .eq. 1) 
  if (tau_x(j) .lt. 0) then
    sign_change = 0
    bounds(2) = j-1
  endif
  j = j + 1
  if (lat(j) .ge. 60) then
    sign_change = 0
    bounds(2) = j-1
  endif 
enddo 

end subroutine easterly_bounds
!#######################################################################
subroutine centered_difference(data,derivative,forward,backward)

real, intent(in), dimension(1:lat_max)    :: data
real, intent(out), dimension(1:lat_max)   :: derivative
integer :: j
integer, intent(in) :: forward, backward

derivative(:) = 0

do j=forward+1,backward-1
  derivative(j) = .5*(data(j+1) - data(j-1))
enddo
derivative(forward) = data(forward+1) - data(forward)
derivative(backward) = data(backward) - data(backward-1)

end subroutine centered_difference
!#######################################################################
end module simple_surface_mod

