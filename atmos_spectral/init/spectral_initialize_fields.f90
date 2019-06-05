module spectral_initialize_fields_mod

use              fms_mod, only: mpp_pe, mpp_root_pe, write_version_number, read_data

use        constants_mod, only: rdgas

use       transforms_mod, only: trans_grid_to_spherical, trans_spherical_to_grid, vor_div_from_uv_grid, &
                                uv_grid_from_vor_div, get_grid_domain, get_spec_domain, area_weighted_global_mean

implicit none
private

public :: spectral_initialize_fields

character(len=128), parameter :: version = &
'$Id: spectral_initialize_fields.f90,v 10.0 2003/10/24 22:00:59 fms Exp $'

character(len=128), parameter :: tagname = &
'$Name: lima $'

logical :: entry_to_logfile_done = .false.

contains

!-------------------------------------------------------------------------------------------------
subroutine spectral_initialize_fields(reference_sea_level_press, triang_trunc, choice_of_init, initial_temperature, &
                        surf_geopotential, ln_ps, vors, divs, ts, psg, ug, vg, tg, vorg, divg, zonal_perturbation)

real,    intent(in) :: reference_sea_level_press
logical, intent(in) :: triang_trunc
integer, intent(in) :: choice_of_init
real,    intent(in) :: initial_temperature

real,    intent(in),  dimension(:,:    ) :: surf_geopotential
complex, intent(out), dimension(:,:    ) :: ln_ps
complex, intent(out), dimension(:,:,:  ) :: vors, divs, ts
real,    intent(out), dimension(:,:    ) :: psg
real,    intent(out), dimension(:,:,:  ) :: ug, vg, tg
real,    intent(out), dimension(:,:,:  ) :: vorg, divg

real, allocatable, dimension(:,:) :: ln_psg

logical, intent(in) :: zonal_perturbation

real :: initial_sea_level_press, global_mean_psg
real :: initial_perturbation   = 1.e-7

integer :: ms, me, ns, ne, is, ie, js, je, num_levels, num_lat, num_lon

if(.not.entry_to_logfile_done) then
  call write_version_number(version, tagname)
  entry_to_logfile_done = .true.
endif

num_levels = size(ug,3)

call get_grid_domain(is, ie, js, je)
call get_spec_domain(ms, me, ns, ne)
allocate(ln_psg(is:ie, js:je))

initial_sea_level_press = reference_sea_level_press  

ug      = 0.
vg      = 0.
tg      = initial_temperature
psg     = 0.
vorg    = 0.
divg    = 0.

vors  = (0.,0.)
divs  = (0.,0.)
ts    = (0.,0.)
ln_ps = (0.,0.)

tg     = initial_temperature
ln_psg = log(initial_sea_level_press) - surf_geopotential/(rdgas*initial_temperature)
psg    = exp(ln_psg)

if (zonal_perturbation) then

  if(choice_of_init == 1) then  ! perturb temperature field
    if(is <= 1 .and. ie >= 1 .and. js <= 1 .and. je >= 1) then
      tg(1,1,:) = tg(1,1,:) + 1.0
    endif
  endif

  if(choice_of_init == 2) then   ! initial vorticity perturbation used in benchmark code
    if(ms <= 1 .and. me >= 1 .and. ns <= 3 .and. ne >= 3) then
      vors(2-ms,4-ns,num_levels  ) = initial_perturbation
      vors(2-ms,4-ns,num_levels-1) = initial_perturbation
      vors(2-ms,4-ns,num_levels-2) = initial_perturbation
    endif
    if(ms <= 5 .and. me >= 5 .and. ns <= 3 .and. ne >= 3) then
      vors(6-ms,4-ns,num_levels  ) = initial_perturbation
      vors(6-ms,4-ns,num_levels-1) = initial_perturbation
      vors(6-ms,4-ns,num_levels-2) = initial_perturbation
    endif
    if(ms <= 1 .and. me >= 1 .and. ns <= 2 .and. ne >= 2) then
      vors(2-ms,3-ns,num_levels  ) = initial_perturbation
      vors(2-ms,3-ns,num_levels-1) = initial_perturbation
      vors(2-ms,3-ns,num_levels-2) = initial_perturbation
    endif
    if(ms <= 5 .and. me >= 5 .and. ns <= 2 .and. ne >= 2) then
      vors(6-ms,3-ns,num_levels  ) = initial_perturbation
      vors(6-ms,3-ns,num_levels-1) = initial_perturbation
      vors(6-ms,3-ns,num_levels-2) = initial_perturbation
    endif

    if(is <= 1 .and. ie >= 1 .and. js <= 1 .and. je >= 1) then
      tg(1,1,:) = tg(1,1,:) + 1.0
    endif
    
    call uv_grid_from_vor_div(vors, divs, ug, vg)
  endif
else
   if (js<=1 .and. je>=1) then
       ug(:,js,:) = initial_perturbation
       vg(:,js,:) = initial_perturbation
   endif
endif

call trans_grid_to_spherical(tg, ts)
call trans_spherical_to_grid(ts, tg)

call trans_grid_to_spherical(ln_psg, ln_ps)
call trans_spherical_to_grid(ln_ps,  ln_psg)
psg = exp(ln_psg)

call vor_div_from_uv_grid(ug, vg, vors, divs, triang=triang_trunc)
call uv_grid_from_vor_div(vors, divs, ug, vg)
call trans_spherical_to_grid(vors, vorg)
call trans_spherical_to_grid(divs, divg)

!  compute and print mean surface pressure
global_mean_psg = area_weighted_global_mean(psg)
if(mpp_pe() == mpp_root_pe()) then
  print '("mean surface pressure=",f9.4," mb")',.01*global_mean_psg
endif

return
end subroutine spectral_initialize_fields
!================================================================================

end module spectral_initialize_fields_mod
