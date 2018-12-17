module grey_radiation_mod

! ==================================================================================
! ==================================================================================

!  use utilities_mod,         only: error_mesg, open_file, file_exist, &
!                                   check_nml_error, FATAL, get_my_pe, & 
!                                   close_file, set_domain, get_num_pes, &
!                                   read_data, write_data, &
!                                   check_system_clock, NOTE, &
!                                   get_domain_decomp, check_system_clock

   use             fms_mod,   only: open_namelist_file, check_nml_error,  &
                                    mpp_pe, mpp_root_pe, close_file, &
                                    write_version_number, stdlog

   use       constants_mod,   only: stefan, cp_air, grav

   use    diag_manager_mod,   only: register_diag_field, send_data

   use    time_manager_mod,   only: time_type, set_date, set_time,  &
                                    get_time,    operator(+),       &
                                    operator(-), operator(/=), get_date



!==================================================================================
implicit none
private
!==================================================================================

! version information 

character(len=128), parameter :: version = &
'$Id: grey_radiation.f90,v 1.1.2.2 2005/05/21 02:01:56 pjp Exp $'

character(len=128), parameter :: tagname = &
'$Name:  $'

!==================================================================================

! public interfaces

public :: grey_radiation_init, grey_radiation, grey_radiation_end
!==================================================================================


! module variables
logical :: initialized =.false.
real, parameter :: p00 = 1000.e2

real    :: solar_constant  = 1360.0
real    :: del_sol         = 0.0
! modif omp: winter/summer hemisphere
real    :: del_sw          = 0.0
real    :: ir_tau_eq       = 4.0
real    :: ir_tau_pole     = 4.0
real    :: atm_abs         = 0.2
real    :: frac_oz         = 0.0
real    :: sw_diff         = 0.0
real    :: long_pert       = 180.
real    :: del_long        =  30.
real    :: size_pert       =  0.
real    :: linear_tau      = 0.1

real    :: lat_pert        = 0.0
real    :: lon_pert        = 180.0
real    :: del_lat         = 30.0
real    :: del_lon         = 90.0
real    :: fcng_pert       = 0.0

real    :: wave_amp        = 0.0
real    :: wave_lon        = 180.0
real    :: wave_lat        = 0.0
real    :: wave_del_lon    = 30.0
real    :: wave_del_lat    = 20.0
real    :: wave_period     = 20.0
real    :: wave_env        = 80.0
logical :: wave_source     = .FALSE.

real    :: forcing_height_center  = 0.0
real    :: forcing_depth          = 0.0
real    :: forcing_lat_width      = 0.0
real    :: forcing_lon_width      = 0.0
real    :: forcing_lat_center     = 0.0
real    :: forcing_magnitude      = 0.0 
logical :: forcing_asymmetry       = .FALSE.
real    :: forcing_wavenum        = 0.0
real    :: forcing_horiz_sharp    = 4.0
real    :: forcing_vert_sharp     = 10.0




real, save :: pi, deg_to_rad , rad_to_deg

namelist/grey_radiation_nml/ solar_constant, del_sol, &
           ir_tau_eq, ir_tau_pole, atm_abs, sw_diff, long_pert, del_long, &
           size_pert, linear_tau, del_sw,                    &
           lat_pert, lon_pert, del_lat, del_lon, fcng_pert, &
           wave_amp, wave_lon, wave_lat, wave_del_lon,      &
           wave_del_lat, wave_period, wave_env, wave_source, frac_oz
		   
		   !forcing_depth, forcing_lat_width, forcing_asymmetry, &
		   !forcing_lon_width, forcing_magnitude, forcing_height_center &
		   !forcing_lat_center, forcing_wavenum, forcing_horiz_sharp, &
		   !forcing_vert_sharp, ozone_abs

!==================================================================================
!-------------------- diagnostics fields -------------------------------

integer :: id_olr, id_swdn_sfc, id_swdn_toa, id_lwdn_sfc, id_lwup_sfc, &
           id_tdt_rad, id_flux_rad, id_flux_lw, id_flux_sw, id_entrop_rad, &
           id_tdt_sw, id_tdt_lw

character(len=14), parameter :: mod_name = 'grey_radiation'

real :: missing_value = -999.


contains



! ==================================================================================
! ==================================================================================


subroutine grey_radiation_init(axes, Time)

!-------------------------------------------------------------------------------------
integer, intent(in), dimension(4) :: axes
type(time_type), intent(in)       :: Time
!-------------------------------------------------------------------------------------
integer, dimension(3) :: half = (/1,2,4/)
integer :: ierr, io, unit
!-----------------------------------------------------------------------------------------
! read namelist and copy to logfile

unit = open_namelist_file ( )
ierr=1
do while (ierr /= 0)
   read  (unit, nml=grey_radiation_nml, iostat=io, end=10)
   ierr = check_nml_error (io, 'grey_radiation_nml')
enddo
10 call close_file (unit)

call write_version_number ( version, tagname )
if ( mpp_pe() == mpp_root_pe() ) then
  write (stdlog(), nml=grey_radiation_nml)
endif
call close_file (unit)

pi    = 4.0*atan(1.)
deg_to_rad = 2.*pi/360.
rad_to_deg = 360.0/2./pi

initialized = .true.

!-----------------------------------------------------------------------
!------- initialize quantities for integral package -------

!       call diag_integral_field_init ('olr',    'f8.3')
!       call diag_integral_field_init ('abs_sw', 'f8.3')

!-----------------------------------------------------------------------
!------------ initialize diagnostic fields ---------------

    id_olr = &
    register_diag_field ( mod_name, 'olr', axes(1:2), Time, &
               'outgoing longwave radiation', &
               'watts/m2', missing_value=missing_value               )

    id_swdn_sfc = &
    register_diag_field ( mod_name, 'swdn_sfc', axes(1:2), Time, &
               'SW flux down at surface', &
               'watts/m2', missing_value=missing_value               )
    id_swdn_toa = &
    register_diag_field ( mod_name, 'swdn_toa', axes(1:2), Time, &
               'SW flux down at TOA', &
               'watts/m2', missing_value=missing_value               )
    id_lwup_sfc = &
    register_diag_field ( mod_name, 'lwup_sfc', axes(1:2), Time, &
               'LW flux up at surface', &
               'watts/m2', missing_value=missing_value               )

    id_lwdn_sfc = &
    register_diag_field ( mod_name, 'lwdn_sfc', axes(1:2), Time, &
               'LW flux down at surface', &
               'watts/m2', missing_value=missing_value               )

    id_tdt_rad = &
        register_diag_field ( mod_name, 'tdt_rad', axes(1:3), Time, &
               'Temperature tendency due to radiation', &
               'K/s', missing_value=missing_value               )

    id_tdt_lw = &
        register_diag_field ( mod_name, 'tdt_lw', axes(1:3), Time, &
                'Temperature tendency due to LW radiation', &
                'K/s', missing_value=missing_value               )
    id_tdt_sw = &
        register_diag_field ( mod_name, 'tdt_sw', axes(1:3), Time, &
                'Temperature tendency due to SW radiation', &
                'K/s', missing_value=missing_value               )

    id_flux_rad = &
        register_diag_field ( mod_name, 'flux_rad', axes(half), Time, &
               'Total radiative flux (positive up)', &
               'W/m^2', missing_value=missing_value               )
    id_flux_lw = &
        register_diag_field ( mod_name, 'flux_lw', axes(half), Time, &
               'Net longwave radiative flux (positive up)', &
               'W/m^2', missing_value=missing_value               )
    id_flux_sw = &
        register_diag_field ( mod_name, 'flux_sw', axes(half), Time, &
               'Net shortwave radiative flux (positive up)', &
               'W/m^2', missing_value=missing_value               )
    id_entrop_rad = &
            register_diag_field ( mod_name, 'entrop_rad', axes(1:3), Time, &
               'Entropy production by radiation', &
               '1/s', missing_value=missing_value               )





return
end subroutine grey_radiation_init


! ==================================================================================

subroutine grey_radiation (is, js, Time_diag, lat, lon, p_half, albedo, t_surf, t, tdt, net_surf_sw_down, surf_lw_down)

integer, intent(in)                 :: is, js
type(time_type), intent(in)         :: Time_diag
real, intent(in) , dimension(:,:)   :: lat, lon, albedo
real, intent(in) , dimension(:,:)   :: t_surf
real, intent(in) , dimension(:,:,:) :: t, p_half
real, intent(inout), dimension(:,:,:) :: tdt
real, intent(out), dimension(:,:)   :: net_surf_sw_down, surf_lw_down

real, dimension(size(t,2)) :: ss, ss2, ss4, ss6, ss8, solar, tau_0, solar_tau_0, p2, p_center, p_width
real, dimension(size(t,1), size(t,2))              :: b_surf 
real, dimension(size(t,1), size(t,2), size(t,3))   :: b, tdt_rad, entrop_rad, tdt_rad_lw, tdt_rad_sw
real, dimension(size(t,1), size(t,2), size(t,3)+1) :: up, down, net, solar_down, flux_rad, flux_sw 
!forcing_perturbation
real, dimension(size(t,2), size(t,3)  )   :: dtrans
real, dimension(size(t,2), size(t,3)+1)   :: tau, solar_tau
real, dimension(size(t,1), size(t,2))   :: olr, swin

!real, dimension(size(t,1), size(t,2))              :: walker_forcing, wave_forcing

integer :: i, j, k, n

real :: surf_albedo
real :: dist
logical :: used

n = size(t,3)
!forcing_perturbation = 0.

! forcing_depth, forcing_lat_width, forcing_symmetry, &
!		   forcing_lon_width, forcing_magnitude, forcing_height_center &
!		   forcing_lat_center, forcing_wavenum, forcing_horiz_sharp, &
		   !forcing_vert_sharp
		   

		   
		   
!forcing perturbation:

!do i=1,size(t,1)
 !  do j=1,size(t,2)
 !     do k=1,size(t,3)
 !        forcing_perturbation(i,j,k)=forcing_magnitude*exp(-(rad_to_deg* &
!	     lat(i,j)-forcing_lat_center)**forcing_horiz_sharp/(forcing_lat_width** &
!	     forcing_horiz_sharp))*exp(-(p_half(i,j,k)-forcing_height_center)**forcing_vert_sharp &
!	     /(forcing_depth**forcing_vert_sharp))/86400
!	 end do
!   end do
!end do
   
!if (forcing_asymmetry) then
!   do i=1,size(t,1)
!      do j=1,size(t,2)
!		forcing_perturbation(i,j,:)=forcing_perturbation(i,j,:)*(cosd(forcing_wavenum*lon(i,j)/2)^4)
!	  end do
!	end do
!end if 


!do i=1,size(t,1)
!   do j=1,size(t,2)
!      long_forcing(i,j) = size_pert*exp(-(rad_to_deg*lon(i,j)-long_pert)**2./ &
!                         del_long**2.)
!   end do
!end do

!do i =  1,size(t,1)
!  do j = 1,size(t,2)
!    dist = SQRT(((rad_to_deg*lat(i,j) - lat_pert)/del_lat)**2                   &
!         + ((rad_to_deg*lon(i,j) - lon_pert)/del_lon)**2)
!    if (dist .lt. 1) then
!        walker_forcing(i,j) = fcng_pert*cos(pi * dist*0.5)**2;
!       else
!        walker_forcing(i,j) = 0.0
!     end if
!  end do
!end do

ss  = sin(lat(1,:))
ss2 = ss*ss
p2 = (1. - 3.*ss*ss)/4.

solar = 0.25*solar_constant*(1.0 + del_sol*p2 + del_sw * ss)

tau_0 = ir_tau_eq +(ir_tau_pole - ir_tau_eq)*ss*ss

solar_tau_0 = (1.0 - sw_diff*ss*ss)*atm_abs

b = stefan*t*t*t*t
b_surf = stefan*t_surf*t_surf*t_surf*t_surf

p_center = (130- 90 * cos(0.9*lat(1,:))**2)*100
p_width = (20+40*sin(1*lat(1,:))**2)*100

do k = 1, n+1

!   tau(:,k)       = tau_0(:)      *(p_half(1,:,k)/p_half(1,:,n+1))**2

! modif omp: optical thickness is obtained by combination of linear and quadratic function. 
! linear_tau controls the fraction of the optical depth that varies linearily. The 
! surface pressure has been replaced by a reference value (p00 = 1000mb) from constant_mod.
 
!  tau(:,k)       = tau_0(:) * (linear_tau * p_half(1,:,k)/p00 + (1.0 - linear_tau) &
!       * (p_half(1,:,k)/p00)**2)
! modif df 1-23-04: changing profile of IR absorber
   tau(:,k)       = tau_0(:) * (linear_tau * p_half(1,:,k)/p00 + (1.0 - linear_tau) &
       * (p_half(1,:,k)/p00)**4)

  ! Split the solar optical depth equally into water vapor absorption and idealized ozone absorption
  ! Ozone absorption peaks in the lower stratosphere and is parameterized as the error function 
  ! Water vapor is as in Frierson et al. 2006 and Dargan Frierson's Ph.D. dissertation 
  do j=1,size(lat,2)
    solar_tau(j,k) = (1-frac_oz)*solar_tau_0(j)*(p_half(1,j,k)/p00)**4 &
        + solar_tau_0(j)*frac_oz*(.25*sqrt(pi)*(erf(((p_half(1,j,k)-140)/60))+1)+&
        (cos(0.3*lat(1,j))**2)*0.3333*sqrt(pi)*(erf((((p_half(1,j,k)-p_center(j))/p_width(j))))+1))
  end do
end do

do k = 1, n
  dtrans(:,k) = exp(-(tau(:,k+1)-tau(:,k)))
end do

up(:,:,n+1) = b_surf
do k = n,1,-1
  do j = 1, size(t,2)
    up(:,j,k) = up(:,j,k+1)*dtrans(j,k) + b(:,j,k)*(1.0 - dtrans(j,k))
  end do
end do

down(:,:,1) = 0.0
do k = 1,n
  do j =1, size(t,2)
    down(:,j,k+1) = down(:,j,k)*dtrans(j,k) + b(:,j,k)*(1.0 - dtrans(j,k))
  end do
end do

do i = 1, size(t,1)
do j = 1, size(t,2)
  do k = 1,n+1
!    solar_down(i,j,k) = (solar(j)+long_forcing(i,j))*exp(-solar_tau(j,k))
!         solar_down(i,j,k) = (walker_forcing(i,j)+wave_forc(i,j)  &
!              + solar(j))*exp(-solar_tau(j,k))
         solar_down(i,j,k) = (solar(j))*exp(-solar_tau(j,k))
  end do
end do
end do

do k = 1,n+1
  net(:,:,k) = up(:,:,k)-down(:,:,k)
  flux_sw(:,:,k) = albedo(:,:)*solar_down(:,:,n+1) - solar_down(:,:,k)
  flux_rad(:,:,k) = net(:,:,k) + flux_sw(:,:,k)
end do

do k = 1,n
  tdt_rad(:,:,k) = (net(:,:,k+1) - net(:,:,k) - solar_down(:,:,k+1) + solar_down(:,:,k))  &
             *grav/(cp_air*(p_half(:,:,k+1)-p_half(:,:,k)))
  tdt_rad_sw(:,:,k) = (- solar_down(:,:,k+1) + solar_down(:,:,k))*grav/(cp_air*(p_half(:,:,k+1)-        p_half(:,:,k)))
  tdt_rad_lw(:,:,k) = (net(:,:,k+1) - net(:,:,k))*grav/(cp_air*(p_half(:,:,k+1)-p_half(:,:,k)))
  tdt(:,:,k) = tdt(:,:,k) + tdt_rad(:,:,k)
end do



surf_lw_down     = down(:,:,n+1)
net_surf_sw_down = solar_down(:,:,n+1)*(1. - albedo(:,:))
olr = up(:,:,1)
swin = solar_down(:,:,1)

    

!------- outgoing lw flux toa (olr) -------
      if ( id_olr > 0 ) then
          used = send_data ( id_olr, olr, Time_diag, is, js )
      endif

!------- downward sw flux surface -------
      if ( id_swdn_sfc > 0 ) then
          used = send_data ( id_swdn_sfc, net_surf_sw_down, Time_diag, is, js )
      endif
!------- incoming sw flux toa -------
      if ( id_swdn_toa > 0 ) then
          used = send_data ( id_swdn_toa, swin, Time_diag, is, js )
      endif
!------- upward lw flux surface -------
      if ( id_lwup_sfc > 0 ) then
          used = send_data ( id_lwup_sfc, b_surf, Time_diag, is, js )
      endif

!------- downward lw flux surface -------
      if ( id_lwdn_sfc > 0 ) then
          used = send_data ( id_lwdn_sfc, surf_lw_down, Time_diag, is, js )
      endif
!------- temperature tendency due to radiation ------------
      if ( id_tdt_rad > 0 ) then
         used = send_data ( id_tdt_rad, tdt_rad, Time_diag, is, js, 1 )
      endif
    if ( id_tdt_sw > 0 ) then
        used = send_data ( id_tdt_sw, tdt_rad_sw, Time_diag, is, js, 1 )
    endif
    if ( id_tdt_lw > 0 ) then
        used = send_data ( id_tdt_lw, tdt_rad_lw, Time_diag, is, js, 1 )
    endif
!------- total radiative flux (at half levels) -----------
      if ( id_flux_rad > 0 ) then
         used = send_data ( id_flux_rad, flux_rad, Time_diag, is, js, 1 )
      endif
!------- longwave radiative flux (at half levels) --------
      if ( id_flux_lw > 0 ) then 
         used = send_data ( id_flux_lw, net, Time_diag, is, js, 1 )
      endif
      if ( id_flux_sw > 0 ) then
         used = send_data ( id_flux_sw, flux_sw, Time_diag, is, js, 1 )
      endif
      if ( id_entrop_rad > 0 ) then
         do k=1,n 
            entrop_rad(:,:,k) =tdt_rad(:,:,k)/t(:,:,k)*p_half(:,:,n+1)/1.e5
         end do
         used = send_data ( id_entrop_rad, entrop_rad, Time_diag, is, js, 1 )
      endif

return
end subroutine grey_radiation

! ==================================================================================

real function erf(x)

   implicit none
   
   real, intent(in) :: x
   
   erf=(2/sqrt(pi))*(x/abs(x))*sqrt(1-exp(-(x**2)))*((sqrt(pi)/2)+ (31/200)*exp(-(x**2)) - (341/8000)*(exp(-2*(x**2))))

end function erf
     
! ==================================================================================

subroutine grey_radiation_end()
end subroutine grey_radiation_end

! ==================================================================================

end module grey_radiation_mod
