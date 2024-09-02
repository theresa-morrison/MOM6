module MOM_SIS_set_ocean_top_stress 

use MOM_domains,       only : pass_var, pass_vector, AGRID, BGRID_NE, CGRID_NE
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_open_boundary, only : OBC_NONE
use MOM_open_boundary, only : OBC_DIRECTION_E, OBC_DIRECTION_W, OBC_DIRECTION_N, OBC_DIRECTION_S
use MOM_time_manager,  only : operator(+), operator(-)
use MOM_time_manager,  only : operator(>), operator(*), operator(/), operator(/=)
use MOM_unit_scaling,  only : unit_scale_type

use MOM_SIS_hor_grid,  only : SIS_hor_grid_type
use SIS_open_boundary, only : OBC_segment_type
use SIS_open_boundary, only : ice_OBC_type, OBC_segment_type
use SIS_types,         only : ice_ocean_flux_type, fast_ice_avg_type
use ice_grid,          only : ice_grid_type

use MOM_SIS_dyn_types,     only : SIS_C_dyn_CS, SIS_dyn_state_2d, FIA_2d  

implicit none ; private

#include <SIS2_memory.h>

public :: set_ocean_top_stress_C2

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Calculate the stresses on the ocean integrated across all the thickness categories with the
!! appropriate staggering, and store them in the public ice data type for use by the ocean
!! model.  This version of the routine uses wind and ice-ocean stresses on a C-grid.
subroutine set_ocean_top_stress_C2(IOF, windstr_x_water, windstr_y_water, &
                                   str_ice_oce_x, str_ice_oce_y, ice_free, ice_cover, G, US, OBC)
  type(ice_ocean_flux_type), intent(inout) :: IOF !< A structure containing fluxes from the ice to
                                                  !! the ocean that are calculated by the ice model.
  type(SIS_hor_grid_type),   intent(inout) :: G   !< The horizontal grid type
  real, dimension(SZIB_(G),SZJ_(G)), &
                             intent(in)    :: windstr_x_water !< The x-direction wind stress over
                                                  !! open water [R Z L T-2 ~> Pa].
  real, dimension(SZI_(G),SZJB_(G)), &
                             intent(in)    :: windstr_y_water !< The y-direction wind stress over
                                                  !! open water [R Z L T-2 ~> Pa].
  real, dimension(SZIB_(G),SZJ_(G)), &
                             intent(in)    :: str_ice_oce_x !< The x-direction ice to ocean stress [R Z L T-2 ~> Pa]
  real, dimension(SZI_(G),SZJB_(G)), &
                             intent(in)    :: str_ice_oce_y !< The y-direction ice to ocean stress [R Z L T-2 ~> Pa]
  real, dimension(SZI_(G),SZJ_(G)), &
                             intent(in)    :: ice_free  !< The fractional open water area coverage [nondim], 0-1
  real, dimension(SZI_(G),SZJ_(G)), &
                             intent(in)    :: ice_cover !< The fractional ice area coverage [nondim], 0-1
  type(unit_scale_type),     intent(in)    :: US  !< A structure with unit conversion factors
  type(ice_OBC_type),        pointer       :: OBC  !< Open boundary structure.

  real    :: ps_ice, ps_ocn ! ice_free and ice_cover interpolated to a velocity point [nondim].
  integer :: i, j, k, isc, iec, jsc, jec
  integer :: l_seg
  logical :: local_open_u_BC, local_open_v_BC
  type(OBC_segment_type), pointer :: segment => NULL()

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  if (IOF%stress_count == 0) then
    IOF%flux_u_ocn(:,:) = 0.0 ; IOF%flux_v_ocn(:,:) = 0.0
  endif

  local_open_u_BC = .false. ; local_open_v_BC = .false.
  if (associated(OBC)) then ; if (OBC%OBC_pe) then
    local_open_u_BC = OBC%open_u_BCs_exist_globally
    local_open_v_BC = OBC%open_v_BCs_exist_globally
  endif ; endif

  if ((local_open_u_BC .or. local_open_v_BC) .and. &
      (IOF%flux_uv_stagger == AGRID) .or. (IOF%flux_uv_stagger == BGRID_NE)) &
        call SIS_error(FATAL, "No open boundaries for given flux staggering")

  !   Copy and interpolate the ice-ocean stress_Cgrid.  This code is slightly
  ! complicated because there are 3 different staggering options supported.

  if (IOF%flux_uv_stagger == AGRID) then
    !$OMP parallel do default(shared) private(ps_ocn, ps_ice)
    do j=jsc,jec ; do i=isc,iec
      ps_ocn = G%mask2dT(i,j) * ice_free(i,j)
      ps_ice = G%mask2dT(i,j) * ice_cover(i,j)
      IOF%flux_u_ocn(i,j) = IOF%flux_u_ocn(i,j) + &
           (ps_ocn * 0.5 * (windstr_x_water(I,j) + windstr_x_water(I-1,j)) + &
            ps_ice * 0.5 * (str_ice_oce_x(I,j) + str_ice_oce_x(I-1,j)) )
      IOF%flux_v_ocn(i,j) = IOF%flux_v_ocn(i,j) + &
           (ps_ocn * 0.5 * (windstr_y_water(i,J) + windstr_y_water(i,J-1)) + &
            ps_ice * 0.5 * (str_ice_oce_y(i,J) + str_ice_oce_y(i,J-1)) )
    enddo ; enddo
  elseif (IOF%flux_uv_stagger == BGRID_NE) then
    !$OMP parallel do default(shared) private(ps_ocn, ps_ice)
    do J=jsc-1,jec ; do I=isc-1,iec
      ps_ocn = 1.0 ; ps_ice = 0.0
      if (G%mask2dBu(I,J)>0.0) then
        ps_ocn = 0.25 * ((ice_free(i+1,j+1) + ice_free(i,j)) + &
                         (ice_free(i+1,j) + ice_free(i,j+1)) )
        ps_ice = 0.25 * ((ice_cover(i+1,j+1) + ice_cover(i,j)) + &
                         (ice_cover(i+1,j) + ice_cover(i,j+1)) )
      endif
      IOF%flux_u_ocn(I,J) = IOF%flux_u_ocn(I,J) + &
          (ps_ocn * 0.5 * (windstr_x_water(I,j) + windstr_x_water(I,j+1)) + &
           ps_ice * 0.5 * (str_ice_oce_x(I,j) + str_ice_oce_x(I,j+1)) )
      IOF%flux_v_ocn(I,J) = IOF%flux_v_ocn(I,J) + &
          (ps_ocn * 0.5 * (windstr_y_water(i,J) + windstr_y_water(i+1,J)) + &
           ps_ice * 0.5 * (str_ice_oce_y(i,J) + str_ice_oce_y(i+1,J)) )
    enddo ; enddo
  elseif (IOF%flux_uv_stagger == CGRID_NE) then
    if (local_open_u_BC) then
      !$OMP parallel do default(shared) private(ps_ocn, ps_ice)
      do j=jsc,jec ; do I=Isc-1,iec
        ps_ocn = 1.0 ; ps_ice = 0.0
        l_seg = OBC%segnum_u(I,j)
        if (G%mask2dCu(I,j)>0.0) then
          ps_ocn = 0.5*(ice_free(i+1,j) + ice_free(i,j))
          ps_ice = 0.5*(ice_cover(i+1,j) + ice_cover(i,j))
        endif
        if (l_seg /= OBC_NONE) then
          if (OBC%segment(l_seg)%open) then
            if (OBC%segment(l_seg)%direction == OBC_DIRECTION_E) then
              ps_ocn = ice_free(i,j)
              ps_ice = ice_cover(i,j)
            else
              ps_ocn = ice_free(i+1,j)
              ps_ice = ice_cover(i+1,j)
            endif
        endif ; endif
        IOF%flux_u_ocn(I,j) = IOF%flux_u_ocn(I,j) + &
            (ps_ocn * windstr_x_water(I,j) + ps_ice * str_ice_oce_x(I,j))
      enddo ; enddo
    else
      !$OMP parallel do default(shared) private(ps_ocn, ps_ice)
      do j=jsc,jec ; do I=Isc-1,iec
        ps_ocn = 1.0 ; ps_ice = 0.0
        if (G%mask2dCu(I,j)>0.0) then
          ps_ocn = 0.5*(ice_free(i+1,j) + ice_free(i,j))
          ps_ice = 0.5*(ice_cover(i+1,j) + ice_cover(i,j))
        endif
        IOF%flux_u_ocn(I,j) = IOF%flux_u_ocn(I,j) + &
            (ps_ocn * windstr_x_water(I,j) + ps_ice * str_ice_oce_x(I,j))
      enddo ; enddo
    endif
    if (local_open_v_BC) then
      !$OMP parallel do default(shared) private(ps_ocn, ps_ice)
      do J=jsc-1,jec ; do i=isc,iec
        l_seg = OBC%segnum_v(i,J)
        ps_ocn = 1.0 ; ps_ice = 0.0
        if (G%mask2dCv(i,J)>0.0) then
          ps_ocn = 0.5*(ice_free(i,j+1) + ice_free(i,j))
          ps_ice = 0.5*(ice_cover(i,j+1) + ice_cover(i,j))
        endif
        if (l_seg /= OBC_NONE) then
          if (OBC%segment(l_seg)%open) then
            if (OBC%segment(l_seg)%direction == OBC_DIRECTION_N) then
              ps_ocn = ice_free(i,j)
              ps_ice = ice_cover(i,j)
            else
              ps_ocn = ice_free(i,j+1)
              ps_ice = ice_cover(i,j+1)
            endif
        endif ; endif
        IOF%flux_v_ocn(i,J) = IOF%flux_v_ocn(i,J) + &
            (ps_ocn * windstr_y_water(i,J) + ps_ice * str_ice_oce_y(i,J))
      enddo ; enddo
    else
      !$OMP parallel do default(shared) private(ps_ocn, ps_ice)
      do J=jsc-1,jec ; do i=isc,iec
        ps_ocn = 1.0 ; ps_ice = 0.0
        if (G%mask2dCv(i,J)>0.0) then
          ps_ocn = 0.5*(ice_free(i,j+1) + ice_free(i,j))
          ps_ice = 0.5*(ice_cover(i,j+1) + ice_cover(i,j))
        endif
        IOF%flux_v_ocn(i,J) = IOF%flux_v_ocn(i,J) + &
            (ps_ocn * windstr_y_water(i,J) + ps_ice * str_ice_oce_y(i,J))
      enddo ; enddo
    endif
  else
    call SIS_error(FATAL, "set_ocean_top_stress_C2: Unrecognized flux_uv_stagger.")
  endif

  IOF%stress_count = IOF%stress_count + 1

end subroutine set_ocean_top_stress_C2

end module MOM_SIS_set_ocean_top_stress 
