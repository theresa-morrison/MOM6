module MOM_SIS_continuity

use ice_grid,          only : ice_grid_type
use MOM_cpu_clock,     only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_obsolete_params, only : obsolete_logical
use MOM_unit_scaling,  only : unit_scale_type
use SIS_diag_mediator, only : time_type, SIS_diag_ctrl
use MOM_SIS_hor_grid,      only : SIS_hor_grid_type

implicit none ; private

#include <SIS2_memory.h>

public summed_continuity, SIS_continuity_CS  

integer :: id_clock_update  !< A CPU time clock ID
integer :: id_clock_correct !< A CPU time clock ID

!> The control structure with parameters regulating the continuity solver
type :: SIS_continuity_CS ; ! private
  type(SIS_diag_ctrl), pointer :: diag => NULL() !< A structure that is used to regulate the
                             !! timing of diagnostic output.
  logical :: use_upwind2d    !< If true, use the non-split upwind scheme that was
                             !! used in older versions of SIS.
  logical :: upwind_1st      !< If true, use a directionally-split first-order upwind scheme.
  logical :: monotonic       !< If true, use the Colella & Woodward monotonic limiter;
                             !! otherwise use a simple positive definite limiter.
  logical :: simple_2nd      !< If true, use a simple second order (arithmetic mean) interpolation
                             !! of the edge values instead of the higher order interpolation.
  logical :: vol_CFL         !< If true, use the ratio of the open face lengths to the tracer
                             !! cell areas when estimating CFL numbers.
  real :: h_neglect_cont     !< The category ice mass per ocean cell area below which the
                             !! transport within this thickness category of out of a cell is
                             !! set to zero [R Z ~> kg m-2]
  real :: frac_neglect       !< When the total fluxes are distributed between categories, any
                             !! category whose ice is less than this fraction of the total mass
                             !! contributes no flux [nondim]
end type SIS_continuity_CS

!> This type is used to specify the active loop bounds
type :: loop_bounds_type ; private
  !>@{ The active index range
  integer :: ish, ieh, jsh, jeh
  !!@}
end type loop_bounds_type


contains

!> summed_continuity time steps the total ice, water, and snow mass changes summed across all the
!! thickness categories due to advection, using a monotonically limited, directionally split PPM
!! scheme or simple upwind 2-d scheme.  It may also update the ice thickness, using fluxes that are
!! proportional to the total fluxes times the ice mass divided by the total mass in the upwind cell.
subroutine summed_continuity(u, v, h_in, h, uh, vh, dt, G, US, IG, CS, h_ice)
  type(SIS_hor_grid_type),           intent(inout) :: G  !< The horizontal grid type
  type(ice_grid_type),               intent(inout) :: IG !< The sea-ice specific grid type
  real, dimension(SZIB_(G),SZJ_(G)), intent(in)    :: u  !< Zonal ice velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G)), intent(in)    :: v  !< Meridional ice velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)    :: h_in !< Initial total ice and snow mass per
                                                         !! unit cell area [R Z ~> kg m-2].
  real, dimension(SZI_(G),SZJ_(G)),  intent(inout) :: h  !< Total ice and snow mass per unit cell
                                                         !! area [R Z ~> kg m-2].
  real, dimension(SZIB_(G),SZJ_(G)), intent(out)   :: uh !< Total mass flux through zonal faces
                                                         !! = u*h*dy [R Z L2 T-1 ~> kg s-1]
  real, dimension(SZI_(G),SZJB_(G)), intent(out)   :: vh !< Total mass flux through meridional faces
                                                         !! = v*h*dx [R Z L2 T-1 ~> kg s-1]
  real,                              intent(in)    :: dt !< Time increment [T ~> s]
  type(unit_scale_type),             intent(in)    :: US !< A structure with unit conversion factors
  type(SIS_continuity_CS),           pointer       :: CS !< The control structure returned by a
                                                         !! previous call to SIS_continuity_init.
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(inout) :: h_ice  !< Total ice mass per unit cell
                                                         !! area [R Z ~> kg m-2].  h_ice must not exceed h.

  ! Local variables
  type(loop_bounds_type) :: LB  ! A structure with the active loop bounds.
  real, dimension(SZIB_(G),SZJ_(G)) :: uh_ice ! Ice mass flux through zonal faces = u*h*dy
                                              ! [R Z L2 T-1 ~> kg s-1].
  real, dimension(SZI_(G),SZJB_(G)) :: vh_ice ! Ice mass flux through meridional faces = v*h*dx
                                              ! [R Z L2 T-1 ~> kg s-1].
  real    :: h_up
  integer :: is, ie, js, je, stencil
  integer :: i, j

  logical :: x_first
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  if (.not.associated(CS)) call SIS_error(FATAL, &
         "SIS_continuity: Module must be initialized before it is used.")
  x_first = (MOD(G%first_direction,2) == 0)

  stencil = 3 ; if (CS%simple_2nd) stencil = 2 ; if (CS%upwind_1st) stencil = 1

  do j=js,je ; do i=is,ie ; if (h_in(i,j) < 0.0) then
    call SIS_error(FATAL, 'Negative mass input to summed_continuity().')
  endif ; enddo ; enddo

  if (present(h_ice)) then ; do j=js,je ; do i=is,ie ; if (h_ice(i,j) > h_in(i,j)) then
    call SIS_error(FATAL, 'ice mass exceeds total mass in summed_continuity().')
  endif ; enddo ; enddo ; endif

  if (CS%use_upwind2d) then
    ! This reproduces the scheme that was originally used in SIS1.
    !$OMP parallel default(shared) private(h_up)
    !$OMP do
    do j=js,je ; do I=is-1,ie
      if (u(I,j) >= 0.0) then ; h_up = h_in(i,j)
      else ; h_up = h_in(i+1,j) ; endif
      if (h_up < IG%CatIce*CS%h_neglect_cont) h_up = 0.0
      uh(I,j) = G%dy_Cu(I,j) * u(I,j) * h_up
    enddo ; enddo
    !$OMP do
    do J=js-1,je ; do i=is,ie
      if (v(i,J) >= 0.0) then ; h_up = h_in(i,j)
      else ; h_up = h_in(i,j+1) ; endif
      if (h_up < IG%CatIce*CS%h_neglect_cont) h_up = 0.0
      vh(i,J) = G%dx_Cv(i,J) * v(i,J) * h_up
    enddo ; enddo
    if (present(h_ice)) then
      !$OMP do
      do j=js,je ; do I=is-1,ie
        if (uh(I,j) < 0.0) then ; uh_ice(I,j) = uh(I,j) * (h_ice(i+1,j) / h_in(i+1,j))
        elseif (uh(I,j) > 0.0) then ; uh_ice(I,j) = uh(I,j) * (h_ice(i,j) / h_in(i,j))
        else ; uh_ice(I,j) = 0.0 ; endif
      enddo ; enddo
      !$OMP do
      do J=js-1,je ; do i=is,ie
        if (vh(i,J) < 0.0) then ; vh_ice(i,J) = vh(i,J) * (h_ice(i,j+1) / h_in(i,j+1))
        elseif (vh(i,J) > 0.0) then ; vh_ice(i,J) = vh(i,J) * (h_ice(i,j) / h_in(i,j))
        else ; vh_ice(i,J) = 0.0 ; endif
      enddo ; enddo
      !$OMP do
      do j=js,je ; do i=is,ie
        h_ice(i,j) = h_ice(i,j) - (dt * G%IareaT(i,j)) * &
             ((uh_ice(I,j) - uh_ice(I-1,j)) + (vh_ice(i,J) - vh_ice(i,J-1)))
      enddo ; enddo
    endif
    !$OMP do
    do j=js,je ; do i=is,ie
      h(i,j) = h_in(i,j) - (dt * G%IareaT(i,j)) * &
           ((uh(I,j) - uh(I-1,j)) + (vh(i,J) - vh(i,J-1)))
      ! if (h(i,j) < 0.0) call SIS_error(FATAL, &
      !   'Negative thickness encountered in ice_total_continuity().')
      ! if (present(h_ice)) then ; if (h_ice(i,j) > h(i,j)) then
      !   call SIS_error(FATAL, 'ice mass exceeds total mass in ice_total_continuity() 2d.')
      ! endif ; endif
    enddo ; enddo
    !$OMP end parallel
  elseif (x_first) then
    ! First, advect zonally.
    LB%ish = G%isc ; LB%ieh = G%iec ; LB%jsh = G%jsc-stencil ; LB%jeh = G%jec+stencil
    call zonal_mass_flux(u, dt, G, US, IG, CS, LB, htot_in=h_in, uh_tot=uh, h_mobilize=CS%h_neglect_cont)

    call cpu_clock_begin(id_clock_update)

    if (present(h_ice)) then
      !$OMP parallel do default(shared)
      do j=LB%jsh,LB%jeh
        do I=LB%ish-1,LB%ieh
          if (uh(I,j) < 0.0) then ; uh_ice(I,j) = uh(I,j) * (h_ice(i+1,j) / h_in(i+1,j))
          elseif (uh(I,j) > 0.0) then ; uh_ice(I,j) = uh(I,j) * (h_ice(i,j) / h_in(i,j))
          else ; uh_ice(I,j) = 0.0 ; endif
        enddo
        do i=LB%ish,LB%ieh
          h_ice(i,j) = h_ice(i,j) - (dt * G%IareaT(i,j)) * (uh_ice(I,j) - uh_ice(I-1,j))
        enddo
      enddo
    endif

    !$OMP parallel do default(shared)
    do j=LB%jsh,LB%jeh ; do i=LB%ish,LB%ieh
      h(i,j) = h_in(i,j) - (dt * G%IareaT(i,j)) * (uh(I,j) - uh(I-1,j))
      ! if (h(i,j) < 0.0) call SIS_error(FATAL, &
      !   'Negative thickness encountered in u-pass of ice_total_continuity().')
      ! if (present(h_ice)) then ; if (h_ice(i,j) > h(i,j)) then
      !   call SIS_error(FATAL, 'ice mass exceeds total mass in ice_total_continuity() x-1.')
      ! endif ; endif
    enddo ; enddo
    call cpu_clock_end(id_clock_update)

    LB%ish = G%isc ; LB%ieh = G%iec ; LB%jsh = G%jsc ; LB%jeh = G%jec
    ! Now advect meridionally, using the updated thicknesses to determine the fluxes.
    call meridional_mass_flux(v, dt, G, US, IG, CS, LB, htot_in=h, vh_tot=vh, h_mobilize=CS%h_neglect_cont)

    call cpu_clock_begin(id_clock_update)
    if (present(h_ice)) then
      !$OMP parallel do default(shared)
      do J=LB%jsh-1,LB%jeh ; do i=LB%ish,LB%ieh
        if (vh(i,J) < 0.0) then ; vh_ice(i,J) = vh(i,J) * (h_ice(i,j+1) / h(i,j+1))
        elseif (vh(i,J) > 0.0) then ; vh_ice(i,J) = vh(i,J) * (h_ice(i,j) / h(i,j))
        else ; vh_ice(i,J) = 0.0 ; endif
      enddo ; enddo
      !$OMP parallel do default(shared)
      do j=LB%jsh,LB%jeh ; do i=LB%ish,LB%ieh
        h_ice(i,j) = h_ice(i,j) - (dt * G%IareaT(i,j)) * (vh_ice(i,J) - vh_ice(i,J-1))
      enddo ; enddo
    endif

    !$OMP parallel do default(shared)
    do j=LB%jsh,LB%jeh ; do i=LB%ish,LB%ieh
      h(i,j) = h(i,j) - (dt * G%IareaT(i,j)) * (vh(i,J) - vh(i,J-1))
      if (h(i,j) < 0.0) call SIS_error(FATAL, &
        'Negative thickness encountered in v-pass of summed_continuity().')
      ! if (present(h_ice)) then ; if (h_ice(i,j) > h(i,j)) then
      !   call SIS_error(FATAL, 'ice mass exceeds total mass in ice_total_continuity() x-2.')
      ! endif ; endif
    enddo ; enddo
    call cpu_clock_end(id_clock_update)

  else  ! .not. x_first
    !  First, advect meridionally, so set the loop bounds accordingly.
    LB%ish = G%isc-stencil ; LB%ieh = G%iec+stencil ; LB%jsh = G%jsc ; LB%jeh = G%jec
    call meridional_mass_flux(v, dt, G, US, IG, CS, LB, htot_in=h_in, vh_tot=vh, h_mobilize=CS%h_neglect_cont)

    call cpu_clock_begin(id_clock_update)
    if (present(h_ice)) then
      !$OMP parallel do default(shared)
      do J=LB%jsh-1,LB%jeh ; do i=LB%ish,LB%ieh
        if (vh(i,J) < 0.0) then ; vh_ice(i,J) = vh(i,J) * (h_ice(i,j+1) / h_in(i,j+1))
        elseif (vh(i,J) > 0.0) then ; vh_ice(i,J) = vh(i,J) * (h_ice(i,j) / h_in(i,j))
        else ; vh_ice(i,J) = 0.0 ; endif
      enddo ; enddo
      !$OMP parallel do default(shared)
      do j=LB%jsh,LB%jeh ; do i=LB%ish,LB%ieh
        h_ice(i,j) = h_ice(i,j) - (dt * G%IareaT(i,j)) * (vh_ice(i,J) - vh_ice(i,J-1))
      enddo ; enddo
    endif

    !$OMP parallel do default(shared)
    do j=LB%jsh,LB%jeh ; do i=LB%ish,LB%ieh
      h(i,j) = h_in(i,j) - (dt * G%IareaT(i,j)) * (vh(i,J) - vh(i,J-1))
      ! if (h(i,j) < 0.0) call SIS_error(FATAL, &
      !   'Negative thickness encountered in v-pass of ice_total_continuity().')
      ! if (present(h_ice)) then ; if (h_ice(i,j) > h(i,j)) then
      !   call SIS_error(FATAL, 'ice mass exceeds total mass in ice_total_continuity() y-1.')
      ! endif ; endif
    enddo ; enddo
    call cpu_clock_end(id_clock_update)

    !  Now advect zonally, using the updated thicknesses to determine the fluxes.
    LB%ish = G%isc ; LB%ieh = G%iec ; LB%jsh = G%jsc ; LB%jeh = G%jec
    call zonal_mass_flux(u, dt, G, US, IG, CS, LB, htot_in=h, uh_tot=uh, h_mobilize=CS%h_neglect_cont)

    call cpu_clock_begin(id_clock_update)

    if (present(h_ice)) then
      !$OMP parallel do default(shared)
      do j=LB%jsh,LB%jeh
        do I=LB%ish-1,LB%ieh
          if (uh(I,j) < 0.0) then ; uh_ice(I,j) = uh(I,j) * (h_ice(i+1,j) / h(i+1,j))
          elseif (uh(I,j) > 0.0) then ; uh_ice(I,j) = uh(I,j) * (h_ice(i,j) / h(i,j))
          else ; uh_ice(I,j) = 0.0 ; endif
        enddo
        do i=LB%ish,LB%ieh
          h_ice(i,j) = h_ice(i,j) - (dt * G%IareaT(i,j)) * (uh_ice(I,j) - uh_ice(I-1,j))
        enddo
      enddo
    endif

    !$OMP parallel do default(shared)
    do j=LB%jsh,LB%jeh ; do i=LB%ish,LB%ieh
      h(i,j) = h(i,j) - (dt * G%IareaT(i,j)) * (uh(I,j) - uh(I-1,j))
      if (h(i,j) < 0.0) call SIS_error(FATAL, &
        'Negative thickness encountered in u-pass of summed_continuity().')
      ! if (present(h_ice)) then ; if (h_ice(i,j) > h(i,j)) then
      !   call SIS_error(FATAL, 'ice mass exceeds total mass in ice_total_continuity() y-2.')
      ! endif ; endif
    enddo ; enddo
    call cpu_clock_end(id_clock_update)

  endif  ! End of x_first block.

end subroutine summed_continuity

!> Calculates the mass or volume fluxes through the zonal
!! faces, and other related quantities.
subroutine zonal_mass_flux(u, dt, G, US, IG, CS, LB, h_in, uh, htot_in, uh_tot, h_mobilize, &
                           masking_uh, masking_uhtot)
  type(SIS_hor_grid_type), intent(inout) :: G   !< The horizontal grid type
  type(ice_grid_type),     intent(inout) :: IG  !< The sea-ice specific grid type
  real, dimension(SZIB_(G),SZJ_(G)), &
                           intent(in)    :: u   !< Zonal ice velocity [L T-1 ~> m s-1].
  real,                    intent(in)    :: dt  !< Time increment [T ~> s]
  type(unit_scale_type),   intent(in)    :: US  !< A structure with unit conversion factors
  type(SIS_continuity_CS), pointer       :: CS  !< The control structure returned by a
                                                !! previous call to SIS_continuity_init.
  type(loop_bounds_type),  intent(in)    :: LB  !< A structure with the active loop bounds.
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)), &
                 optional, intent(in)    :: h_in !< Category thickness used to calculate the fluxes [R Z ~> kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZCAT_(IG)), &
                 optional, intent(out)   :: uh  !< Category volume flux through zonal faces = u*h*dy
                                                !! [R Z L2 T-1 ~> kg s-1].
  real, dimension(SZI_(G),SZJ_(G)), &
                 optional, intent(in)    :: htot_in !< Total thicknesses used to calculate the fluxes [R Z ~> kg m-2].
  real, dimension(SZIB_(G),SZJ_(G)), &
                 optional, intent(out)   :: uh_tot !< Total mass flux through zonal faces = u*htot*dy
                                                !! [R Z L2 T-1 ~> kg s-1].
  real,          optional, intent(in)    :: h_mobilize !< The minimum ice thickness per category that
                                                !! is able to move out of a cell [R Z ~> kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZCAT_(IG)), &
                 optional, intent(in)    :: masking_uh  !< If this is 0, uh = 0.  Often this is another
                                                !! zonal mass flux, e.g. of ice [R Z L2 T-1 ~> kg s-1].
  real, dimension(SZIB_(G),SZJ_(G)), &
                 optional, intent(in)    :: masking_uhtot !< If this is 0, uh_tot = 0.  Often this is another
                                                !! total zonal mass flux, e.g. of ice [R Z L2 T-1 ~> kg s-1].
!   This subroutine calculates the mass or volume fluxes through the zonal
! faces, and other related quantities.

  ! Local variables
!  real, dimension(SZIB_(G)) :: &
!    duhdu      ! Partial derivative of uh with u [R Z L ~> kg m-1].
  real, dimension(SZI_(G),SZJ_(G)) :: &
    htot, &    ! The total thickness summed across categories [R Z ~> kg m-2].
    I_htot, &  ! The inverse of htot or 0 [R-1 Z-1 ~> m2 kg-1].
    hl, hr      ! Left and right face thicknesses [R Z ~> kg m-2].
  real, dimension(SZIB_(G)) :: &
    uhtot      ! The total transports [R Z L2 T-1 ~> kg s-1].
  real :: CFL  ! The CFL number based on the local velocity and grid spacing [nondim].
  real :: curv_3 ! A measure of the thickness curvature over a grid length,
                 ! with the same units as h_in.
!  real :: h_marg ! The marginal thickness of a flux [R Z ~> kg m-2].
!  real :: dx_E, dx_W ! Effective x-grid spacings to the east and west [L ~> m].
  integer :: i, j, k, ish, ieh, jsh, jeh, nCat

  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh ; nCat = IG%CatIce

  call cpu_clock_begin(id_clock_update)

  htot(:,:) = 0.0
  if (present(htot_in)) then
    if (present(h_mobilize)) then
      !$OMP parallel do default(shared)
      do j=jsh,jeh ; do i=G%isd,G%ied
        if (htot_in(i,j) >= nCat*h_mobilize) htot(i,j) = htot(i,j) + htot_in(i,j)
      enddo ; enddo
    else
      !$OMP parallel do default(shared)
      do j=jsh,jeh ; do i=G%isd,G%ied
        htot(i,j) = htot(i,j) + htot_in(i,j)
      enddo ; enddo
    endif
  elseif (present(h_in)) then
    if (present(h_mobilize)) then
      !$OMP parallel do default(shared)
      do j=jsh,jeh ; do k=1,nCat ; do i=G%isd,G%ied
        if (h_in(i,j,k) >= h_mobilize) htot(i,j) = htot(i,j) + h_in(i,j,k)
      enddo ; enddo ; enddo
    else
      !$OMP parallel do default(shared)
      do j=jsh,jeh ; do k=1,nCat ; do i=G%isd,G%ied
        htot(i,j) = htot(i,j) + h_in(i,j,k)
      enddo ; enddo ; enddo
    endif
  else
    call SIS_error(FATAL, "Either h_in or htot_in must be present in call to zonal_mass_flux.")
  endif

  ! This sets hl and hr.
  if (CS%upwind_1st) then
    do j=jsh,jeh ; do i=ish-1,ieh+1
      hl(i,j) = htot(i,j) ; hr(i,j) = htot(i,j)
    enddo ; enddo
  else
    call PPM_reconstruction_x(htot, hl, hr, G, LB, 0.0, CS%monotonic, &
                              simple_2nd=CS%simple_2nd)
  endif
  call cpu_clock_end(id_clock_update)

  call cpu_clock_begin(id_clock_correct)

  !$OMP parallel do default(shared) private(uhtot)
  do j=jsh,jeh
    ! Set uhtot and duhdu.
    do I=ish-1,ieh
      ! Set new values of uh and duhdu.
      if (u(I,j) > 0.0) then
        if (CS%vol_CFL) then ; CFL = (u(I,j) * dt) * (G%dy_Cu(I,j) * G%IareaT(i,j))
        else ; CFL = u(I,j) * dt * G%IdxT(i,j) ; endif
        curv_3 = hL(i,j) + hR(i,j) - 2.0*htot(i,j)
        uhtot(I) = G%dy_Cu(I,j) * u(I,j) * &
            (hR(i,j) + CFL * (0.5*(hL(i,j) - hR(i,j)) + curv_3*(CFL - 1.5)))
!        h_marg = hR(i,j) + CFL * ((hL(i,j) - hR(i,j)) + 3.0*curv_3*(CFL - 1.0))
      elseif (u(I,j) < 0.0) then
        if (CS%vol_CFL) then ; CFL = (-u(I,j) * dt) * (G%dy_Cu(I,j) * G%IareaT(i+1,j))
        else ; CFL = -u(I,j) * dt * G%IdxT(i+1,j) ; endif
        curv_3 = hL(i+1,j) + hR(i+1,j) - 2.0*htot(i+1,j)
        uhtot(I) = G%dy_Cu(I,j) * u(I,j) * &
            (hL(i+1,j) + CFL * (0.5*(hR(i+1,j)-hL(i+1,j)) + curv_3*(CFL - 1.5)))
!        h_marg = hL(i+1) + CFL * ((hR(i+1,j)-hL(i+1,j)) + 3.0*curv_3*(CFL - 1.0))
      else
        uhtot(I) = 0.0
!        h_marg = 0.5 * (hl(i+1,j) + hr(i,j))
      endif
!      duhdu(I,j) = G%dy_Cu(I,j) * h_marg ! * visc_rem(I)
    enddo

    ! Partition the transports by category in proportion to their relative masses.
    if (present(uh)) then
      do i=ish-1,ieh+1
        I_htot(i,j) = 0.0 ; if (htot(i,j) > 0.0) I_htot(i,j) = 1.0 / htot(i,j)
      enddo
      if (present(h_mobilize)) then
        do k=1,nCat ; do I=ish-1,ieh
          uh(I,j,k) = 0.0
          if (u(I,j) >= 0.0) then
            if (h_in(i,j,k) >= h_mobilize) uh(I,j,k) = uhtot(I) * (h_in(i,j,k) * I_htot(i,j))
          else
            if (h_in(i+1,j,k) >= h_mobilize) uh(I,j,k) = uhtot(I) * (h_in(i+1,j,k) * I_htot(i+1,j))
          endif
        enddo ; enddo
      else
        do k=1,nCat ; do I=ish-1,ieh
          if (u(I,j) >= 0.0) then
            uh(I,j,k) = uhtot(I) * (h_in(i,j,k) * I_htot(i,j))
          else
            uh(I,j,k) = uhtot(I) * (h_in(i+1,j,k) * I_htot(i+1,j))
          endif
        enddo ; enddo
      endif
    endif

    ! Block mass fluxes in categories where a related flux (e.g. of ice) is zero.
    if (present(masking_uh) .and. present(uh)) then
      do k=1,nCat ; do I=ish-1,ieh
        if (masking_uh(I,j,k) == 0.0) uh(I,j,k) = 0.0
        if (abs(uh(I,j,k)) < CS%frac_neglect*abs(masking_uh(I,j,k))) uh(I,j,k) = 0.0
      enddo ; enddo
    endif

    if (present(uh_tot) .and. present(masking_uhtot)) then
      do I=ish-1,ieh
        uh_tot(I,j) = uhtot(I)
        if (masking_uhtot(I,j) == 0.0) uh_tot(I,j) = 0.0
      enddo
    elseif (present(uh_tot)) then
      do I=ish-1,ieh
        uh_tot(I,j) = uhtot(I)
      enddo
    endif

  enddo ! j-loop
  call cpu_clock_end(id_clock_correct)

end subroutine zonal_mass_flux

!> Calculates the mass or volume fluxes through the meridional
!! faces, and other related quantities.
subroutine meridional_mass_flux(v, dt, G, US, IG, CS, LB, h_in, vh, htot_in, vh_tot, &
                                h_mobilize, masking_vh, masking_vhtot)
  type(SIS_hor_grid_type), intent(inout) :: G   !< The horizontal grid type
  type(ice_grid_type),     intent(inout) :: IG  !< The sea-ice specific grid type
  real, dimension(SZI_(G),SZJB_(G)), &
                           intent(in)    :: v   !< Meridional ice velocity [L T-1 ~> m s-1].
  real,                    intent(in)    :: dt  !< Time increment [T ~> s]
  type(unit_scale_type),   intent(in)    :: US  !< A structure with unit conversion factors
  type(SIS_continuity_CS), pointer       :: CS  !< The control structure returned by a
                                                !! previous call to SIS_continuity_init.
  type(loop_bounds_type),  intent(in)    :: LB  !< A structure with the active loop bounds.
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)), &
                 optional, intent(in)    :: h_in !< Category thickness used to calculate the fluxes [R Z ~> kg m-2].
  real, dimension(SZI_(G),SZJB_(G),SZCAT_(IG)), &
                 optional, intent(out)   :: vh  !< Category volume flux through meridional faces = v*h*dx
                                                !! [R Z L2 T-1 ~> kg s-1].
  real, dimension(SZI_(G),SZJ_(G)), &
                 optional, intent(in)    :: htot_in !< Total thicknesses used to calculate the fluxes [R Z ~> kg m-2].
  real, dimension(SZI_(G),SZJB_(G)), &
                 optional, intent(out)   :: vh_tot !< Total mass flux through meridional faces = v*htot*dx
                                                !! [R Z L2 T-1 ~> kg s-1].
  real,          optional, intent(in)    :: h_mobilize !< The minimum ice thickness per category that
                                                !! is able to move out of a cell [R Z ~> kg m-2].
  real, dimension(SZI_(G),SZJB_(G),SZCAT_(IG)), &
                 optional, intent(in)    :: masking_vh  !< If this is 0, vh = 0.  Often this is another
                                                !! meridional mass flux, e.g. of ice [R Z L2 T-1 ~> kg s-1].
  real, dimension(SZI_(G),SZJB_(G)), &
                 optional, intent(in)    :: masking_vhtot !< If this is 0, vh_tot = 0.  Often this is another
                                                !! total meridional mass flux, e.g. of ice [R Z L2 T-1 ~> kg s-1].

!   This subroutine calculates the mass or volume fluxes through the meridional
! faces, and other related quantities.

  ! Local variables
  real, dimension(SZI_(G)) :: &
    dvhdv      ! Partial derivative of vh with v [R Z L ~> kg m-1].
  real, dimension(SZI_(G),SZJ_(G)) :: &
    htot, &    ! The total thickness summed across categories [R Z ~> kg m-2].
    I_htot, &  ! The inverse of htot or 0 [R-1 Z-1 ~> m2 kg-1].
    hl, hr     ! Left and right face thicknesses [R Z ~> kg m-2].
  real, dimension(SZI_(G)) :: &
    vhtot      ! The total transports [R Z L2 T-1 ~> kg s-1].
  real :: CFL ! The CFL number based on the local velocity and grid spacing [nondim].
  real :: curv_3 ! A measure of the thickness curvature over a grid length,
                 ! with the same units as h_in.
  real :: h_marg ! The marginal thickness of a flux [R Z ~> kg m-2].
!  real :: dy_N, dy_S ! Effective y-grid spacings to the north and south [L ~> m].
  integer :: i, j, k, ish, ieh, jsh, jeh, nCat

  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh ; nCat = IG%CatIce

  call cpu_clock_begin(id_clock_update)

  htot(:,:) = 0.0

  if (present(htot_in)) then
    if (present(h_mobilize)) then
      !$OMP parallel do default(shared)
      do j=G%jsd,G%jed ; do i=ish,ieh
        if (htot_in(i,j) >= nCat*h_mobilize) htot(i,j) = htot(i,j) + htot_in(i,j)
      enddo ; enddo
    else
      !$OMP parallel do default(shared)
      do j=G%jsd,G%jed ; do i=ish,ieh
        htot(i,j) = htot(i,j) + htot_in(i,j)
      enddo ; enddo
    endif
  elseif (present(h_in)) then
    if (present(h_mobilize)) then
      !$OMP parallel do default(shared)
      do j=G%jsd,G%jed ; do k=1,nCat ; do i=ish,ieh
        if (h_in(i,j,k) >= h_mobilize) htot(i,j) = htot(i,j) + h_in(i,j,k)
      enddo ; enddo ; enddo
    else
      !$OMP parallel do default(shared)
      do j=G%jsd,G%jed ; do k=1,nCat ; do i=ish,ieh
        htot(i,j) = htot(i,j) + h_in(i,j,k)
      enddo ; enddo ; enddo
    endif
  else
    call SIS_error(FATAL, "Either h_in or htot_in must be present in call to meridional_mass_flux.")
  endif
  if (present(vh)) then ; do j=jsh-1,jeh+1 ; do i=ish,ieh
    I_htot(i,j) = 0.0 ; if (htot(i,j) > 0.0) I_htot(i,j) = 1.0 / htot(i,j)
  enddo ; enddo ; endif

  ! This sets hl and hr.
  if (CS%upwind_1st) then
    do j=jsh-1,jeh+1 ; do i=ish,ieh
      hl(i,j) = htot(i,j) ; hr(i,j) = htot(i,j)
    enddo ; enddo
  else
    call PPM_reconstruction_y(htot, hl, hr, G, LB, 0.0, CS%monotonic, &
                              simple_2nd=CS%simple_2nd)
  endif
  call cpu_clock_end(id_clock_update)

  call cpu_clock_begin(id_clock_correct)
  !$OMP parallel do default(shared) private(vhtot)
  do J=jsh-1,jeh
    ! This sets vh and dvhdv.
    do i=ish,ieh
      if (v(i,J) > 0.0) then
        if (CS%vol_CFL) then ; CFL = (v(i,J) * dt) * (G%dx_Cv(i,J) * G%IareaT(i,j))
        else ; CFL = v(i,J) * dt * G%IdyT(i,j) ; endif
        curv_3 = hL(i,j) + hR(i,j) - 2.0*htot(i,j)
        vhtot(i) = G%dx_Cv(i,J) * v(i,J) * ( hR(i,j) + CFL * &
            (0.5*(hL(i,j) - hR(i,j)) + curv_3*(CFL - 1.5)) )
       ! h_marg = hR(i,j) + CFL * ((hL(i,j) - hR(i,j)) + 3.0*curv_3*(CFL - 1.0))
      elseif (v(i,J) < 0.0) then
        if (CS%vol_CFL) then ; CFL = (-v(i,J) * dt) * (G%dx_Cv(i,J) * G%IareaT(i,j+1))
        else ; CFL = -v(i,J) * dt *  G%IdyT(i,j+1) ; endif
        curv_3 = hL(i,j+1) + hR(i,j+1) - 2.0*htot(i,j+1)
        vhtot(i) = G%dx_Cv(i,J) * v(i,J) * ( hL(i,j+1) + CFL * &
            (0.5*(hR(i,j+1)-hL(i,j+1)) + curv_3*(CFL - 1.5)) )
       ! h_marg = hL(i,j+1) + CFL * ((hR(i,j+1)-hL(i,j+1)) + 3.0*curv_3*(CFL - 1.0))
      else
        vhtot(i) = 0.0
        ! h_marg = 0.5 * (hl(i,j+1) + hr(i,j))
      endif
      ! dvhdv(i) = G%dx_Cv(i,J) * h_marg ! * visc_rem(i)
    enddo

    ! Partition the transports by category in proportion to their relative masses.
    if (present(vh)) then
      if (present(h_mobilize)) then
        do k=1,nCat ; do i=ish,ieh
          vh(i,J,k) = 0.0
          if (v(i,J) >= 0.0) then
            if (h_in(i,j,k) >= h_mobilize) vh(i,J,k) = vhtot(i) * (h_in(i,j,k) * I_htot(i,j))
          else
            if (h_in(i,j+1,k) >= h_mobilize) vh(i,J,k) = vhtot(i) * (h_in(i,j+1,k) * I_htot(i,j+1))
          endif
        enddo ; enddo
      else
        do k=1,nCat ; do i=ish,ieh
          if (v(i,J) >= 0.0) then
            vh(i,J,k) = vhtot(i) * (h_in(i,j,k) * I_htot(i,j))
          else
            vh(i,J,k) = vhtot(i) * (h_in(i,j+1,k) * I_htot(i,j+1))
          endif
        enddo ; enddo
      endif
    endif

    ! Block mass fluxes in categories where a related flux (e.g. of ice) is zero.
    if (present(masking_vh) .and. present(vh)) then
      do k=1,nCat ; do i=ish,ieh
        if (masking_vh(i,J,k) == 0.0) vh(i,J,k) = 0.0
        if (abs(vh(i,J,k)) < CS%frac_neglect*abs(masking_vh(i,J,k))) vh(i,J,k) = 0.0
      enddo ; enddo
    endif

    if (present(vh_tot) .and. present(masking_vhtot)) then
      do i=ish,ieh
        vh_tot(i,J) = vhtot(I)
        if (masking_vhtot(i,J) == 0.0) vh_tot(i,J) = 0.0
      enddo
    elseif (present(vh_tot)) then
      do i=ish,ieh
        vh_tot(i,J) = vhtot(i)
      enddo
    endif

  enddo ! j-loop
  call cpu_clock_end(id_clock_correct)

end subroutine meridional_mass_flux

!> Calculate a piecewise parabolic thickness reconstruction in the x-direction.
subroutine PPM_reconstruction_x(h_in, h_l, h_r, G, LB, h_min, monotonic, simple_2nd)
  type(SIS_hor_grid_type),          intent(in)  :: G   !< The horizontal grid type
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: h_in !< Initial thickness of a category [R Z ~> kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: h_l !< Left edge value of thickness reconstruction [R Z ~> kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: h_r !< Right edge value of thickness reconstruction [R Z ~> kg m-2]
  type(loop_bounds_type),           intent(in)  :: LB  !< A structure with the active loop bounds.
  real,                             intent(in)  :: h_min !< The minimum thickness that can be
                                                       !! obtained by a concave parabolic fit [R Z ~> kg m-2].
  logical, optional,                intent(in)  :: monotonic !< If true, use the Colella & Woodward monotonic limiter.
                                                       !! Otherwise use a simple positive-definite limiter.
  logical, optional,                intent(in)  :: simple_2nd !< If true, use the arithmetic mean thicknesses as the
                                                       !! default edge values for a simple 2nd order scheme.
! This subroutine calculates left/right edge values for PPM reconstruction.

  ! Local variables with useful mnemonic names.
  real, dimension(SZI_(G),SZJ_(G))  :: slp ! The slopes.
  real, parameter :: oneSixth = 1./6.
  real :: h_ip1, h_im1
  real :: dMx, dMn
  logical :: use_CW84, use_2nd
  character(len=256) :: mesg
  integer :: i, j, isl, iel, jsl, jel, stencil

  use_CW84 = .false. ; if (present(monotonic)) use_CW84 = monotonic
  use_2nd = .false. ; if (present(simple_2nd)) use_2nd = simple_2nd
  isl = LB%ish-1 ; iel = LB%ieh+1 ; jsl = LB%jsh ; jel = LB%jeh

  ! This is the stencil of the reconstruction, not the scheme overall.
  stencil = 2 ; if (use_2nd) stencil = 1

  if ((isl-stencil < G%isd) .or. (iel+stencil > G%ied)) then
    write(mesg,'("In SIS_continuity, PPM_reconstruction_x called with a ", &
               & "x-halo that needs to be increased by ",i2,".")') &
               stencil + max(G%isd-isl,iel-G%ied)
    call SIS_error(FATAL,mesg)
  endif
  if ((jsl < G%jsd) .or. (jel > G%jed)) then
    write(mesg,'("In SIS_continuity, PPM_reconstruction_x called with a ", &
               & "y-halo that needs to be increased by ",i2,".")') &
               max(G%jsd-jsl,jel-G%jed)
    call SIS_error(FATAL,mesg)
  endif

  if (use_2nd) then
    !$OMP parallel do default(shared) private(h_im1,h_ip1)
    do j=jsl,jel ; do i=isl,iel
      h_im1 = G%mask2dT(i-1,j) * h_in(i-1,j) + (1.0-G%mask2dT(i-1,j)) * h_in(i,j)
      h_ip1 = G%mask2dT(i+1,j) * h_in(i+1,j) + (1.0-G%mask2dT(i+1,j)) * h_in(i,j)
      h_l(i,j) = 0.5*( h_im1 + h_in(i,j) )
      h_r(i,j) = 0.5*( h_ip1 + h_in(i,j) )
    enddo ; enddo
  else
    !$OMP parallel do default(shared) private(dMx,dMn,h_im1,h_ip1)
    do j=jsl,jel
      do i=isl-1,iel+1
        if ((G%mask2dT(i-1,j) * G%mask2dT(i,j) * G%mask2dT(i+1,j)) == 0.0) then
          slp(i,j) = 0.0
        else
          ! This uses a simple 2nd order slope.
          slp(i,j) = 0.5 * (h_in(i+1,j) - h_in(i-1,j))
          ! Monotonic constraint, see Eq. B2 in Lin 1994, MWR (132)
          dMx = max(h_in(i+1,j), h_in(i-1,j), h_in(i,j)) - h_in(i,j)
          dMn = h_in(i,j) - min(h_in(i+1,j), h_in(i-1,j), h_in(i,j))
          slp(i,j) = sign(1.,slp(i,j)) * min(abs(slp(i,j)), 2. * min(dMx, dMn))
                  ! * (G%mask2dT(i-1,j) * G%mask2dT(i,j) * G%mask2dT(i+1,j))
        endif
      enddo

      do i=isl,iel
        ! Neighboring values should take into account any boundaries.  The 3
        ! following sets of expressions are equivalent.
      ! h_im1 = h_in(i-1,j,k) ; if (G%mask2dT(i-1,j) < 0.5) h_im1 = h_in(i,j)
      ! h_ip1 = h_in(i+1,j,k) ; if (G%mask2dT(i+1,j) < 0.5) h_ip1 = h_in(i,j)
        h_im1 = G%mask2dT(i-1,j) * h_in(i-1,j) + (1.0-G%mask2dT(i-1,j)) * h_in(i,j)
        h_ip1 = G%mask2dT(i+1,j) * h_in(i+1,j) + (1.0-G%mask2dT(i+1,j)) * h_in(i,j)
        ! Left/right values following Eq. B2 in Lin 1994, MWR (132)
        h_l(i,j) = 0.5*( h_im1 + h_in(i,j) ) + oneSixth*( slp(i-1,j) - slp(i,j) )
        h_r(i,j) = 0.5*( h_ip1 + h_in(i,j) ) + oneSixth*( slp(i,j) - slp(i+1,j) )
      enddo
    enddo
  endif

  if (use_CW84) then
    call PPM_limit_CW84(h_in, h_l, h_r, G, isl, iel, jsl, jel)
  else
    call PPM_limit_pos(h_in, h_l, h_r, h_min, G, isl, iel, jsl, jel)
  endif

  return
end subroutine PPM_reconstruction_x

!> Calculate a piecewise parabolic thickness reconstruction in the y-direction.
subroutine PPM_reconstruction_y(h_in, h_l, h_r, G, LB, h_min, monotonic, simple_2nd)
  type(SIS_hor_grid_type),          intent(in)  :: G   !< The horizontal grid type
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: h_in !< Initial thickness of a category [R Z ~> kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: h_l !< Left edge value of thickness reconstruction [R Z ~> kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: h_r !< Right edge value of thickness reconstruction [R Z ~> kg m-2]
  type(loop_bounds_type),           intent(in)  :: LB  !< A structure with the active loop bounds.
  real,                             intent(in)  :: h_min !< The minimum thickness that can be
                                                       !! obtained by a concave parabolic fit [R Z ~> kg m-2].
  logical, optional,                intent(in)  :: monotonic !< If true, use the Colella & Woodward monotonic limiter.
                                                       !! Otherwise use a simple positive-definite limiter.
  logical, optional,                intent(in)  :: simple_2nd !< If true, use the arithmetic mean thicknesses as the
                                                       !! default edge values for a simple 2nd order scheme.
! This subroutine calculates left/right edge values for PPM reconstruction.

! Local variables with useful mnemonic names.
  real, dimension(SZI_(G),SZJ_(G))  :: slp ! The slopes.
  real, parameter :: oneSixth = 1./6.
  real :: h_jp1, h_jm1
  real :: dMx, dMn
  logical :: use_CW84, use_2nd
  character(len=256) :: mesg
  integer :: i, j, isl, iel, jsl, jel, stencil

  use_CW84 = .false. ; if (present(monotonic)) use_CW84 = monotonic
  use_2nd = .false. ; if (present(simple_2nd)) use_2nd = simple_2nd
  isl = LB%ish ; iel = LB%ieh ; jsl = LB%jsh-1 ; jel = LB%jeh+1

  ! This is the stencil of the reconstruction, not the scheme overall.
  stencil = 2 ; if (use_2nd) stencil = 1

  if ((isl < G%isd) .or. (iel > G%ied)) then
    write(mesg,'("In SIS_continuity, PPM_reconstruction_y called with a ", &
               & "x-halo that needs to be increased by ",i2,".")') &
               max(G%isd-isl,iel-G%ied)
    call SIS_error(FATAL,mesg)
  endif
  if ((jsl-stencil < G%jsd) .or. (jel+stencil > G%jed)) then
    write(mesg,'("In SIS_continuity, PPM_reconstruction_y called with a ", &
                 & "y-halo that needs to be increased by ",i2,".")') &
                 stencil + max(G%jsd-jsl,jel-G%jed)
    call SIS_error(FATAL,mesg)
  endif

  if (use_2nd) then
    !$OMP parallel do default(shared) private(h_jm1,h_jp1)
    do j=jsl,jel ; do i=isl,iel
      h_jm1 = G%mask2dT(i,j-1) * h_in(i,j-1) + (1.0-G%mask2dT(i,j-1)) * h_in(i,j)
      h_jp1 = G%mask2dT(i,j+1) * h_in(i,j+1) + (1.0-G%mask2dT(i,j+1)) * h_in(i,j)
      h_l(i,j) = 0.5*( h_jm1 + h_in(i,j) )
      h_r(i,j) = 0.5*( h_jp1 + h_in(i,j) )
    enddo ; enddo
  else
    !$OMP parallel do default(shared)  private(dMx,dMn)
    do j=jsl-1,jel+1 ; do i=isl,iel
      if ((G%mask2dT(i,j-1) * G%mask2dT(i,j) * G%mask2dT(i,j+1)) == 0.0) then
        slp(i,j) = 0.0
      else
        ! This uses a simple 2nd order slope.
        slp(i,j) = 0.5 * (h_in(i,j+1) - h_in(i,j-1))
        ! Monotonic constraint, see Eq. B2 in Lin 1994, MWR (132)
        dMx = max(h_in(i,j+1), h_in(i,j-1), h_in(i,j)) - h_in(i,j)
        dMn = h_in(i,j) - min(h_in(i,j+1), h_in(i,j-1), h_in(i,j))
        slp(i,j) = sign(1.,slp(i,j)) * min(abs(slp(i,j)), 2. * min(dMx, dMn))
                ! * (G%mask2dT(i,j-1) * G%mask2dT(i,j) * G%mask2dT(i,j+1))
      endif
    enddo ; enddo
    !$OMP parallel do default(shared) private(h_jm1,h_jp1)
    do j=jsl,jel ; do i=isl,iel
      ! Neighboring values should take into account any boundaries.
      h_jm1 = G%mask2dT(i,j-1) * h_in(i,j-1) + (1.0-G%mask2dT(i,j-1)) * h_in(i,j)
      h_jp1 = G%mask2dT(i,j+1) * h_in(i,j+1) + (1.0-G%mask2dT(i,j+1)) * h_in(i,j)
      ! Left/right values following Eq. B2 in Lin 1994, MWR (132)
      h_l(i,j) = 0.5*( h_jm1 + h_in(i,j) ) + oneSixth*( slp(i,j-1) - slp(i,j) )
      h_r(i,j) = 0.5*( h_jp1 + h_in(i,j) ) + oneSixth*( slp(i,j) - slp(i,j+1) )
    enddo ; enddo
  endif

  if (use_CW84) then
    call PPM_limit_CW84(h_in, h_l, h_r, G, isl, iel, jsl, jel)
  else
    call PPM_limit_pos(h_in, h_l, h_r, h_min, G, isl, iel, jsl, jel)
  endif

  return
end subroutine PPM_reconstruction_y

end module MOM_SIS_continuity
