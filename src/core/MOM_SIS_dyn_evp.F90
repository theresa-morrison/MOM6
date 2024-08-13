!> Update sea-ice dynamics using elastic-viscous-plastic rheology with a C-grid discretization
module MOM_SIS_dyn_evp

! This file is a part of SIS2. See LICENSE.md for the license.

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!                                                                              !
! C-grid SEA ICE DYNAMICS using ELASTIC-VISCOUS-PLASTIC RHEOLOGY adapted from  !
! Hunke and Dukowicz (JPO 1997, H&D hereafter) with some derivation from SIS1  !
! and with guidance from the C-grid implementations of sea-ice in MITgcm as    !
! documented in MITgcm user notes by Martin Losch and in LIM3 by S. Bouillon   !
! et al. (Ocean Modelling, 2009 & 2013). This code initially written by        !
! Robert Hallberg in 2013.                                                     !
!                                                                              !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

use ice_grid,          only : ice_grid_type

use MOM_error_handler, only : MOM_error, FATAL, WARNING, NOTE, MOM_mesg
use MOM_domains,       only : pass_var, pass_vector, CGRID_NE, CORNER, pe_here
use MOM_time_manager,  only : time_type, real_to_time, operator(+), operator(-)
use MOM_unit_scaling,  only : unit_scale_type
use MOM_diag_mediator, only : disable_averaging, post_data, safe_alloc_ptr

!use SIS_diag_mediator, only : post_data, SIS_diag_ctrl
!use SIS_diag_mediator, only : query_SIS_averaging_enabled, enable_SIS_averaging
!use SIS_diag_mediator, only : register_diag_field=>register_SIS_diag_field
!use SIS_debugging,     only : chksum, Bchksum, hchksum, uvchksum
!use SIS_debugging,     only : check_redundant_B, check_redundant_C
!use SIS_restart,       only : register_restart_field, only_read_from_restarts, SIS_restart_CS
!use SIS_restart,       only : query_initialized=>query_inited
!use SIS_framework,     only : safe_alloc
!use SIS_hor_grid,      only : SIS_hor_grid_type
!use SIS_types,         only : ice_state_type
!use SIS2_ice_thm,      only : get_SIS2_thermo_coefs

use MOM_grid, only : ocean_grid_type
!use MOM_forcing_type, only : mech_forcing
!use MOM_forcing_type,  only : SIS_C_EVP_state
use MOM_SIS_C_dyn_CS_type, only : SIS_C_dyn_CS

use MOM_debugging,     only : hchksum, uvchksum

!use combined_ice_ocean_driver, only : direct_copy_from_EVPT, direct_copy_to_EVPT

implicit none ; private

#include <SIS2_memory.h>

public :: EVP_step_loop, SIS_C_EVP_state

contains

!> The control structure with the state that is used and updated in the EVP loop
type, public :: SIS_C_EVP_state
  type(SIS_C_dyn_CS), pointer     :: SIS_C_dyn_CSp => NULL()
  real ::  dt_slow 
  real, allocatable, dimension(:,:) :: ci  !< Sea ice concentration [nondim]
  real, allocatable, dimension(:,:) :: ui    !< Zonal ice velocity [L T-1 ~> m s-1]
  real, allocatable, dimension(:,:) :: vi    !< Meridional ice velocity [L T-1 ~> m s-1]
  real, allocatable, dimension(:,:) :: mice  !< Mass per unit ocean area of sea ice [R Z ~> kg m-2]
  
  real, allocatable, dimension(:,:) :: Cor_u
  real, allocatable, dimension(:,:) :: Cor_v
  
  real, allocatable, dimension(:,:) :: &  
    fxat  !< Zonal air stress on ice [R Z L T-2 ~> Pa]
    fxoc  !<
    fxlf  !<
    fxic  !<
    fxic_d
    fxic_t
    fxic_s
  
  real, allocatable, dimension(:,:) :: &  
    fyat  !< Meridional air stress on ice [R Z L T-2 ~> Pa]
    fyoc  !<
    fylf  !<
    fyic  !<
    fyic_d
    fyic_t
    fyic_s

  real, allocatable, dimension(:,:) :: &
    pres_mice, & ! The ice internal pressure per unit column mass [L2 T-2 ~> N m kg-1].
    diag_val, & ! A temporary diagnostic array.
    del_sh_min_pr     ! When multiplied by pres_mice, this gives the minimum
                ! value of del_sh that is used in the calculation of zeta [T-1 ~> s-1].
                ! This is set based on considerations of numerical stability,
                ! and varies with the grid spacing.  

  real, allocatable, dimension(:,:) :: &
    ui_min_trunc, &  ! The range of v-velocities beyond which the velocities
    ui_max_trunc, &  ! are truncated [L T-1 ~> m s-1], or 0 for land cells
    mi_u, &  ! The total ice and snow mass interpolated to u points [R Z ~> kg m-2].
    f2dt_u, &! The squared effective Coriolis parameter at u-points times a
             ! time step [T-1 ~> s-1].
    PFu, &   ! Zonal hydrostatic pressure driven acceleration [L T-2 ~> m s-2].
    I1_f2dt2_u  ! 1 / ( 1 + f^2 dt^2) at u-points [nondim].
             
  real, allocatable, dimension(:,:) :: &
    vi_min_trunc, &  ! The range of v-velocities beyond which the velocities
    vi_max_trunc, &  ! are truncated [L T-1 ~> m s-1], or 0 for land cells.
    mi_v, &  ! The total ice and snow mass interpolated to v points [R Z ~> kg m-2].
    f2dt_v, &! The squared effective Coriolis parameter at v-points times a
             ! time step [T-1 ~> s-1].
    PFv, &   !  hydrostatic pressure driven acceleration [L T-2 ~> m s-2].
    I1_f2dt2_v  ! 1 / ( 1 + f^2 dt^2) at v-points [nondim].
    
  real, allocatable, dimension(:,:) :: &
    azon, bzon, & !  _zon & _mer are the values of the Coriolis force which
    czon, dzon, & ! are applied to the neighboring values of vi & ui,
    amer, bmer, & ! respectively to get the barotropic inertial rotation,
    cmer, dmer    ! in units of [T-1 ~> s-1].  azon and amer couple the same pair of
                  ! velocities, but with the influence going in opposite
                  ! directions.
                  
  real, allocatable, dimension(:,:) :: &
    mi_ratio_A_q    ! A ratio of the masses interpolated to the faces around a
             ! vorticity point that ranges between (4 mi_min/mi_max) and 1,
             ! divided by the sum of the ocean areas around a point [L-2 ~> m-2].                              
end type SIS_C_EVP_state
! subroutine EVP_step_loop(dt_slow, ci, ui, vi, mice, uo, vo, &
!                     fxat, fyat, fxoc, fyoc, pres_mice, diag_val, del_sh_min_pr, &
!                     ui_min_trunc, ui_max_trunc, vi_min_trunc, vi_max_trunc, &
!                     mi_u, f2dt_u, I1_f2dt2_u, PFu, mi_v, f2dt_v, I1_f2dt2_v, PFv, &
!                     azon, bzon, czon, dzon, amer, bmer, cmer, dmer, &
!                     mi_ratio_A_q,  &
!                     G, US, CS)
! The EVP function from the SIS2 model
subroutine EVP_step_loop(EVPT, G, uo, vo, PFu, PFv, fxoc, fyoc)
! In need from sea-ice: ci, ui, vi, mice, pres_mice, diag_val (?), del_sh_min_pr (?),
!                       ui_min_trunc, ui_max_trunc, vi_min_trunc, vi_max_trunc,
!                       mi_u, f2dt_u, I1_f2dt2_u, mi_v, f2dt_v, I1_f2dt2_v, 
!                       azon, bzon, czon, dzon, amer, bmer, cmer, dmer, mi_ratio_A_q
! Could get somewhere else: fxat, fyat, dt_slow (?)
! need new version from btstep: uo, vo, PFu, PFv (PF ~ sea level) 
! Out: pass to SIS2: ui, vi
! Out: pass to bt_step: fxoc, fyoc
  
  type(SIS_C_EVP_state) , intent(in) :: EVPT
 
  type(ocean_grid_type),  intent(   in) :: G       !< The ocean's grid structure.  
                    
  real, dimension(SZIB_(G),SZJ_( G)), intent(in   ) :: uo    !< Zonal ocean velocity [L T-1 ~> m s-1]
  real, dimension(SZI_( G),SZJB_(G)), intent(in   ) :: vo    !< Meridional ocean velocity [L T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJ_( G)), intent(inout) :: fxoc  !< Zonal ice stress on ocean [R Z L T-2 ~> Pa]
  real, dimension(SZI_( G),SZJB_(G)), intent(inout) :: fyoc  !< Meridional ice stress on ocean [R Z L T-2 ~> Pa]
  
  ! Local Vars                               
  !type(SIS_hor_grid_type)            :: G   !< The horizontal grid type
  type(unit_scale_type)              :: US    !< A structure with unit conversion factors
  type(SIS_C_dyn_CS)                 :: CS    !< The control structure for this module
                     
  real  :: dt_slow !< The amount of time over which the ice
                                                            !! dynamics are to be advanced [T ~> s].

  real, dimension(SZI_(G),SZJ_(G))    :: ci  !< Sea ice concentration [nondim]
  real, dimension(SZIB_(G),SZJ_(G))   :: ui    !< Zonal ice velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G))   :: vi    !< Meridional ice velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G))    :: mice  !< Mass per unit ocean area of sea ice [R Z ~> kg m-2]
  real, dimension(SZIB_(G),SZJ_(G))   :: fxat  !< Zonal air stress on ice [R Z L T-2 ~> Pa]
  real, dimension(SZI_(G),SZJB_(G))   :: fyat  !< Meridional air stress on ice [R Z L T-2 ~> Pa]

  real, dimension(SZI_(G),SZJ_(G))  :: &
    pres_mice, & ! The ice internal pressure per unit column mass [L2 T-2 ~> N m kg-1].
    diag_val, & ! A temporary diagnostic array.
    del_sh_min_pr     ! When multiplied by pres_mice, this gives the minimum
                ! value of del_sh that is used in the calculation of zeta [T-1 ~> s-1].
                ! This is set based on considerations of numerical stability,
                ! and varies with the grid spacing.  

  real, dimension(SZIB_(G),SZJ_(G))  :: &
    ui_min_trunc, &  ! The range of v-velocities beyond which the velocities
    ui_max_trunc, &  ! are truncated [L T-1 ~> m s-1], or 0 for land cells
    mi_u, &  ! The total ice and snow mass interpolated to u points [R Z ~> kg m-2].
    f2dt_u, &! The squared effective Coriolis parameter at u-points times a
             ! time step [T-1 ~> s-1].
    PFu, &   ! Zonal hydrostatic pressure driven acceleration [L T-2 ~> m s-2].
    I1_f2dt2_u  ! 1 / ( 1 + f^2 dt^2) at u-points [nondim].
             
  real, dimension(SZI_(G),SZJB_(G))  :: &
    vi_min_trunc, &  ! The range of v-velocities beyond which the velocities
    vi_max_trunc, &  ! are truncated [L T-1 ~> m s-1], or 0 for land cells.
    mi_v, &  ! The total ice and snow mass interpolated to v points [R Z ~> kg m-2].
    f2dt_v, &! The squared effective Coriolis parameter at v-points times a
             ! time step [T-1 ~> s-1].
    PFv, &   !  hydrostatic pressure driven acceleration [L T-2 ~> m s-2].
    I1_f2dt2_v  ! 1 / ( 1 + f^2 dt^2) at v-points [nondim].
    
  real, dimension(SZIB_(G),SZJ_(G))  :: &
    azon, bzon, & !  _zon & _mer are the values of the Coriolis force which
    czon, dzon, & ! are applied to the neighboring values of vi & ui,
    amer, bmer, & ! respectively to get the barotropic inertial rotation,
    cmer, dmer    ! in units of [T-1 ~> s-1].  azon and amer couple the same pair of
                  ! velocities, but with the influence going in opposite
                  ! directions.
                  
  real, dimension(SZIB_(G),SZJB_(G))  :: &
    mi_ratio_A_q    ! A ratio of the masses interpolated to the faces around a
             ! vorticity point that ranges between (4 mi_min/mi_max) and 1,
             ! divided by the sum of the ocean areas around a point [L-2 ~> m-2].

  real, dimension(SZI_(G),SZJ_(G)) :: &
    sh_Dt, &    ! sh_Dt is the horizontal tension (du/dx - dv/dy) including
                ! all metric terms [T-1 ~> s-1].
    sh_Dd       ! sh_Dd is the flow divergence (du/dx + dv/dy) including all
                ! metric terms [T-1 ~> s-1].
  real, dimension(SZIB_(G),SZJB_(G)) :: &
    sh_Ds       ! sh_Ds is the horizontal shearing strain (du/dy + dv/dx)
                ! including all metric terms [T-1 ~> s-1].
    
  real, dimension(SZI_(G),SZJ_(G)) :: &
    ci_proj, &  ! The projected ice concentration [nondim]. 
    zeta, &     ! The ice bulk viscosity [R Z L2 T-1 ~> Pa m s] (i.e., [N s m-1]).
    del_sh, &   ! The magnitude of the shear rates [T-1 ~> s-1].
    dx2T, dy2T, &   ! dx^2 or dy^2 at T points [L2 ~> m2].
    dx_dyT, dy_dxT  ! dx/dy or dy_dx at T points [nondim].

  real, dimension(SZIB_(G),SZJ_(G)) :: &
    fxic, &   ! Zonal force due to internal stresses [R Z L T-2 ~> Pa].
    fxic_d, & ! Zonal force due to divergence internal stress [R Z L T-2 ~> Pa].
    fxic_t, & ! Zonal force due to tension internal stress [R Z L T-2 ~> Pa].
    fxic_s, & ! Zonal force due to shearing internal stress [R Z L T-2 ~> Pa].
    fxlf, &   ! Zonal landfast ice stress [R Z L T-2 ~> Pa]
    Cor_u, & ! Zonal Coriolis acceleration [L T-2 ~> m s-2]. 
    u_tmp    ! A temporary copy of the old values of ui [L T-1 ~> m s-1].
    
  real, dimension(SZI_(G),SZJB_(G)) :: &
    fyic, &   ! Meridional force due to internal stresses [R Z L T-2 ~> Pa].
    fyic_d, & ! Meridional force due to divergence internal stress [R Z L T-2 ~> Pa].
    fyic_t, & ! Meridional force due to tension internal stress [R Z L T-2 ~> Pa].
    fyic_s, & ! Meridional force due to shearing internal stress [R Z L T-2 ~> Pa].
    fylf, &   ! Meridional landfast ice stress [R Z L T-2 ~> Pa]
    Cor_v     ! Meridional Coriolis acceleration [L T-2 ~> m s-2].

  real, dimension(SZIB_(G),SZJB_(G)) :: &
    dx2B, dy2B, &   ! dx^2 or dy^2 at B points [L2 ~> m2].
    dx_dyB, dy_dxB  ! dx/dy or dy_dx at B points [nondim].
             
  real :: Cor       ! A Coriolis acceleration [L T-2 ~> m s-2].       
  real :: fxic_now  ! Zonal ice internal stress convergence [R Z L T-2 ~> Pa].
  real :: fyic_now  ! Meridional ice internal stress convergence [R Z L T-2 ~> Pa].  
  real :: drag_u, drag_v ! Drag rates with the ocean at u & v points [R Z T-1 ~> kg m-2 s-1].    
  real :: drag_LFu  ! Drag rates to the land for landfast ice at u points [R Z T-1 ~> kg m-2 s-1].
  real :: drag_LFv  ! Drag rates to the land for landfast ice at v points [R Z T-1 ~> kg m-2 s-1].
  real :: drag_max  ! A maximum drag rate allowed in the ocean [R Z T-1 ~> kg m-2 s-1].
  
  real :: v2_at_u     ! The squared v-velocity interpolated to u points [L2 T-2 ~> m2 s-2].
  real :: u2_at_v     ! The squared u-velocity interpolated to v points [L2 T-2 ~> m2 s-2].
  real :: v2_at_u_min ! The squared v-velocity interpolated to u points [L2 T-2 ~> m2 s-2].
  real :: u2_at_v_min ! The squared u-velocity interpolated to v points [L2 T-2 ~> m2 s-2].
  real :: uio_init    ! Ice-ocean velocity differences [L T-1 ~> m s-1]
  real :: vio_init    ! Ice-ocean velocity differences [L T-1 ~> m s-1]
  real :: m_uio_explicit ! Ice-ocean x-velocity differences times the ice mass [R Z L T-1 ~> kg m-1 s-1]
  real :: m_vio_explicit ! Ice-ocean y-velocity differences times the ice mass [R Z L T-1 ~> kg m-1 s-1]
  real :: uio_pred    ! Ice-ocean x-velocity differences [L T-1 ~> m s-1]
  real :: vio_pred    ! Ice-ocean y-velocity differences [L T-1 ~> m s-1]
  real :: I_cdRhoDt   ! The inverse of the product of the drag coefficient, ocean density and
                      ! timestep [L Z-1 R-1 T-1 ~> m3 kg-1 s-1].
  real :: cdRho       ! The ice density times the drag coefficient and rescaling factors [R Z L-1 ~> kg m-3]
  real :: b_vel0      ! The initial difference between the velocity magnitude
                      ! and the absolute value of the u- or v- component, plus
                      ! the ice thickness divided by the time step and the drag
                      ! coefficient [L T-1 ~> m s-1].
  real :: uio_C   ! A u-velocity difference between the ocean and ice [L T-1 ~> m s-1].
  real :: vio_C   ! A v-velocity difference between the ocean and ice [L T-1 ~> m s-1].
  
  real :: I_1pdt_T    ! 1.0 / (1.0 + dt_2Tdamp) [nondim].
  real :: I_1pE2dt_T  ! 1.0 / (1.0 + EC^2 * dt_2Tdamp) [nondim].
  
  real :: EC2     ! EC^2, where EC is the yield curve axis ratio.
  real :: I_EC2   ! 1/EC^2, where EC is the yield curve axis ratio.
  real :: I_EC    ! 1/EC, where EC is the yield curve axis ratio.
  real :: I_2EC   ! 1/(2*EC), where EC is the yield curve axis ratio.
  real, parameter :: H_subroundoff = 1e-30 ! A negligible ice thickness [m].
  
  real :: m_neglect  ! A tiny mass per unit area [R Z ~> kg m-2].
  real :: m_neglect2 ! A tiny mass per unit area squared [R2 Z2 ~> kg2 m-4].
  real :: m_neglect4 ! A tiny mass per unit area to the 4th power [R4 Z4 ~> kg4 m-8]. 
 
  real :: Tdamp   ! The damping timescale of the stress tensor components
                  ! toward their equilibrium solution due to the elastic terms [T ~> s].
  real :: dt      ! The short timestep associated with the EVP dynamics [T ~> s].
  real :: dt_2Tdamp ! The ratio of the timestep to the elastic damping timescale [nondim].
  real :: dt_cumulative ! The elapsed time within this call to EVP dynamics [T ~> s].
  integer :: EVP_steps ! The number of EVP sub-steps that will actually be taken.
  
  type(time_type) :: &
    time_it_start, &  ! The starting time of the iterative steps.
    time_step_end, &  ! The end time of an iterative step.
    time_end_in       ! The end time for diagnostics when this routine started.
  real :: time_int_in ! The diagnostics' time interval when this routine started.
  logical :: do_hifreq_output  ! If true, output occurs every iterative step.
  logical :: do_trunc_its  ! If true, overly large velocities in the iterations are truncated.

  integer :: halo_sh_Ds  ! The halo size that can be used in calculating sh_Ds.
  integer :: i, j, isc, iec, jsc, jec, n
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  ! first: extract the info from the EVPT 
  call direct_copy_from_EVPT(forces%EVPT, CS, dt_slow, G, ci, ui, vi, mice,  &
                        fxat, fyat, pres_mice, diag_val, del_sh_min_pr, &
                        ui_min_trunc, ui_max_trunc, vi_min_trunc, vi_max_trunc, &
                        mi_u, f2dt_u, I1_f2dt2_u, PFu, mi_v, f2dt_v, I1_f2dt2_v, PFv, &
                        azon, bzon, czon, dzon, amer, bmer, cmer, dmer, &
                        mi_ratio_A_q)

  if (CS%dt_Rheo > 0.0) then
    EVP_steps = max(CEILING(dt_slow/CS%dt_Rheo - 0.0001), 1)
  else
    EVP_steps = CS%evp_sub_steps
  endif
  dt = dt_slow/EVP_steps
  
  Tdamp = CS%Tdamp
  if (CS%Tdamp == 0.0) then
    ! Hunke (2001) chooses a specified multiple (0.36) of dt_slow for Tdamp, and shows that
    ! stability requires Tdamp > 2*dt.  Here 0.2 is used instead for greater stability.
    Tdamp = max(0.2*dt_slow, 3.0*dt)
  elseif (CS%Tdamp < 0.0) then
    Tdamp = max(-CS%Tdamp*dt_slow, 3.0*dt)
  endif
  dt_2Tdamp = dt / (2.0 * Tdamp)
  
  
  do J=jsc-1,jec ; do I=isc-1,iec
    dx2B(I,J) = G%dxBu(I,J)*G%dxBu(I,J) ; dy2B(I,J) = G%dyBu(I,J)*G%dyBu(I,J)
  enddo ; enddo
  do J=jsc-2,jec+1 ; do I=isc-2,iec+1
    dx_dyB(I,J) = G%dxBu(I,J)*G%IdyBu(I,J) ; dy_dxB(I,J) = G%dyBu(I,J)*G%IdxBu(I,J)
  enddo ; enddo
  do j=jsc-1,jec+1 ; do i=isc-1,iec+1
    dx2T(i,j) = G%dxT(i,j)*G%dxT(i,j) ; dy2T(i,j) = G%dyT(i,j)*G%dyT(i,j)
    dx_dyT(i,j) = G%dxT(i,j)*G%IdyT(i,j) ; dy_dxT(i,j) = G%dyT(i,j)*G%IdxT(i,j)
  enddo ; enddo  

  halo_sh_Ds = min(isc-G%isd, jsc-G%jsd, 2)

  ! Zero these arrays to accumulate sums.
  fxoc(:,:) = 0.0 ; fyoc(:,:) = 0.0
  fxlf(:,:) = 0.0 ; fylf(:,:) = 0.0
  fxic(:,:) = 0.0 ; fyic(:,:) = 0.0
  Cor_u(:,:) = 0.0 ; Cor_v(:,:) = 0.0
  fxic_d(:,:) = 0.0 ; fyic_d(:,:) = 0.0 
  fxic_t(:,:) = 0.0 ; fyic_t(:,:) = 0.0
  fxic_s(:,:) = 0.0 ; fyic_s(:,:) = 0.0

  drag_max = CS%Rho_ocean * CS%min_ocn_inertial_h / dt_slow
  I_cdRhoDt = 1.0 / (CS%cdw * US%L_to_Z*CS%Rho_ocean * dt)
  do_trunc_its = (CS%CFL_check_its .and. (CS%CFL_trunc > 0.0) .and. (dt_slow > 0.0))

  EC2 = CS%EC**2
  I_EC = 0.0 ; if (CS%EC > 0.0) I_EC = 1.0 / CS%EC
  I_2EC = 0.0 ; if (CS%EC > 0.0) I_2EC = 0.5 / CS%EC
  I_EC2 = 0.0 ; if (EC2 > 0.0) I_EC2 = 1.0 / EC2
  
  do_hifreq_output = .false.
  if ((CS%id_ui_hifreq > 0) .or. (CS%id_vi_hifreq > 0) .or. &
      (CS%id_str_d_hifreq > 0) .or. (CS%id_str_t_hifreq > 0) .or. &
      (CS%id_str_s_hifreq > 0) .or. (CS%id_sh_d_hifreq > 0) .or. &
      (CS%id_sh_t_hifreq > 0) .or. (CS%id_sh_s_hifreq > 0) .or. &
      (CS%id_ci_hifreq > 0) .or. (CS%id_stren_hifreq > 0)) then
    ! do_hifreq_output = query_SISinMOM_averaging_enabled(CS%diag, time_int_in, time_end_in)
    if (do_hifreq_output) &
      time_it_start = time_end_in - real_to_time(US%T_to_s*dt_slow)
  endif
  
  m_neglect = H_subroundoff*US%m_to_Z*CS%Rho_ice
  m_neglect2 = m_neglect**2 ; m_neglect4 = m_neglect**4
  
  dt_cumulative = 0.0

  ! Do the iterative time steps.
  do n=1,EVP_steps

    dt_cumulative = dt_cumulative + dt
    ! If there is a 2-point wide halo and symmetric memory, this is the only
    ! halo update that is needed per iteration.  With a 1-point wide halo and
    ! symmetric memory, an update is also needed for sh_Ds.
    call pass_vector(ui, vi, G%Domain, stagger=CGRID_NE)

    !    Calculate the strain tensor for viscosities and forcing elastic eqn.
    !  The following are the forms of the horizontal tension and horizontal
    !  shearing strain advocated by Smagorinsky (1993) and discussed in
    !  Griffies and Hallberg (MWR, 2000).  Similar forms are used in the sea
    !  ice model of Bouillon et al. (Ocean Modelling, 2009).

    !   The calculation of sh_Ds has the widest halo. The logic below avoids
    ! a halo update when possible.
    !   With a halo of >= 2 this is:  do J=jsc-2,jec+1 ; do I=isc-2,iec+1
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,halo_sh_Ds,sh_Ds,G, &
!$OMP                                  dx_dyB,dy_dxB,ui,vi)
    do J=jsc-halo_sh_Ds,jec+halo_sh_Ds-1 ; do I=isc-halo_sh_Ds,iec+halo_sh_Ds-1
      ! This uses a no-slip boundary condition.
      sh_Ds(I,J) = (2.0-G%mask2dBu(I,J)) * &
          (dx_dyB(I,J)*(ui(I,j+1)*G%IdxCu(I,j+1) - ui(I,j)*G%IdxCu(I,j)) + &
           dy_dxB(I,J)*(vi(i+1,J)*G%IdyCv(i+1,J) - vi(i,J)*G%IdyCv(i,J)))
    enddo ; enddo
    if (halo_sh_Ds < 2) call pass_var(sh_Ds, G%Domain, position=CORNER)
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,sh_Dt,sh_Dd,dy_dxT,dx_dyT,G,ui,vi)
    do j=jsc-1,jec+1 ; do i=isc-1,iec+1
      sh_Dt(i,j) = (dy_dxT(i,j)*(G%IdyCu(I,j) * ui(I,j) - &
                                 G%IdyCu(I-1,j)*ui(I-1,j)) - &
                    dx_dyT(i,j)*(G%IdxCv(i,J) * vi(i,J) - &
                                 G%IdxCv(i,J-1)*vi(i,J-1)))
      sh_Dd(i,j) = (G%IareaT(i,j)*(G%dyCu(I,j) * ui(I,j) - &
                                   G%dyCu(I-1,j)*ui(I-1,j)) + &
                    G%IareaT(i,j)*(G%dxCv(i,J) * vi(i,J) - &
                                   G%dxCv(i,J-1)*vi(i,J-1)))
    enddo ; enddo

   if (CS%project_ci) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ci_proj,ci,dt_cumulative, &
!$OMP                                  sh_Dd,pres_mice,CS)
     do j=jsc-1,jec+1 ; do i=isc-1,iec+1
       ! Estimate future ice concentrations from the approximate expression
       !   d ci / dt = - ci * sh_Dt
       ! The choice to base this on the final velocity, the initial concentration
       ! and the elapsed time is because it is that final velocity that will drive
       ! ice convergence.
       ci_proj(i,j) = ci(i,j) * exp(-dt_cumulative*sh_Dd(i,j))
       ! Recompute pres_mice.
       pres_mice(i,j) = CS%p0_rho*exp(-CS%c0*max(1.0-ci_proj(i,j),0.0))
     enddo ; enddo
   endif

   ! calculate viscosities - how often should we do this ?
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,del_sh,zeta,sh_Dd,sh_Dt, &
!$OMP                                  I_EC2,sh_Ds,pres_mice,mice,del_sh_min_pr)
    do j=jsc-1,jec+1 ; do i=isc-1,iec+1
      ! Averaging the squared shearing strain is larger than squaring
      ! the averaged strain.  I don't know what is better. -RWH
      del_sh(i,j) = sqrt(sh_Dd(i,j)**2 + I_EC2 * (sh_Dt(i,j)**2 + &
                   (0.25 * ((sh_Ds(I-1,J-1) + sh_Ds(I,J)) + &
                            (sh_Ds(I-1,J) + sh_Ds(I,J-1))))**2 ) ) ! H&D eqn 9

      if (max(del_sh(i,j), del_sh_min_pr(i,j)*pres_mice(i,j)) /= 0.) then
        zeta(i,j) = 0.5*pres_mice(i,j)*mice(i,j) / &
           max(del_sh(i,j), del_sh_min_pr(i,j)*pres_mice(i,j))
      else
        zeta(i,j) = 0.
      endif
    enddo ; enddo

    ! Step the stress component equations semi-implicitly.
    I_1pdt_T = 1.0 / (1.0 + dt_2Tdamp)
    I_1pE2dt_T = 1.0 / (1.0 + EC2*dt_2Tdamp)
    if (CS%weak_low_shear) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,CS,I_1pdt_T,dt_2Tdamp,zeta, &
!$OMP                                  sh_Dd,del_sh,I_EC2,sh_Dt)
      do j=jsc-1,jec+1 ; do i=isc-1,iec+1
        ! This expression uses that Pres=2*del_sh*zeta with an elliptic yield curve.
        CS%str_d(i,j) = I_1pdt_T * ( CS%str_d(i,j) + dt_2Tdamp * &
                    ( zeta(i,j) * (sh_Dd(i,j) - del_sh(i,j)) ) )
        CS%str_t(i,j) = I_1pdt_T * ( CS%str_t(i,j) + (I_EC2 * dt_2Tdamp) * &
                    ( zeta(i,j) * sh_Dt(i,j) ) )
      enddo ; enddo
    else
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,CS,I_1pdt_T,dt_2Tdamp,zeta, &
!$OMP                                  sh_Dd,I_EC2,sh_Dt,pres_mice,mice)
      do j=jsc-1,jec+1 ; do i=isc-1,iec+1
        ! This expression uses that Pres=2*del_sh*zeta with an elliptic yield curve.
        CS%str_d(i,j) = I_1pdt_T * ( CS%str_d(i,j) + dt_2Tdamp * &
                    ( zeta(i,j) * sh_Dd(i,j) - 0.5*pres_mice(i,j)*mice(i,j) ) )
        CS%str_t(i,j) = I_1pdt_T * ( CS%str_t(i,j) + (I_EC2 * dt_2Tdamp) * &
                    ( zeta(i,j) * sh_Dt(i,j) ) )
      enddo ; enddo
    endif
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,CS,I_1pdt_T,I_EC2,dt_2Tdamp, &
!$OMP                                  G,zeta,mi_ratio_A_q,sh_Ds)
    do J=jsc-1,jec ; do I=isc-1,iec
      ! zeta is already set to 0 over land.
      CS%str_s(I,J) = I_1pdt_T * ( CS%str_s(I,J) + (I_EC2 * dt_2Tdamp) * &
                  ( ((G%areaT(i,j)*zeta(i,j) + G%areaT(i+1,j+1)*zeta(i+1,j+1)) + &
                     (G%areaT(i+1,j)*zeta(i+1,j) + G%areaT(i,j+1)*zeta(i,j+1))) * &
                   mi_ratio_A_q(I,J) * sh_Ds(I,J) ) )
    enddo ; enddo

    if (CS%str_underflow > 0.0) then
      !$OMP parallel do default(shared)
      do j=jsc-1,jec+1 ; do i=isc-1,iec+1
        if (abs(CS%str_d(i,j)) < CS%str_underflow) CS%str_d(i,j) = 0.0
        if (abs(CS%str_t(i,j)) < CS%str_underflow) CS%str_t(i,j) = 0.0
      enddo ; enddo
      !$OMP parallel do default(shared)
      do J=jsc-1,jec ; do I=isc-1,iec
        if (abs(CS%str_s(I,J)) < CS%str_underflow) CS%str_s(I,J) = 0.0
      enddo ; enddo
    endif

    cdRho = CS%cdw * US%L_to_Z*CS%Rho_ocean
    ! Save the current values of u for later use in updating v.
    do I=isc-1,iec
      u_tmp(I,jsc-1) = ui(I,jsc-1) ; u_tmp(I,jec+1) = ui(I,jec+1) ;
    enddo
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,u_tmp,ui,vi,azon,bzon,czon,dzon, &
!$OMP                                  G,CS,dy2T,dx2B,vo,uo,Cor_u,f2dt_u,I1_f2dt2_u,    &
!$OMP                                  mi_u,dt,PFu,fxat,I_cdRhoDt,cdRho,m_neglect,fxoc, &
!$OMP                                  fxlf,fxic,fxic_d,fxic_t,fxic_s,do_trunc_its,drag_max) &
!$OMP                          private(Cor,fxic_now,v2_at_u,v2_at_u_min,uio_init,drag_u,drag_LFu,b_vel0, &
!$OMP                                  m_uio_explicit,uio_pred,uio_C)
    do j=jsc,jec ; do I=isc-1,iec
      ! Save the current values of u for later use in updating v.
      u_tmp(I,j) = ui(I,j)

      Cor = ((azon(I,j) * vi(i+1,J) + czon(I,j) * vi(i,J-1)) + &
             (bzon(I,j) * vi(i,J) + dzon(I,j) * vi(i+1,J-1))) ! - Cor_ref_u(I,j)
      !  Evaluate 1/m x.Div(m strain).  This expressions include all metric terms
      !  for an orthogonal grid.  The str_d term integrates out to no curl, while
      !  str_s & str_t terms impose no divergence and do not act on solid body rotation.
      fxic_now = G%IdxCu(I,j) * (CS%str_d(i+1,j) - CS%str_d(i,j)) + &
            (G%IdyCu(I,j)*(dy2T(i+1,j)*CS%str_t(i+1,j) - &
                           dy2T(i,j)  *CS%str_t(i,j)) + &
             G%IdxCu(I,j)*(dx2B(I,J)  *CS%str_s(I,J) - &
                           dx2B(I,J-1)*CS%str_s(I,J-1)) ) * G%IareaCu(I,j)
      v2_at_u =  CS%drag_bg_vel2 + 0.25 * &
                     (((vi(i,J)-vo(i,J))**2 + (vi(i+1,J-1)-vo(i+1,J-1))**2) + &
                      ((vi(i+1,J)-vo(i+1,J))**2 + (vi(i,J-1)-vo(i,J-1))**2))
      if (CS%lemieux_landfast .or. CS%itd_landfast) &
               v2_at_u_min = min(abs(vi(I,j)), abs(vi(i+1,J-1)), &
                                 abs(vi(i+1,J)), abs(vi(i,J-1)))**2

      uio_init = (ui(I,j)-uo(I,j))

      ! Determine the Coriolis acceleration and sum for averages...
      Cor_u(I,j) = Cor_u(I,j) + (Cor - f2dt_u(I,j) * ui(I,j)) * I1_f2dt2_u(I,j)

      if (CS%project_drag_vel) then
      ! Project the new u-velocity using a quasi-analytic implicit treatment for
      ! drag, but explicit treatments for everything else, to estimate the drag
      ! coefficient, then take the larger of the two estimates of
      ! the ice-ocean drag.
        drag_u = 0.0
        if (G%mask2dCu(I,j) > 0.0) then
          m_uio_explicit = uio_init*mi_u(I,j) + dt * &
               ((Cor + PFu(I,j))*mi_u(I,j) + (fxic_now + fxat(I,j)))
          b_vel0 = mi_u(I,j) * I_cdRhoDt + &
                   ( sqrt(uio_init**2 + v2_at_u) - abs(uio_init) )
          if (b_vel0**2 > 1e8*I_cdRhoDt*abs(m_uio_explicit)) then
            uio_pred = m_uio_explicit * I_cdRhoDt / b_vel0
          else
            uio_pred = 0.5 * (sqrt(b_vel0**2 + 4.0*I_cdRhoDt*abs(m_uio_explicit)) - b_vel0)
          endif
          drag_u = cdRho * sqrt(max(uio_init**2, uio_pred**2) + v2_at_u )
        endif
      else
        drag_u = cdRho * sqrt(uio_init**2 + v2_at_u )
      endif
      if (drag_max>0.) drag_u = min( drag_u, drag_max )
      drag_LFu = 0.0
      if (CS%lemieux_landfast .or. CS%itd_landfast) then
        drag_LFu = CS%Tb_u(I,j) / (sqrt(ui(I,j)**2 + v2_at_u_min ) + CS%lemieux_u0)
      endif

      !   This is a quasi-implicit timestep of Coriolis, followed by an explicit
      ! update of the other terms and an implicit bottom drag calculation.
      uio_C =  G%mask2dCu(I,j) * ( mi_u(I,j) * &
               ((ui(I,j) + dt * Cor) * I1_f2dt2_u(I,j) - uo(I,j)) + &
                dt * ((mi_u(I,j) * PFu(I,j) + (fxic_now + fxat(I,j))) - drag_LFu*uo(I,j)) ) / &
               (mi_u(I,j) + m_neglect + dt * (drag_u + drag_LFu))

      ui(I,j) = (uio_C + uo(I,j)) * G%mask2dCu(I,j)

      ! Note that fxoc is the stress felt by the ocean.
      fxoc(I,j) = fxoc(I,j) + drag_u*uio_C

      ! Here fxlf is the stress felt by the landfast ice.
      fxlf(I,j) = fxlf(I,j) - drag_LFu*ui(I,j)

      ! sum accelerations to take averages.
      fxic(I,j) = fxic(I,j) + fxic_now

      if (CS%id_fix_d>0) fxic_d(I,j) = fxic_d(I,j) + G%mask2dCu(I,j) * &
                 G%IdxCu(I,j) * (CS%str_d(i+1,j) - CS%str_d(i,j))
      if (CS%id_fix_t>0) fxic_t(I,j) = fxic_t(I,j) + G%mask2dCu(I,j) * &
                  G%IdyCu(I,j)*(dy2T(i+1,j)* CS%str_t(i+1,j) - &
                                dy2T(i,j)  * CS%str_t(i,j) ) * G%IareaCu(I,j)
      if (CS%id_fix_s>0) fxic_s(I,j) = fxic_s(I,j) + G%mask2dCu(I,j) * &
                  G%IdxCu(I,j)*(dx2B(I,J)  *CS%str_s(I,J) - &
                                dx2B(I,J-1)*CS%str_s(I,J-1)) * G%IareaCu(I,j)

    enddo ; enddo

!$OMP parallel do default(none) shared(isc,iec,jsc,jec,amer,bmer,cmer,dmer,u_tmp,G,CS, &
!$OMP                                  dx2T,dy2B,uo,vo,vi,Cor_v,f2dt_v,I1_f2dt2_v,mi_v, &
!$OMP                                  dt,PFv,fyat,I_cdRhoDt,cdRho,m_neglect,fyoc,fyic, &
!$OMP                                  fylf,fyic_d,fyic_t,fyic_s,do_trunc_its,vi_min_trunc,  &
!$OMP                                  vi_max_trunc,drag_max) &
!$OMP                          private(Cor,fyic_now,u2_at_v,vio_init,drag_v,drag_LFv,u2_at_v_min, &
!$OMP                                  m_vio_explicit,b_vel0,vio_pred,vio_C)
    do J=jsc-1,jec ; do i=isc,iec
      Cor = -1.0*((amer(I-1,j) * u_tmp(I-1,j) + cmer(I,j+1) * u_tmp(I,j+1)) + &
                  (bmer(I,j) * u_tmp(I,j) + dmer(I-1,j+1) * u_tmp(I-1,j+1)))
      !  Evaluate 1/m y.Div(m strain).  This expressions include all metric terms
      !  for an orthogonal grid.  The str_d term integrates out to no curl, while
      !  str_s & str_t terms impose no divergence and do not act on solid body rotation.
      fyic_now = G%IdyCv(i,J) * (CS%str_d(i,j+1)-CS%str_d(i,j)) + &
            (-G%IdxCv(i,J)*(dx2T(i,j+1)*CS%str_t(i,j+1) - &
                            dx2T(i,j)  *CS%str_t(i,j)) + &
              G%IdyCv(i,J)*(dy2B(I,J)  *CS%str_s(I,J) - &
                            dy2B(I-1,J)*CS%str_s(I-1,J)) )*G%IareaCv(i,J)
      u2_at_v = CS%drag_bg_vel2 + 0.25 * &
                (((u_tmp(I,j)-uo(I,j))**2 + (u_tmp(I-1,j+1)-uo(I-1,j+1))**2) + &
                 ((u_tmp(I,j+1)-uo(I,j+1))**2 + (u_tmp(I-1,j)-uo(I-1,j))**2))
      if (CS%lemieux_landfast .or. CS%itd_landfast) &
                u2_at_v_min = min(abs(u_tmp(i,J)), abs(u_tmp(I-1,j+1)), &
                                  abs(u_tmp(I,j+1)), abs(u_tmp(I-1,j)))**2

      vio_init = (vi(i,J)-vo(i,J))

      ! Determine the Coriolis acceleration and sum for averages...
      Cor_v(I,J) = Cor_v(I,J) + (Cor - f2dt_v(i,J) * vi(i,J)) * I1_f2dt2_v(i,J)

      if (CS%project_drag_vel) then
      ! Project the new v-velocity using a quasi-analytic implicit treatment for
      ! drag, but explicit treatments for everything else, to estimate the drag
      ! coefficient, then take the larger of the two estimates of
      ! the ice-ocean drag.

        drag_v = 0.0
        if (G%mask2dCv(i,J) > 0.0) then
          m_vio_explicit = vio_init*mi_v(i,J) + dt * &
               ((Cor + PFv(i,J))*mi_v(i,J) + (fyic_now + fyat(i,J)))
          b_vel0 = mi_v(i,J) * I_cdRhoDt + (sqrt(vio_init**2 + u2_at_v) - abs(vio_init))
          if (b_vel0**2 > 1e8*I_cdRhoDt*abs(m_vio_explicit)) then
            vio_pred = m_vio_explicit * I_cdRhoDt / b_vel0
          else
            vio_pred = 0.5 * (sqrt(b_vel0**2 + 4.0*I_cdRhoDt*abs(m_vio_explicit)) - b_vel0)
          endif
          drag_v = cdRho * sqrt(max(vio_init**2, vio_pred**2) + u2_at_v )
        endif
      else
        drag_v = cdRho * sqrt(vio_init**2 + u2_at_v )
      endif
      if (drag_max>0.) drag_v = min( drag_v, drag_max )
      drag_LFv = 0.0
      if (CS%lemieux_landfast .or. CS%itd_landfast) then
        drag_LFv = CS%Tb_v(i,J) / (sqrt(vi(i,J)**2 + u2_at_v_min ) + CS%lemieux_u0)
      endif

      !   This is a quasi-implicit timestep of Coriolis, followed by an explicit
      ! update of the other terms and an implicit bottom drag calculation.
      vio_C =  G%mask2dCv(i,J) * ( mi_v(i,J) * &
               ((vi(i,J) + dt * Cor) * I1_f2dt2_v(i,J) - vo(i,J)) + &
                dt * ((mi_v(i,J) * PFv(i,J) + (fyic_now + fyat(i,J))) - drag_LFv*vo(i,J)) ) / &
               (mi_v(i,J) + m_neglect + dt * (drag_v + drag_LFv))

      vi(i,J) = (vio_C + vo(i,J)) * G%mask2dCv(i,J)

      ! Note that fyoc is the stress felt by the ocean.
      fyoc(i,J) = fyoc(i,J) + drag_v*vio_C

      ! Here fylf is the stress felt by the landfast ice.
      fylf(I,j) = fylf(I,j) - drag_LFv*vi(I,j)

      ! sum accelerations to take averages.
      fyic(i,J) = fyic(i,J) + fyic_now

      if (CS%id_fiy_d>0) fyic_d(i,J) = fyic_d(i,J) + G%mask2dCv(i,J) * &
                 G%IdyCv(i,J) * (CS%str_d(i,j+1)-CS%str_d(i,j))
      if (CS%id_fiy_t>0) fyic_t(i,J) = fyic_t(i,J) + G%mask2dCv(i,J) * &
                 (G%IdxCv(i,J)*(dx2T(i,j+1)*(-CS%str_t(i,j+1)) - &
                                dx2T(i,j)  *(-CS%str_t(i,j))) ) * G%IareaCv(i,J)
      if (CS%id_fiy_s>0) fyic_s(i,J) = fyic_s(i,J) + G%mask2dCv(i,J) * &
                 (G%IdyCv(i,J)*(dy2B(I,J)  *CS%str_s(I,J) - &
                                dy2B(I-1,J)*CS%str_s(I-1,J)) ) * G%IareaCv(i,J)

    enddo ; enddo

    ! Apply appropriate limits on the magnitude of the velocies, both to handle
    ! underflow and to keep failing runs going so that they can be diagnosed.
    if (do_trunc_its .or. (CS%vel_underflow > 0.0)) then
      !$OMP parallel do default(shared)
      do j=jsc,jec ; do I=isc-1,iec
        if (abs(ui(I,j)) < CS%vel_underflow) ui(I,j) = 0.0
        if (do_trunc_its) then
          if (ui(I,j) < ui_min_trunc(I,j)) then
            ui(I,j) = ui_min_trunc(I,j)
          elseif (ui(I,j) > ui_max_trunc(I,j)) then
            ui(I,j) = ui_max_trunc(I,j)
          endif
        endif
      enddo ; enddo

      !$OMP parallel do default(shared)
      do J=jsc-1,jec ; do i=isc,iec
        if (abs(vi(i,J)) < CS%vel_underflow) vi(i,J) = 0.0
        if (do_trunc_its) then
          if (vi(i,J) < vi_min_trunc(i,J)) then
            vi(i,J) = vi_min_trunc(i,J)
          elseif (vi(i,J) > vi_max_trunc(i,J)) then
            vi(i,J) = vi_max_trunc(i,J)
          endif
        endif
      enddo ; enddo
    endif

    if (do_hifreq_output) then
      time_step_end = time_it_start + real_to_time(n*US%T_to_s*dt)
      !call enable_SIS_averaging(US%T_to_s*dt, time_step_end, CS%diag)
      !if (CS%id_ui_hifreq > 0) call post_data(CS%id_ui_hifreq, ui, CS%diag)
      !if (CS%id_vi_hifreq > 0) call post_data(CS%id_vi_hifreq, vi, CS%diag)
      !if (CS%id_str_d_hifreq > 0) call post_data(CS%id_str_d_hifreq, CS%str_d, CS%diag)
      !if (CS%id_str_t_hifreq > 0) call post_data(CS%id_str_t_hifreq, CS%str_t, CS%diag)
      !if (CS%id_str_s_hifreq > 0) call post_data(CS%id_str_s_hifreq, CS%str_s, CS%diag)
      !if (CS%id_sh_d_hifreq > 0) call post_data(CS%id_sh_d_hifreq, sh_Dd, CS%diag)
      !if (CS%id_sh_t_hifreq > 0) call post_data(CS%id_sh_t_hifreq, sh_Dt, CS%diag)
      !if (CS%id_sh_s_hifreq > 0) call post_data(CS%id_sh_s_hifreq, sh_Ds, CS%diag)
      !if (CS%id_sigi_hifreq>0) then
      !  call find_sigI(mice, ci_proj, CS%str_d, diag_val, G, US, CS)
      !  call post_data(CS%id_sigi_hifreq, diag_val, CS%diag)
      !endif
      !if (CS%id_sigii_hifreq>0) then
      !  call find_sigII(mice, ci_proj, CS%str_t, CS%str_s, diag_val, G, US, CS)
      !  call post_data(CS%id_sigii_hifreq, diag_val, CS%diag)
      !endif
      !if (CS%id_ci_hifreq>0) call post_data(CS%id_ci_hifreq, ci_proj, CS%diag)
      !if (CS%id_stren_hifreq>0) then
      !  do j=jsc,jec ; do i=isc,iec
      !    diag_val(i,j) = pres_mice(i,j)*mice(i,j)
      !  enddo ; enddo
      !  call post_data(CS%id_stren_hifreq, diag_val, CS%diag)
      !endif
    endif

    if (CS%debug_EVP .and. CS%debug) then
      call hchksum(CS%str_d, "str_d in SIS_C_dynamics", G%HI, haloshift=1, scale=US%RZ_to_kg_m2*US%L_T_to_m_s**2)
      call hchksum(CS%str_t, "str_t in SIS_C_dynamics", G%HI, haloshift=1, scale=US%RZ_to_kg_m2*US%L_T_to_m_s**2)
      !call Bchksum(CS%str_s, "str_s in SIS_C_dynamics", G%HI, &
      !             haloshift=0, symmetric=.true., scale=US%RZ_to_kg_m2*US%L_T_to_m_s**2)
    endif
    if (CS%debug_EVP .and. (CS%debug .or. CS%debug_redundant)) then
      !call uvchksum("f[xy]ic in SIS_C_dynamics", fxic, fyic, G, scale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
      !call uvchksum("f[xy]oc in SIS_C_dynamics", fxoc, fyoc, G, scale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
      !call uvchksum("f[xy]lf in SIS_C_dynamics", fxlf, fylf, G, scale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
      !call uvchksum("Cor_[uv] in SIS_C_dynamics", Cor_u, Cor_v, G, scale=US%L_T_to_m_s*US%s_to_T)
      !call uvchksum("[uv]i in SIS_C_dynamics", ui, vi, G, scale=US%L_T_to_m_s)
    endif

  enddo ! l=1,EVP_steps

end subroutine EVP_step_loop


subroutine direct_copy_from_EVPT(EVPT, CS, dt_slow, G, ci, ui, vi, mice,  &
                        fxat, fyat, pres_mice, diag_val, del_sh_min_pr, &
                        ui_min_trunc, ui_max_trunc, vi_min_trunc, vi_max_trunc, &
                        mi_u, f2dt_u, I1_f2dt2_u, PFu, mi_v, f2dt_v, I1_f2dt2_v, PFv, &
                        azon, bzon, czon, dzon, amer, bmer, cmer, dmer, &
                        mi_ratio_A_q)

  type(SIS_C_EVP_state), intent(in) :: EVPT 
  
  type(SIS_C_dyn_CS),   intent(out)     :: CS
  !type(SIS_hor_grid_type),     intent(inout)    :: G
  type(ocean_grid_type),     intent(in)    :: G
  real, intent(out) :: dt_slow
  
  real, dimension(SZI_(G),SZJ_(G)),  intent(out  ) :: ci  !< Sea ice concentration [nondim]
  real, dimension(SZI_(G),SZJ_(G)),  intent(out  ) :: mice  !< Mass per unit ocean area of sea ice [R Z ~> kg m-2]
  real, dimension(SZIB_(G),SZJ_(G)),  intent(out  ) :: ui    !< Zonal ice velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G)),  intent(out  ) :: vi    !< Meridional ice velocity [L T-1 ~> m s-1]
  
  real, dimension(SZIB_(G),SZJ_(G)),  intent(out  ) :: fxat  !< Zonal air stress on ice [R Z L T-2 ~> Pa]
  real, dimension(SZI_(G),SZJB_(G)),  intent(out  ) :: fyat  !< Meridional air stress on ice [R Z L T-2 ~> Pa]

  real, dimension(SZI_(G),SZJ_(G)),  intent(out  ) :: &
    pres_mice, & ! The ice internal pressure per unit column mass [L2 T-2 ~> N m kg-1].
    diag_val, & ! A temporary diagnostic array.
    del_sh_min_pr     ! When multiplied by pres_mice, this gives the minimum
                ! value of del_sh that is used outthe calculation of zeta [T-1 ~> s-1].
                ! This is set based on considerations of numerical stability,
                ! and varies with the grid spacing.  

  real, dimension(SZIB_(G),SZJ_(G)),  intent(out  )  :: &
    ui_min_trunc, &  ! The range of v-velocities beyond which the velocities
    ui_max_trunc, &  ! are truncated [L T-1 ~> m s-1], or 0 for land cells
    mi_u, &  ! The total ice and snow mass interpolated to u points [R Z ~> kg m-2].
    f2dt_u, &! The squared effective Coriolis parameter at u-points times a
             ! time step [T-1 ~> s-1].
    PFu, &   ! Zonal hydrostatic pressure driven acceleration [L T-2 ~> m s-2].
    I1_f2dt2_u  ! 1 / ( 1 + f^2 dt^2) at u-points [nondim].

  real, dimension(SZI_(G),SZJB_(G)),  intent(out  )  :: &
    vi_min_trunc, &  ! The range of v-velocities beyond which the velocities
    vi_max_trunc, &  ! are truncated [L T-1 ~> m s-1], or 0 for land cells.
    mi_v, &  ! The total ice and snow mass interpolated to v points [R Z ~> kg m-2].
    f2dt_v, &! The squared effective Coriolis parameter at v-points times a
             ! time step [T-1 ~> s-1].
    PFv, &   !  hydrostatic pressure driven acceleration [L T-2 ~> m s-2].
    I1_f2dt2_v  ! 1 / ( 1 + f^2 dt^2) at v-points [nondim].

  real, dimension(SZIB_(G),SZJ_(G)),  intent(out  )  :: &
    azon, bzon, & !  _zon & _mer are the values of the Coriolis force which
    czon, dzon, & ! are applied to the neighboring values of vi & ui,
    amer, bmer, & ! respectively to get the barotropic inertial rotation,
    cmer, dmer    ! outunits of [T-1 ~> s-1].  azon and amer couple the same pair of
                  ! velocities, but with the influence going outopposite
                  ! directions.
                  
  real, dimension(SZIB_(G),SZJB_(G)),  intent(out  )  :: &
    mi_ratio_A_q    ! A ratio of the masses interpolated to the faces around a
             ! vorticity point that ranges between (4 mi_min/mi_max) and 1,
             ! divided by the sum of the ocean areas around a point [L-2 ~> m-2].  

  CS = EVPT%SIS_C_dyn_CSp  

  dt_slow = EVPT%dt_slow

  ci(:,:)   = EVPT%ci(:,:)
  mice(:,:) = EVPT%mice(:,:)
  ui(:,:)   = EVPT%ui(:,:)
  vi(:,:)   = EVPT%vi(:,:)
  fxat(:,:) = EVPT%fxat(:,:)
  fyat(:,:) = EVPT%fyat(:,:)

  pres_mice(:,:)    = EVPT%pres_mice(:,:)
  diag_val(:,:)     = EVPT%diag_val(:,:)
  del_sh_min_pr(:,:)= EVPT%del_sh_min_pr(:,:)

  ui_min_trunc(:,:) = EVPT%ui_min_trunc(:,:)
  ui_max_trunc(:,:) = EVPT%ui_max_trunc(:,:)
  mi_u(:,:)         = EVPT%mi_u(:,:)
  f2dt_u(:,:)       = EVPT%f2dt_u(:,:)
  PFu(:,:)          = EVPT%PFu(:,:)
  I1_f2dt2_u(:,:)   = EVPT%I1_f2dt2_u(:,:)
  
  vi_min_trunc(:,:) = EVPT%vi_min_trunc(:,:)
  vi_max_trunc(:,:) = EVPT%vi_max_trunc(:,:)
  mi_v(:,:)         = EVPT%mi_v(:,:)
  f2dt_v(:,:)       = EVPT%f2dt_v(:,:)
  PFv(:,:)          = EVPT%PFv(:,:)
  I1_f2dt2_v(:,:)   = EVPT%I1_f2dt2_v(:,:)

  azon(:,:) = EVPT%azon(:,:)
  bzon(:,:) = EVPT%bzon(:,:)
  czon(:,:) = EVPT%czon(:,:)
  dzon(:,:) = EVPT%dzon(:,:)

  amer(:,:) = EVPT%amer(:,:)
  bmer(:,:) = EVPT%bmer(:,:)
  cmer(:,:) = EVPT%cmer(:,:)
  dmer(:,:) = EVPT%dmer(:,:)

  mi_ratio_A_q(:,:) = EVPT%mi_ratio_A_q(:,:)

end subroutine direct_copy_from_EVPT


!!> Indicate whether averaging diagnostics is currently enabled
!logical function query_SISinMOM_averaging_enabled(diag_cs, time_int, time_end)
!  type(SIS_diag_ctrl),           intent(in)  :: diag_cs !< A structure that is used to regulate diagnostic output
!  real,            optional, intent(out) :: time_int !< The current setting of diag_cs%time_int [s].
!  type(time_type), optional, intent(out) :: time_end !< The current setting of diag_cs%time_end.
!
!  if (present(time_int)) time_int = diag_cs%time_int
!  if (present(time_end)) time_end = diag_cs%time_end
!  query_SIS_averaging_enabled = diag_cs%ave_enabled
!end function query_SISinMOM_averaging_enabled





end module MOM_SIS_dyn_evp
