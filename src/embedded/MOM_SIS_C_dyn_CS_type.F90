!> Update sea-ice dynamics using elastic-viscous-plastic rheology with a C-grid discretization
module MOM_SIS_C_dyn_CS_type

use MOM_time_manager,  only : time_type
use SIS_diag_mediator, only : SIS_diag_ctrl
implicit none ; private

#include <SIS2_memory.h>

!> The control structure with parameters regulating C-grid ice dynamics
type, public :: SIS_C_dyn_CS ; ! private
  real, allocatable, dimension(:,:) :: &
    str_t, &  !< The tension stress tensor component [R Z L2 T-2 ~> Pa m].
    str_d, &  !< The divergence stress tensor component [R Z L2 T-2 ~> Pa m].
    str_s     !< The shearing stress tensor component (cross term) [R Z L2 T-2 ~> Pa m].

  ! parameters for calculating water drag and internal ice stresses
  real :: p0                  !< Pressure constant in the Hibler rheology [R L2 T-2 ~> Pa]
  real :: p0_rho              !< The pressure constant divided by ice density [L2 T-2 ~> N m kg-1].
  real :: c0                  !< another pressure constant, c* in Hunke & Dukowicz 1997 [nondim]
  real :: cdw                 !< The drag coefficient between the sea ice and water. [nondim]
  real :: EC = 2.0            !< yield curve axis ratio [nondim]
  real :: Rho_ocean           !< The nominal density of sea water [R ~> kg m-3].
  real :: Rho_ice             !< The nominal density of sea ice [R ~> kg m-3].
  real :: drag_bg_vel2 = 0.0  !< A background (subgridscale) velocity for drag with the ocean
                              !< squared [L2 T-2 ~> m2 s-2].  This is always 0 for now.
  real :: min_ocn_inertial_h  !< A minimum ocean thickness used to limit the viscous coupling
                              !! rate implied for the ocean by the ice-ocean stress [Z ~> m].
  real :: Tdamp               !< The damping timescale of the stress tensor components toward their
                              !! equilibrium solution due to the elastic terms [T ~> s] or [nondim].
  real :: del_sh_min_scale    !< A scaling factor for the minimum permitted value of minimum
                              !! shears used in the denominator of the stress equations [nondim].
                              !  I suspect that this needs to be greater than 1. -RWH
  real    :: vel_underflow    !< Velocity components smaller than vel_underflow
                              !! are set to 0 [L T-1 ~> m s-1].
  real    :: str_underflow    !< Stress tensor components smaller than str_underflow
                              !! are set to 0 [R Z L2 T-2 ~> Pa m].
  real    :: CFL_trunc        !< Velocity components will be truncated when they are large enough
                              !! that the corresponding CFL number exceeds this value [nondim].
  logical :: CFL_check_its    !< If true, check the CFL number for every iteration
                              !! of the rheology solver; otherwise only check the
                              !! final velocities that are used for transport.
  logical :: do_embedded      !< If true,
  logical :: embedded_setup   !< If true,
  logical :: embedded_finish  !< If true,
  logical :: debug            !< If true, write verbose checksums for debugging purposes.
  logical :: debug_EVP        !< If true, write out verbose debugging data for each of
                              !! the steps within the EVP solver.
  logical :: debug_redundant  !< If true, debug redundant points.
  logical :: project_drag_vel !< If true, project forward the ice velocity used in the drag
                              !! calculation to avoid an instability that can occur when an finite
                              !! stress is applied to thin ice moving with the velocity of the ocean.
  logical :: project_ci       !< If true, project the ice concentration and related ice strength
                              !! changes due to the convergent or divergent ice flow.
  logical :: weak_coast_stress = .false. !< If true, do not use land masks in determining the area
                              !! for stress convergence, which acts to weaken the stress-driven
                              !! acceleration in coastal points.
  logical :: weak_low_shear = .false. !< If true, the divergent stresses go toward 0 in the C-grid
                              !! dynamics when the shear magnitudes are very weak.
                              !! Otherwise they go to -P_ice.  This setting is temporary.
  integer :: evp_sub_steps    !< The number of iterations in the EVP dynamics
                              !! for each slow time step.
  real    :: dt_Rheo          !< The maximum sub-cycling time step for the EVP dynamics [T ~> s].
  type(time_type), pointer :: Time => NULL() !< A pointer to the ice model's clock.
  type(SIS_diag_ctrl), pointer :: diag => NULL() !< A structure that is used to regulate the
                              !! timing of diagnostic output.
  integer, pointer :: ntrunc => NULL() !< The number of times the velocity has been truncated
                              !! since the last call to write_ice_statistics.
  character(len = 200) :: u_trunc_file !< The complete path to the file in which a column's worth
                              !! of u-accelerations are written if velocity truncations occur.
  character(len = 200) :: v_trunc_file !< The complete path to the file in which a column's worth
                              !! of v-accelerations are written if velocity truncations occur.
  integer :: u_file = -1      !< The unit number for an opened u-truncation file, or -1 if it has
                              !! not been opened.
  integer :: v_file = -1      !< The unit number for an opened v-truncation file, or -1 if it has
                              !! not been opened.
  integer :: cols_written     !< The number of columns whose output has been
                              !! written by this PE during the current run.
  integer :: max_writes       !< The maximum number of times any PE can write out
                              !! a column's worth of accelerations during a run.
  logical :: lemieux_landfast !< If true, use the Lsemieux landfast ice parameterization.
  real :: lemieux_k1          !< 1st free parameter for landfast parameterization [nondim]
  real :: lemieux_k2          !< second free parameter (N/m^3) for landfast parametrization [R L T-2 ~> N m-3]
  real :: lemieux_alphab      !< Cb factor in Lemieux et al 2015 [nondim]
  real :: lemieux_threshold_hw !< max water depth for grounding [Z ~> m]
                              !! see keel data from Amundrud et al. 2004 (JGR)
  real :: lemieux_u0          !< residual velocity for basal stress [L T-1 ~> m s-1]
  logical :: itd_landfast     !< If true, use the probabilistic landfast ice parameterization.
  real :: basal_stress_min_thick !< min ice thickness for grounding [Z ~> m]
  real :: basal_stress_max_depth !< max water depth for grounding [Z ~> m]
  real :: basal_stress_mu_s   !< bottom drag parameter [L Z-1 ~> nondim]
  real :: bathy_roughness_min !< minimum bathymetric roughness [z ~> m]
  real :: bathy_roughness_max !< maximum bathymetric roughness [z ~> m]
  real :: puny                !< small number [nondim]
  real :: onemeter            !< make the units work out (hopefully) [Z ~> m]
  real :: basal_stress_cutoff !< tunable parameter for the bottom drag [nondim]
  integer :: ncat_b           ! number of bathymetry categories
  integer :: ncat_i           ! number of ice thickness categories (log-normal)

  real, pointer, dimension(:,:) :: Tb_u=>NULL() !< Basal stress component at u-points
                                                !! [R Z L T-2 -> kg m-1 s-2]
  real, pointer, dimension(:,:) :: Tb_v=>NULL() !< Basal stress component at v-points
                                                !! [R Z L T-2 -> kg m-1 s-2]
  real, pointer, dimension(:,:) :: sigma_b=>NULL()   !< !< Bottom depth variance [Z ~> m].

  logical :: FirstCall = .true. !< If true, this module has not been called before
  !>@{ Diagnostic IDs
  integer :: id_fix = -1, id_fiy = -1, id_fcx = -1, id_fcy = -1
  integer :: id_fwx = -1, id_fwy = -1, id_sigi = -1, id_sigii = -1
  integer :: id_flfx = -1, id_flfy = -1, id_stren = -1, id_stren0 = -1
  integer :: id_ui = -1, id_vi = -1, id_Coru = -1, id_Corv = -1
  integer :: id_PFu = -1, id_PFv = -1, id_fpx = -1, id_fpy = -1
  integer :: id_fix_d = -1, id_fix_t = -1, id_fix_s = -1
  integer :: id_fiy_d = -1, id_fiy_t = -1, id_fiy_s = -1
  integer :: id_str_d = -1, id_str_t = -1, id_str_s = -1
  integer :: id_sh_d = -1, id_sh_t = -1, id_sh_s = -1
  integer :: id_del_sh = -1, id_del_sh_min = -1
  integer :: id_mis = -1, id_ci = -1, id_ci0 = -1, id_miu = -1, id_miv = -1
  integer :: id_ui_hifreq = -1, id_vi_hifreq = -1
  integer :: id_str_d_hifreq = -1, id_str_t_hifreq = -1, id_str_s_hifreq = -1
  integer :: id_sh_d_hifreq = -1, id_sh_t_hifreq = -1, id_sh_s_hifreq = -1
  integer :: id_sigi_hifreq = -1, id_sigii_hifreq = -1
  integer :: id_stren_hifreq = -1, id_ci_hifreq = -1
  integer :: id_siu = -1, id_siv = -1, id_sispeed = -1 ! SIMIP diagnostics
  integer :: id_ui_east = -1, id_vi_north= -1
  !!@}
end type SIS_C_dyn_CS

end module MOM_SIS_C_dyn_CS_type
