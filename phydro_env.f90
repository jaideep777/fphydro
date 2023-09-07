module md_phydro_env
  use md_precision
  use md_phydro_physical
  use md_photosynth, only: calc_viscosity_h2o, calc_density_h2o
  implicit none

  ! list of methods to calculate gs
  integer (kind = int4), parameter :: GS_IGF = 0, GS_QNG = 1, GS_APX = 2, GS_APX2 = 3
  
  integer (kind = int4), parameter :: ET_DIFFUSION = 0, ET_PM = 1

  ! Define the data type for ParEnv
  type par_env_type
    real(kind = dbl8) :: tc          ! Temperature [degC]
    real(kind = dbl8) :: patm        ! Atmospheric pressure [Pa]
    real(kind = dbl8) :: vpd         ! VPD [Pa]
    real(kind = dbl8) :: Rn          ! Net radiation [W m-2]
    real(kind = dbl8) :: v_wind      ! Wind speed [m s-1]
    real(kind = dbl8) :: viscosity_water  ! [Pa s]
    real(kind = dbl8) :: density_water    ! [kg m-3]
    real(kind = dbl8) :: rho         ! Density of air [kg m-3]
    real(kind = dbl8) :: cp          ! Specific heat capacity of moist air [J kg-1 K-1]
    real(kind = dbl8) :: gamma       ! Psychrometric constant [Pa K-1]
    real(kind = dbl8) :: epsilon     ! Slope of saturation-pressure - temp curve [Pa K-1]
    real(kind = dbl8) :: lv          ! Latent heat of vaporization of water [J kg-1]
    integer(kind = int4) :: gs_method = GS_IGF ! GsMethod
    integer(kind = int4) :: et_method = ET_DIFFUSION  ! ETMethod
  end type par_env_type

  ! Interface for member subroutines
  interface par_env_type_interface
    module procedure :: create_par_env
    ! module procedure :: calc_temp_dependencies
    ! module procedure :: print_par_env
  end interface

  contains

  ! Constructor for ParEnv
  subroutine create_par_env(this, tc, patm, vpd, Rn, v_wind) 
    type(par_env_type), intent(inout) :: this
    real(kind = dbl8), intent(in) :: tc, patm, vpd, Rn, v_wind
    this%tc = tc
    this%vpd = vpd
    this%patm = patm
    this%Rn = Rn
    this%v_wind = v_wind
    this%gs_method = GS_IGF 
    this%et_method = ET_DIFFUSION
    call calc_temp_dependencies(this)
  end subroutine create_par_env

  ! Separate constructor without v_wind as a parameter
  subroutine create_par_env_no_wind(this, tc, patm, vpd, Rn)
    type(par_env_type), intent(out) :: this
    real(kind = dbl8), intent(in) :: tc, patm, vpd, Rn
    call create_par_env(this, tc, patm, vpd, Rn, 3.0d0) ! Default v_wind
  end subroutine create_par_env_no_wind

  ! Calculate temperature dependencies
  subroutine calc_temp_dependencies(this)
    type(par_env_type), intent(inout) :: this
    this%viscosity_water = calc_viscosity_h2o(real(this%tc), real(this%patm))
    this%density_water = calc_density_h2o(real(this%tc), real(this%patm))
    this%rho = calc_density_air(this%tc, this%patm, this%vpd, .true.)
    this%cp = calc_cp_moist_air(this%tc)
    this%gamma = calc_psychro(this%tc, this%patm)
    this%epsilon = calc_sat_slope(this%tc) / this%gamma
    this%lv = calc_enthalpy_vap(this%tc)
  end subroutine calc_temp_dependencies

  ! Print ParEnv information
  subroutine print_par_env(this)
    type(par_env_type), intent(in) :: this
    write(*, *) "Env:"
    write(*, *) "   tc = ", this%tc, " [degC]"
    write(*, *) "   patm = ", this%patm, " [Pa]"
    write(*, *) "   vpd = ", this%vpd, " [Pa]"
    write(*, *) "   Rn = ", this%Rn, " [W m-2]"
    write(*, *) "   v_wind = ", this%v_wind, " [m s-1]"
    write(*, *) "   viscosity_water = ", this%viscosity_water, " [Pa s]"
    write(*, *) "   density_water = ", this%density_water, " [kg m-3]"
    write(*, *) "   rho = ", this%rho, " [kg m-3]"
    write(*, *) "   cp = ", this%cp, " [J kg-1 K-1]"
    write(*, *) "   gamma = ", this%gamma, " [Pa K-1]"
    write(*, *) "   epsilon = ", this%epsilon, " [Pa K-1]"
    write(*, *) "   lv = ", this%lv, " [J kg-1]"
  end subroutine print_par_env

end module md_phydro_env
  