module md_phydro_pm
  use md_phydro_env
  implicit none

  contains

  function calc_g_aero(h_canopy, v_wind, z_measurement) result(g_aero)
    ! Aerodynamic conductance [m s-1]
    ! To convert to mol m-2 s-1, see this: https://rdrr.io/cran/bigleaf/man/ms.to.mol.html (but not convincing)
    ! Refs: 
    !    Eq 13 in Leuning et al (2008). https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2007WR006562
    !    Eq 7 in Zhang et al (2008): https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2017JD027025
    !    Box 4 in https://www.fao.org/3/x0490e/x0490e06.htm 
    real(kind = dbl8), intent(in) :: h_canopy, v_wind, z_measurement
    real(kind = dbl8) :: g_aero, k_karman, d, z_om, z_ov
    
    k_karman = 0.41        ! von Karman's constant [-]
    d = h_canopy * 2.0 / 3.0   ! zero-plane displacement height [m]
    z_om = 0.123 * h_canopy    ! roughness lengths governing transfer of water and momentum [m]
    z_ov = 0.1 * z_om
    
    g_aero = (k_karman * k_karman * v_wind) / (log((z_measurement - d) / z_om) * log((z_measurement - d) / z_ov))
  end function calc_g_aero
    

  function gs_conv(tc, patm) result(gs_conv_value)
    ! multiplier to convert:
    !   stomatal conductance to CO2 [mol m-2 s-1] ----> stomatal conductance to water [m s-1]
    real(kind = dbl8), intent(in) :: tc, patm
    real(kind = dbl8) :: gs_conv_value, R
    
    R = 8.31446261815324 ! Universal gas constant [J mol-1 K-1]
    
    gs_conv_value = 1.6 * R * (tc + 273.16) / patm
  end function gs_conv
    

  function calc_transpiration_pm(gs, ga, par_env) result(trans)
    ! Calculate PML transpiration [mol m-2 s-1]
    ! gs   Stomatal conductance to CO2 [mol m-2 s-1]
    ! ga   Aerodynamic conductance [m s-1]
    ! Rn   Absorbed net radiation [W m-2]
    real(kind = dbl8), intent(in) :: gs, ga
    type(par_env_type), intent(in) :: par_env
    real(kind = dbl8) :: trans, gw, latent_energy
    
    gw = gs * gs_conv(par_env%tc, par_env%patm)  ! gw in [m s-1]
    
    latent_energy = (par_env%epsilon * par_env%Rn + (par_env%rho * par_env%cp / par_env%gamma) &
                    * ga * par_env%vpd) / (par_env%epsilon + 1 + ga / gw) ! latent energy W m-2 
    trans = latent_energy * (55.5 / par_env%lv) ! W m-2 ---> mol m-2 s-1
  end function calc_transpiration_pm


  function calc_max_transpiration_pm(ga, par_env) result(trans_max)
    ! Calculate maximum possible PML transpiration for a given ga, calculated by setting gs = inf, [mol m-2 s-1]
    ! ga   Aerodynamic conductance [m s-1]
    ! Rn   Absorbed net radiation [W m-2]
    real(kind = dbl8), intent(in) :: ga
    type(par_env_type), intent(in) :: par_env
    real(kind = dbl8) :: trans_max, latent_energy

    latent_energy = (par_env%epsilon * par_env%Rn + (par_env%rho * par_env%cp / par_env%gamma) &
                    * ga * par_env%vpd) / (par_env%epsilon + 1) ! latent energy W m-2 
    trans_max = latent_energy * (55.5 / par_env%lv) ! W m-2 ---> mol m-2 s-1
  end function calc_max_transpiration_pm


  function calc_gs_pm(Q, ga, par_env) result(gs)
    ! Calculate PML stomatal conductance to CO2 [mol m-2 s-1]
    ! Q    Sap flux [mol m-2 s-1]
    ! ga   Aerodynamic conductance [m s-1]
    ! Rn   Absorbed net radiation [W m-2]
    real(kind = dbl8), intent(in) :: Q, ga
    type(par_env_type), intent(in) :: par_env
    real(kind = dbl8) :: gs, Q_energy, den, gw

    Q_energy = Q * (par_env%lv / 55.5)

    den = par_env%epsilon * par_env%Rn + (par_env%rho * par_env%cp / par_env%gamma) &
          * ga * par_env%vpd - (1 + par_env%epsilon) * Q_energy
    !den = fmax(den, 0)

    gw = ga * Q_energy / den ! stomatal conductance to water [m s-1]

    gs = gw / gs_conv(par_env%tc, par_env%patm) ! stomatal conductance to CO2 [mol m-2 s-1]
  end function calc_gs_pm


  function calc_dE_dgs_pm(gs, ga, par_env) result(dE_dgs)
    ! Calculate derivative of transpiration wrt stomatal conductance to CO2 [unitless] - analytical version
    real(kind = dbl8), intent(in) :: gs, ga
    type(par_env_type), intent(in) :: par_env
    real(kind = dbl8) :: dE_dgs, gw, num, den, d_le_dgw

    gw = gs * gs_conv(par_env%tc, par_env%patm)  ! [m s-1]

    num = ga * (par_env%epsilon * par_env%Rn + (par_env%rho * par_env%cp / par_env%gamma) * ga * par_env%vpd)
    den = par_env%epsilon * gw + gw + ga

    d_le_dgw = (num / den / den) ! derivative of latent energy wrt stomatal conductance for water in m s-1

    dE_dgs = d_le_dgw * (55.5 / par_env%lv) * gs_conv(par_env%tc, par_env%patm)
  end function calc_dE_dgs_pm


  function calc_dE_dgs_pm_num(gs, ga, par_env) result(dE_dgs)
    ! Calculate derivative of transpiration wrt stomatal conductance to CO2 [unitless] - numerical version
    real(kind = dbl8), intent(in) :: gs, ga
    type(par_env_type), intent(in) :: par_env
    real(kind = dbl8) :: dE_dgs, E, E_plus

    E = calc_transpiration_pm(gs, ga, par_env)
    E_plus = calc_transpiration_pm(gs + 1.0e-6, ga, par_env)

    dE_dgs = (E_plus - E) / 1.0e-6
  end function calc_dE_dgs_pm_num


end module md_phydro_pm
  