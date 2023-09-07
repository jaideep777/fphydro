module md_phydro_main
  use md_phydro_solver
  use md_phydro_transpiration
  use md_sofunutils

  implicit none

  type phydro_result_type
    real(kind = dbl8) :: a
    real(kind = dbl8) :: e
    real(kind = dbl8) :: gs
    real(kind = dbl8) :: ci
    real(kind = dbl8) :: chi
    real(kind = dbl8) :: vcmax
    real(kind = dbl8) :: jmax
    real(kind = dbl8) :: dpsi
    real(kind = dbl8) :: psi_l
    real(kind = dbl8) :: nfnct
    real(kind = dbl8) :: niter
    real(kind = dbl8) :: mc
    real(kind = dbl8) :: mj
    real(kind = dbl8) :: gammastar
    real(kind = dbl8) :: kmm
    real(kind = dbl8) :: vcmax25
    real(kind = dbl8) :: jmax25
    real(kind = dbl8) :: rd
    real(kind = dbl8) :: isVcmaxLimited
    real(kind = dbl8) :: ac
    real(kind = dbl8) :: aj
    real(kind = dbl8) :: le
    real(kind = dbl8) :: le_s_wet
  end type phydro_result_type

  type par_control_type
    integer(kind = int4) :: gs_method       = GS_IGF
    integer(kind = int4) :: et_method       = ET_DIFFUSION
    integer(kind = int4) :: ftemp_vj_method = FV_kumarathunge19
    integer(kind = int4) :: ftemp_rd_method = FR_heskel16
    integer(kind = int4) :: ftemp_br_method = FB_atkin15
    integer(kind = int4) :: scale_alpha     = 0
  end type par_control_type


  contains

  function phydro_analytical(tc, tg, ppfd, netrad, vpd, co2, elv, fapar, kphio, psi_soil, rdark, vwind, &
                             par_plant, par_cost, par_control) result(res)
    real(kind=dbl8), intent(in) :: tc, tg, ppfd, netrad, vpd, co2, elv, fapar, kphio, psi_soil, rdark, vwind
    type(par_plant_type), intent(in) :: par_plant
    type(par_cost_type), intent(inout) :: par_cost
    type(par_control_type), intent(in) :: par_control
    type(par_env_type) :: par_env
    type(par_photosynth_type) :: par_photosynth
    type(phydro_result_type) :: res

    real(kind=dbl8) :: pa, e, gs, gsprime, x, J, jmax, vcmax, a, dpsi_opt
    type(dpsi_bounds_type) :: bounds
  
    pa = calc_patm(real(elv))
    call create_par_photosynth(par_photosynth, tc, pa, kphio, co2, ppfd, fapar, rdark, tg, par_plant%tchome, &
                               par_control%ftemp_vj_method, par_control%ftemp_rd_method, par_control%ftemp_br_method)
    call create_par_env(par_env, tc, pa, vpd, netrad, vwind)
    
    if (par_control%scale_alpha > 0) par_cost%alpha = par_cost%alpha / par_photosynth%fT_jmax  ! Convert alpha from cost of jmax to cost of jmax25
    par_env%gs_method = par_control%gs_method
    par_env%et_method = par_control%et_method
    
    bounds = calc_dpsi_bound(dble(psi_soil), par_plant, par_env, par_photosynth, par_cost)
    dpsi_opt = zero(bounds%Iabs_bound * 0.001, bounds%Iabs_bound * 0.999, profit_fun, 1.0d-6)
    
    e = calc_sapflux(dpsi_opt, dble(psi_soil), par_plant, par_env)
    gs = calc_gs_from_Q(e, dble(psi_soil), par_plant, par_env)
    gsprime = calc_gsprime(dpsi_opt, gs, dble(psi_soil), par_plant, par_env)
    x = calc_x_from_dpsi(dpsi_opt, gsprime, par_photosynth, par_cost)
    J = calc_J(gs, x, par_photosynth)
    jmax = calc_jmax_from_J(J, par_photosynth)
    vcmax = (J / 4.0d0) * (x * par_photosynth%ca + par_photosynth%kmm) / (x * par_photosynth%ca + 2.0d0 * par_photosynth%gammastar)
    a = gs * (par_photosynth%ca / par_photosynth%patm * 1.0d6) * (1.0d0 - x)
    
    res%a = a
    res%e = e
    res%ci = x * par_photosynth%ca
    res%gs = gs
    res%chi = x
    res%vcmax = vcmax
    res%jmax = jmax
    res%dpsi = dpsi_opt
    res%psi_l = psi_soil - dpsi_opt
    res%nfnct = -999
    res%mc = (x * par_photosynth%ca - par_photosynth%gammastar) / (x * par_photosynth%ca + par_photosynth%kmm)
    res%mj = (x * par_photosynth%ca - par_photosynth%gammastar) / (x * par_photosynth%ca + 2.0d0 * par_photosynth%gammastar)
    res%gammastar = par_photosynth%gammastar
    res%kmm = par_photosynth%kmm
    res%vcmax25 = vcmax / par_photosynth%fT_vcmax
    res%jmax25 = jmax / par_photosynth%fT_jmax
    res%rd = vcmax * par_photosynth%delta
    res%isVcmaxLimited = 0.5d0
    res%ac = a
    res%aj = a
    res%le = e * 0.018015d0 * par_env%lv
    res%le_s_wet = (1.0d0 - fapar) * netrad * (par_env%epsilon / (1.0d0 + par_env%epsilon))
  
    contains

    function profit_fun(dpsi)
      real(kind = dbl8), intent(in) :: dpsi
      real(kind = dbl8) :: profit_fun
      type(dfdx_type) :: dfdx_res
      dfdx_res = dFdx(dpsi, dble(psi_soil), par_plant, par_env, par_photosynth, par_cost)
      profit_fun = dfdx_res%dPdx
    end  

  end function phydro_analytical
  
  function phydro_instantaneous_analytical(vcmax25, jmax25, tc, tg, ppfd, netrad, vpd, co2, elv, &
                                           fapar, kphio, psi_soil, rdark, vwind, par_plant, par_cost, par_control) result(res)
    real(kind=dbl8), intent(in) :: vcmax25, jmax25, tc, tg, ppfd, netrad, vpd, co2, elv, fapar, kphio, psi_soil, rdark, vwind
    type(par_plant_type), intent(in) :: par_plant
    type(par_cost_type), intent(inout) :: par_cost
    type(par_control_type), intent(in) :: par_control
    type(par_env_type) :: par_env
    type(par_photosynth_type) :: par_photosynth
    type(phydro_result_type) :: res
    real(kind=dbl8) :: pa, e, gs
    real(kind=dbl8) :: bound, jmax, vcmax, dpsi_opt
    type(ACi_type) :: Aa, Ac, Aj
    
    pa = calc_patm(real(elv))
    call create_par_photosynth(par_photosynth, tc, pa, kphio, co2, ppfd, fapar, rdark, tg, par_plant%tchome, &
                               par_control%ftemp_vj_method, par_control%ftemp_rd_method, par_control%ftemp_br_method)
    call create_par_env(par_env, tc, pa, vpd, netrad, vwind)
    
    ! call print_par_photosynth(par_photosynth)

    ! optionally convert alpha from cost of jmax to cost of jmax25
    if (par_control%scale_alpha > 0) par_cost%alpha = par_cost%alpha / par_photosynth%fT_jmax  !
    par_env%gs_method = par_control%gs_method
    par_env%et_method = par_control%et_method
    
    vcmax = vcmax25 * par_photosynth%fT_vcmax
    jmax = jmax25 * par_photosynth%fT_jmax
    
    bound = calc_dpsi_bound_inst(psi_soil, par_plant, par_env, par_photosynth, par_cost)
    dpsi_opt = zero(0.0d0, 0.99d0 * bound, profit_fun_inst, 1.0d-6)
    
    e = calc_sapflux(dpsi_opt, psi_soil, par_plant, par_env)
    gs = calc_gs_from_Q(e, psi_soil, par_plant, par_env)
    Aa = calc_assimilation_limiting(vcmax, jmax, gs, par_photosynth)
    Ac = calc_assim_rubisco_limited(gs, vcmax, par_photosynth)
    Aj = calc_assim_light_limited(gs, jmax, par_photosynth)
    
    res%a = Aa%a
    res%e = e
    res%ci = Aa%ci
    res%gs = gs
    res%chi = Aa%ci / par_photosynth%ca
    res%vcmax = vcmax
    res%jmax = jmax
    res%dpsi = dpsi_opt
    res%psi_l = psi_soil - dpsi_opt
    res%mc = (Aa%ci - par_photosynth%gammastar) / (Aa%ci + par_photosynth%kmm)
    res%mj = (Aa%ci - par_photosynth%gammastar) / (Aa%ci + 2.0d0 * par_photosynth%gammastar)
    res%gammastar = par_photosynth%gammastar
    res%kmm = par_photosynth%kmm
    res%vcmax25 = vcmax25
    res%jmax25 = jmax25
    res%rd = vcmax * par_photosynth%delta
    res%isVcmaxLimited = merge(1.d0, 0.d0, Aa%isVcmaxLimited)
    res%ac = Ac%a
    res%aj = Aj%a
    res%le = e * 0.018015d0 * par_env%lv
    res%le_s_wet = (1.0d0 - fapar) * netrad * (par_env%epsilon / (1.0d0 + par_env%epsilon))
  
    contains

    function profit_fun_inst(dpsi)
      real(kind = dbl8), intent(in) :: dpsi
      real(kind = dbl8) :: profit_fun_inst
      profit_fun_inst = calc_dP_ddpsi(dpsi, vcmax, jmax, psi_soil, par_plant, par_env, par_photosynth, par_cost)
    end  


  end function phydro_instantaneous_analytical

end module md_phydro_main
