module md_phydro_solver
  use md_precision
  use md_phydro_photosynthesis
  use md_phydro_transpiration
  
  implicit none

  type par_cost_type
    real (kind = dbl8) :: alpha
    real (kind = dbl8) :: gamma
  end type par_cost_type
 	
  type dpsi_bounds_type
    real (kind = dbl8) ::  exact
    real (kind = dbl8) ::  approx_O2
    real (kind = dbl8) ::  Iabs_bound
  end type dpsi_bounds_type

  type dfdx_type
    real (kind = dbl8) ::  dPdx
    real (kind = dbl8) ::  J
    real (kind = dbl8) ::  djmax_dJ
    real (kind = dbl8) ::  dJ_dchi
  end type dfdx_type


  contains


  !!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !!! Phydro analytical solver (acclimating)
  !!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function calc_J(gs, x, par_photosynth) result(J)
    real (kind = dbl8), intent(in) :: gs, x
    type(par_photosynth_type), intent(in) :: par_photosynth
    real (kind = dbl8) :: g, K, ca, d, J
    g = par_photosynth%gammastar / par_photosynth%ca
    k = par_photosynth%kmm / par_photosynth%ca
    ca = par_photosynth%ca / par_photosynth%patm*1e6
    d = par_photosynth%delta
    J = 4*gs*ca*(1-x)*(x+ 2*g)/(x*(1-d)-(g+d*k))
  end


  function calc_jmax_from_J(J, par_photosynth) result(jmax)
    real (kind = dbl8), intent(in) :: J
    type(par_photosynth_type), intent(in) :: par_photosynth
    real (kind = dbl8) :: pp, pj, jmax
    pp = 4*par_photosynth%phi0 * par_photosynth%Iabs;
    pj = pp/J;
    jmax = pp/sqrt(pj*pj-1);
  end


  function calc_djmax_dJ(J, par_photosynth) result(djdj)
    real (kind = dbl8), intent(in) :: J
    type(par_photosynth_type), intent(in) :: par_photosynth
    real (kind = dbl8) :: pp, sq, psq, djdj
    pp = 4*par_photosynth%phi0 * par_photosynth%Iabs
    sq = sqrt(pp*pp-J*J)
    psq = pp/sq
    djdj = psq*psq*psq
  end


  function calc_dJ_dchi(gs, x, par_photosynth) result(djdx)
    real (kind = dbl8), intent(in) :: gs, x
    type(par_photosynth_type), intent(in) :: par_photosynth
    real (kind = dbl8) :: g, K, ca, d, djdx, d1
    g = par_photosynth%gammastar / par_photosynth%ca
    k = par_photosynth%kmm / par_photosynth%ca
    ca = par_photosynth%ca / par_photosynth%patm*1e6
    d = par_photosynth%delta
    ! gs*ca * ((d*(2*g*(k + 1) + k*(2*x - 1) + x^2) + 2*g^2 + g*(2*x - 3) - x^2)/(d*(k + x) + g - x)^2)
    d1 = d*(k + x) + g - x;
    djdx = 4*gs*ca * ((d*(2*g*(k + 1) + k*(2*x - 1) + x*x) - ((x-g)*(x-g)+3*g*(1-g)))/(d1*d1));
    ! gs*ca*(3*(g-1)*g/(g-x)^2 - 1)
  end


  function calc_dJ_ddpsi(gsprime, x, par_photosynth) result(djdp)
    real (kind = dbl8), intent(in) :: gsprime, x
    type(par_photosynth_type), intent(in) :: par_photosynth
    real (kind = dbl8) :: g, K, ca, d, djdp
    g = par_photosynth%gammastar / par_photosynth%ca
    k = par_photosynth%kmm / par_photosynth%ca
    ca = par_photosynth%ca / par_photosynth%patm*1e6
    d = par_photosynth%delta
    djdp = 4*gsprime*ca*(1-x)*(x+2*g)/(x*(1-d)-(g+d*k))
  end


  function calc_x_from_dpsi(dpsi, gsprime, par_photosynth, par_cost) result(x)
    real (kind = dbl8), intent(in) :: dpsi, gsprime
    type(par_photosynth_type), intent(in) :: par_photosynth
    type(par_cost_type), intent(in) :: par_cost

    real (kind = dbl8) gstar, Km, ca, br, y, ca2, x

    gstar = par_photosynth%gammastar/par_photosynth%patm*1e6
    Km = par_photosynth%kmm/par_photosynth%patm*1e6
    ca = par_photosynth%ca/par_photosynth%patm*1e6
    br = par_photosynth%delta
    y = par_cost%gamma
    
    ca2 = ca*ca;
    x = (-2*ca*dpsi*(gstar + br*Km)*y + &
      ca2*((3 - 2*br)*gstar + br*Km)*gsprime + &
      -sqrt(2.0D+00)*sqrt( &
        ca2*dpsi*((-3 + 2*br)*gstar - br*Km)*((-1 + br)*ca + gstar + &
                                                  br*Km)*y* &
          (-2*dpsi*y + (ca + 2*gstar)* &
              gsprime)))/ &
      (ca2*(2*(-1 + br)*dpsi*y + ((3 - 2*br)*gstar + br*Km)* &
              gsprime))
    
    if (x < (gstar + br*Km)/(ca - br*ca)) x = (gstar + br*Km)/(ca - br*ca)+1e-12
  end
  

  function dFdx(dpsi, psi_soil, par_plant, par_env, par_photosynth, par_cost) result(res)
    real (kind = dbl8), intent(in) :: dpsi, psi_soil
    type(par_plant_type), intent(in) :: par_plant
    type(par_env_type), intent(in) :: par_env
    type(par_photosynth_type), intent(in) :: par_photosynth
    type(par_cost_type), intent(in) :: par_cost

    real (kind = dbl8) :: Q, gs, gsprime, X, J, ca, g, djmax_dJ, dJ_dchi, dP_dx
    type(dfdx_type) :: res

    Q = calc_sapflux(dpsi, psi_soil, par_plant, par_env)
    gs = calc_gs_from_Q(Q, psi_soil, par_plant, par_env)
    gsprime = calc_gsprime(dpsi, gs, psi_soil, par_plant, par_env)
  
    X =  calc_x_from_dpsi(dpsi, gsprime, par_photosynth, par_cost)

    J = calc_J(gs, X, par_photosynth)

    ca = par_photosynth%ca / par_photosynth%patm*1e6
    g = par_photosynth%gammastar / par_photosynth%ca

    djmax_dJ = calc_djmax_dJ(J, par_photosynth)
    dJ_dchi  = calc_dJ_dchi(gs, X, par_photosynth)

    dP_dx = -gs*ca - par_cost%alpha * djmax_dJ * dJ_dchi

    res = dfdx_type(dP_dx, J, djmax_dJ, dJ_dchi)
  end


  function calc_dpsi_bound(psi_soil, par_plant, par_env, par_photosynth, par_cost) result(bounds)
    real (kind = dbl8), intent(in) :: psi_soil
    type(par_plant_type), intent(in) :: par_plant
    type(par_env_type), intent(in) :: par_env
    type(par_photosynth_type), intent(in) :: par_photosynth
    type(par_cost_type), intent(in) :: par_cost

    type(dpsi_bounds_type) :: bounds 

    real (kind = dbl8) :: gstar, ca, y, K, Pox, Ppox, Pppox
    real (kind = dbl8) :: a,b,c,del
    real (kind = dbl8) :: ex, appo2, iabsb, use_bound
    real (kind = dbl8) :: ga, Qmax, max_dpsi

    gstar = par_photosynth%gammastar/par_photosynth%patm*1e6
    ca = par_photosynth%ca/par_photosynth%patm*1e6
    y = par_cost%gamma

    K = scale_conductivity(par_plant%conductivity, par_env)/(1.6*par_env%vpd/par_env%patm);
 
    Pox = P(psi_soil, par_plant%psi50, par_plant%b);
    Ppox = Pprime(psi_soil, par_plant%psi50, par_plant%b);
    Pppox = Pprimeprime(psi_soil, par_plant%psi50, par_plant%b);
    
    a = (ca + 2*gstar)*K*Pppox*4.0d0/8.0d0;
    b = -(2*y + (ca + 2*gstar)*K*Ppox);
    c = (ca + 2*gstar)*K*Pox;
    del = b*b-4*a*c;

    appo2 = (-b-sqrt(del))/(2*a)
    ex = zero(0.0d0, 10.0d0, f2, 1d-6)

    use_bound = ex
    
    iabsb = zero(use_bound*0.001, use_bound*0.99, f1, 1D-6);
       
    ! If using PM, find max dpsi from max possible transpiration 
    if (par_env%et_method == ET_PM) then
      ga = calc_g_aero(par_plant%h_canopy, dble(par_env%v_wind), par_plant%h_wind_measurement);
      Qmax = calc_max_transpiration_pm(ga, par_env);
      max_dpsi = calc_dpsi_from_sapflux(Qmax, psi_soil, par_plant, par_env);
      iabsb = min(max_dpsi, iabsb);
    endif
  

    bounds = dpsi_bounds_type(ex, appo2, iabsb)

    contains

    function f2(dpsi) result(gg)
      real(kind = dbl8), intent(in) :: dpsi
      real(kind = dbl8) :: gg, gsprime
      gsprime = calc_gsprime_from_dpsi(dpsi, psi_soil, par_plant, par_env)
      gg = (-2*dpsi*y + (ca + 2*gstar)*gsprime)
    end

    function f1(dpsi) result(J)
      real(kind = dbl8), intent(in) :: dpsi
      real(kind = dbl8) :: J, gs, x, Q, gsprime
      Q = calc_sapflux(dpsi, psi_soil, par_plant, par_env);
      gs = calc_gs_from_Q(Q, psi_soil, par_plant, par_env);
      gsprime = calc_gsprime(dpsi, gs, psi_soil, par_plant, par_env);
      x = calc_x_from_dpsi(dpsi,gsprime, par_photosynth, par_cost);
      J = calc_J(gs, x, par_photosynth)-4.0d0*par_photosynth%phi0*par_photosynth%Iabs;
    end

  end



  !!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !!! Phydro analytical solver (instantaneous)
  !!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function calc_dP_ddpsi(dpsi, vcmax, jmax, psi_soil, par_plant, par_env, par_photosynth, par_cost) result(dP_ddpsi)
    real(kind=dbl8), intent(in) :: dpsi, vcmax, jmax, psi_soil
    type(par_plant_type), intent(in) :: par_plant
    type(par_env_type), intent(in) :: par_env
    type(par_photosynth_type), intent(in) :: par_photosynth
    type(par_cost_type), intent(in) :: par_cost
    real(kind=dbl8) :: gstar, Km, ca, br, y, Q, gs, P, dpsi1, Q1, gs1, P1
    type(ACi_type) :: Assim, Assim1
    real(kind=dbl8) :: dP_ddpsi
  
    gstar = par_photosynth%gammastar / par_photosynth%patm * 1.0d6
    Km = par_photosynth%kmm / par_photosynth%patm * 1.0d6
    ca = par_photosynth%ca / par_photosynth%patm * 1.0d6
    br = par_photosynth%delta
    y = par_cost%gamma
  
    Q = calc_sapflux(dpsi, psi_soil, par_plant, par_env)
    gs = calc_gs_from_Q(Q, psi_soil, par_plant, par_env)
    Assim = calc_assimilation_limiting(vcmax, jmax, gs, par_photosynth)
    P = Assim%a - y * dpsi * dpsi
  
    dpsi1 = dpsi + 1.0d-6
    Q1 = calc_sapflux(dpsi1, psi_soil, par_plant, par_env)
    gs1 = calc_gs_from_Q(Q1, psi_soil, par_plant, par_env)
    Assim1 = calc_assimilation_limiting(vcmax, jmax, gs1, par_photosynth)
    P1 = Assim1%a - y * (dpsi1) * (dpsi1)
  
    dP_ddpsi = (P1 - P) / 1.0d-6
  
  end function calc_dP_ddpsi
  

  function calc_dpsi_bound_inst(psi_soil, par_plant, par_env, par_photosynth, par_cost) result(bound)
    real(kind=dbl8), intent(in) :: psi_soil
    type(par_plant_type), intent(in) :: par_plant
    type(par_env_type), intent(in) :: par_env
    type(par_photosynth_type), intent(in) :: par_photosynth
    type(par_cost_type), intent(in) :: par_cost
    real(kind=dbl8) :: bound, ga, Qmax, max_dpsi
  
    bound = 100.0d0
  
    ! If using PM, find max dpsi from max possible transpiration 
    if (par_env%et_method == ET_PM) then
      ga = calc_g_aero(par_plant%h_canopy, dble(par_env%v_wind), par_plant%h_wind_measurement)
      Qmax = calc_max_transpiration_pm(ga, par_env)
      max_dpsi = calc_dpsi_from_sapflux(Qmax, psi_soil, par_plant, par_env)
      bound = min(max_dpsi, bound)
    end if
  
  end function calc_dpsi_bound_inst


end module md_phydro_solver

