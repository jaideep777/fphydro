module md_phydro_transpiration
  use md_precision
  use md_phydro_pm
  use md_phydro_env
  use md_sofunutils_plus, only: zero, gammad

  implicit none

  type par_plant_type
    real (kind = dbl8) :: conductivity        ! = ci/ca, leaf-internal to ambient CO2 partial pressure, ci/ca (unitless)
    real (kind = dbl8) :: psi50               ! leaf-internal CO2 partial pressure (Pa)
    real (kind = dbl8) :: b                   ! ci-limitation factor of light-limited assimilation (unitless)

    real (kind = dbl8) :: h_canopy = 20
    real (kind = dbl8) :: h_wind_measurement = 22
    real (kind = dbl8) :: tchome = 25

    integer (kind = int4) :: gs_method = GS_IGF
  end type par_plant_type


  contains

  ! Constructors for par plant
  subroutine init_par_plant(this, cond, psi, b)
    class(par_plant_type), intent(out) :: this
    real(kind=dbl8), intent(in) :: cond, psi, b
    this%conductivity = cond
    this%psi50 = psi
    this%b = b
  end subroutine init_par_plant

  subroutine init_par_plant_6args(this, cond, psi, b, hcanopy, hwind, tchome)
    class(par_plant_type), intent(out) :: this
    real(kind=dbl8), intent(in) :: cond, psi, b, hcanopy, hwind, tchome
    this%conductivity = cond
    this%psi50 = psi
    this%b = b
    this%h_canopy = hcanopy
    this%h_wind_measurement = hwind
    this%tchome = tchome
  end subroutine init_par_plant_6args   


  !!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !!! Vulnerability curve
  !!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function P(psi, psi50, b)
    real (kind = dbl8), intent(in) :: psi
    real (kind = dbl8), intent(in) :: psi50
    real (kind = dbl8), intent(in) :: b
    real (kind = dbl8) :: P
    P = 0.5 ** ((psi/psi50) ** b)
  end

  function Pprime(psi, psi50, b)
    real (kind = dbl8), intent(in) :: psi
    real (kind = dbl8), intent(in) :: psi50
    real (kind = dbl8), intent(in) :: b
    real (kind = dbl8) :: Pprime
    Pprime = log(0.5) * P(psi,psi50,b) * b * ((psi/psi50)**(b-1)) / psi50
  end

  function Pprimeprime(psi, psi50, b)
    real (kind = dbl8), intent(in) :: psi
    real (kind = dbl8), intent(in) :: psi50
    real (kind = dbl8), intent(in) :: b
    real (kind = dbl8) :: Pprimeprime
    Pprimeprime = log(0.5)*b*((psi/psi50)**(b-1))/psi50 * Pprime(psi, psi50, b) &
                + log(0.5)*P(psi, psi50, b)/(psi50*psi50)*b*(b-1)* ((psi/psi50)**(b-2))
  end

  !!! Convert conductivity from m (m3/m2) to mol/m2/s/Mpa
  function scale_conductivity(K, par_env) result(K4)
    real (kind = dbl8), intent(in) :: K
    type(par_env_type), intent(in) :: par_env
    real (kind = dbl8) :: K2, K3, K4
    real (kind = dbl8) :: mol_h20_per_kg_h20 = 55.5

    ! Flow rate in m3/m2/s/Pa
    K2 = K/par_env%viscosity_water
  
    ! Flow rate in mol/m2/s/Pa
    K3 = K2 * par_env%density_water * mol_h20_per_kg_h20;
    
    ! Flow rate in mol/m2/s/Mpa
    K4 = K3 * 1e6;
  end function scale_conductivity
  

  !!! integrate vulnerability curve
  function integral_P_analytical(dpsi, psi_soil, psi50, b) result(I)
  ! int P(p, p50, b) = -(p/b) * (log2)^(-1/b) * G(1/b, (x/p)^b*log2)  <--- G is unnormalized upper incomplete gamma function
  !                  = -(p/b) * (log2)^(-1/b) * G(1/b) * (1 - I((x/p)^b*log2) <--- I is lower incomplete gamma integral
  !                  = -(p/b) * (log2)^(-1/b) * G(1/b) * (- I((pl/p)^b*log2 + I((ps/p)^b*log2) <--- I is lower incomplete gamma integral
  !                  = +(p/b) * (log2)^(-1/b) * G(1/b) * (  I((pl/p)^b*log2 - I((ps/p)^b*log2) <--- I is lower incomplete gamma integral
    real (kind = dbl8), intent(in) :: dpsi, psi_soil, psi50, b
    real (kind = dbl8) :: I, ps, pl, l2
    integer (kind = int4) :: ifault
    ps = psi_soil/psi50;
    pl = (psi_soil-dpsi)/psi50;
    l2 = log(2.0);
    I = (psi50/b) * (l2**(-1/b)) * gamma(1/b) * (gammad(l2*(pl**b), 1/b, ifault) - gammad(l2*(ps**b), 1/b, ifault))
  end


  function integral_P_approx(dpsi, psi_soil, psi50, b) result(I)
    real (kind = dbl8), intent(in) :: dpsi, psi_soil, psi50, b
    real (kind = dbl8) :: I
    I = -P(psi_soil-dpsi/2.0, psi50, b)*dpsi
  end


  function integral_P_approx2(dpsi, psi_soil, psi50, b) result(I)
    real (kind = dbl8), intent(in) :: dpsi, psi_soil, psi50, b
    real (kind = dbl8) :: I
    I = -(P(psi_soil, psi50, b)+P(psi_soil-dpsi, psi50, b))/2 * dpsi
  end


  function integral_P(dpsi, psi_soil, par_plant) result(I)
    real (kind = dbl8), intent(in) :: dpsi, psi_soil
    type(par_plant_type), intent(in) :: par_plant
    real (kind = dbl8) :: I

    ! if      (par_plant%gs_method == GS_QNG)  then; I = integral_P_numerical( dpsi, psi_soil, par_plant%psi50, par_plant%b);
    if      (par_plant%gs_method == GS_IGF)  then; I = integral_P_analytical(dpsi, psi_soil, par_plant%psi50, par_plant%b);
    else if (par_plant%gs_method == GS_APX)  then; I = integral_P_approx(    dpsi, psi_soil, par_plant%psi50, par_plant%b);
    else if (par_plant%gs_method == GS_APX2) then; I = integral_P_approx2(   dpsi, psi_soil, par_plant%psi50, par_plant%b);
    else; error stop "Unsupported gs_method specified"
    end if
  end

  !!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !!! Transpiration and stomatal conductance
  !!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function calc_sapflux(dpsi, psi_soil, par_plant, par_env) result(E)
    real (kind = dbl8), intent(in) :: dpsi, psi_soil
    type(par_plant_type), intent(in) :: par_plant
    type(par_env_type), intent(in) :: par_env
    real (kind = dbl8) :: E, K

    K = scale_conductivity(par_plant%conductivity, par_env)
    E = K * (-integral_P(dpsi, psi_soil, par_plant))
  end

  function calc_max_sapflux(psi_soil, par_plant, par_env) result(E)
    real (kind = dbl8), intent(in) :: psi_soil
    type(par_plant_type), intent(in) :: par_plant
    type(par_env_type), intent(in) :: par_env
    real (kind = dbl8) :: E, K
    
    K = scale_conductivity(par_plant%conductivity, par_env)
    E = K * (-integral_P(1e20_dbl8, psi_soil, par_plant))
  end


  !                                 _ps-dpsi 
  ! Calculate dpsi that solves    _/   K(psi') dpsi' = Q
  !                             ps
  function calc_dpsi_from_sapflux(Q, psi_soil, par_plant, par_env) result(dpsi)
    type(par_plant_type) :: par_plant
    type(par_env_type) :: par_env
    real(kind=dbl8) :: Q, psi_soil, dpsi, Qmax
    
    Qmax = calc_max_sapflux(psi_soil, par_plant, par_env);
    if (Q > Qmax) then
      dpsi = 999999999.0_dbl8
    else
      dpsi = zero(0.0_dbl8, 100.0_dbl8, f, 1e-6_dbl8)
    endif
  
    contains
    
      function f(dpsi) 
        real(kind=dbl8), intent(in) :: dpsi
        real(kind=dbl8) :: f      
        f = calc_sapflux(dpsi, psi_soil, par_plant, par_env) - Q;
      end function f

  end function calc_dpsi_from_sapflux
  

  ! Calculates regulated stomatal conductance given transpiration/sapflux
  ! water balance is assumed
  ! plant hydraulic traits, and the environment.
  function calc_gs_from_Q(Q, psi_soil, par_plant, par_env) result(gs)
    real(dbl8), intent(in) :: Q, psi_soil
    type(par_plant_type), intent(in) :: par_plant
    type(par_env_type), intent(in) :: par_env
    real(dbl8) :: D, gs, ga

    D = (par_env%vpd / par_env%patm)

    if (par_env%et_method == ET_DIFFUSION) then
        gs = Q / (1.6d0 * D)
    else if (par_env%et_method == ET_PM) then
        ga = calc_g_aero(par_plant%h_canopy, dble(par_env%v_wind), par_plant%h_wind_measurement)
        gs = calc_gs_pm(Q, ga, par_env)
    else
        write(*,*) 'Unknown et_method:', par_env%et_method
        stop
    end if
  end function calc_gs_from_Q

  ! Derivative of sapflux wrt dpsi, dQ/ddpsi
  function calc_Qprime_analytical(dpsi, psi_soil, par_plant, par_env) result(Qprime)
    real(dbl8), intent(in) :: dpsi, psi_soil
    type(par_plant_type), intent(in) :: par_plant
    type(par_env_type), intent(in) :: par_env
    real(dbl8) :: K
    real(dbl8) :: Qprime

    K = scale_conductivity(par_plant%conductivity, par_env)
    Qprime = K * P(psi_soil - dpsi, par_plant%psi50, par_plant%b)
  end function calc_Qprime_analytical


  function calc_Qprime_approx(dpsi, psi_soil, par_plant, par_env) result(Qprime)
    type(par_plant_type) :: par_plant
    type(par_env_type) :: par_env
    real(kind=dbl8) :: dpsi, psi_soil, Qprime, K

    K = scale_conductivity(par_plant%conductivity, par_env)
    Qprime = K * (P(psi_soil - dpsi / 2, par_plant%psi50, par_plant%b) - &
                  Pprime(psi_soil - dpsi / 2, par_plant%psi50, par_plant%b) * dpsi / 2)
  end function calc_Qprime_approx

  function calc_Qprime_approx2(dpsi, psi_soil, par_plant, par_env) result(Qprime)
    type(par_plant_type) :: par_plant
    type(par_env_type) :: par_env
    real(kind=dbl8) :: dpsi, psi_soil, Qprime, K

    K = scale_conductivity(par_plant%conductivity, par_env)
    Qprime = K * ((P(psi_soil, par_plant%psi50, par_plant%b) &
                 + P(psi_soil - dpsi, par_plant%psi50, par_plant%b)) / 2 &
                 - Pprime(psi_soil - dpsi, par_plant%psi50, par_plant%b) * dpsi / 2)
  end function calc_Qprime_approx2

  ! Derivative of sapflux wrt dpsi, dQ/ddpsi
  function calc_Qprime(dpsi, psi_soil, par_plant, par_env) result(Qprime)
    type(par_plant_type) :: par_plant
    type(par_env_type) :: par_env
    real(kind=dbl8) :: dpsi, psi_soil, Qprime

    if (par_env%gs_method == GS_APX) then
        Qprime = calc_Qprime_approx(dpsi, psi_soil, par_plant, par_env)
    else if (par_env%gs_method == GS_APX2) then
        Qprime = calc_Qprime_approx2(dpsi, psi_soil, par_plant, par_env)
    else if (par_env%gs_method == GS_IGF) then
        Qprime = calc_Qprime_analytical(dpsi, psi_soil, par_plant, par_env)
    ! else if (par_env%gs_method == GS_QNG) then
    !     Qprime = calc_Qprime_analytical(dpsi, psi_soil, par_plant, par_env)
    else
        write(*,*) "Unsupported gs_method specified"
        stop
    end if
  end function calc_Qprime

  function calc_dE_dgs_dif(par_env) result(dE_dgs)
    type(par_env_type) :: par_env
    real(kind=dbl8) :: dE_dgs, D

    D = dble(par_env%vpd) / dble(par_env%patm)
    dE_dgs = 1.6 * D
  end function calc_dE_dgs_dif

  function calc_dE_dgs_pm_from_gs(gs, par_plant, par_env) result(dE_dgs)
    type(par_plant_type) :: par_plant
    type(par_env_type) :: par_env
    real(kind=dbl8) :: gs, dE_dgs, ga

    ga = calc_g_aero(par_plant%h_canopy, dble(par_env%v_wind), par_plant%h_wind_measurement)
    dE_dgs = calc_dE_dgs_pm(gs, ga, par_env)
  end function calc_dE_dgs_pm_from_gs

  function calc_dE_dgs_pm_from_dpsi(dpsi, psi_soil, par_plant, par_env) result(dE_dgs)
    type(par_plant_type) :: par_plant
    type(par_env_type) :: par_env
    real(kind=dbl8) :: dpsi, psi_soil, dE_dgs, ga, Q, gs

    ga = calc_g_aero(par_plant%h_canopy, dble(par_env%v_wind), par_plant%h_wind_measurement)
    Q = calc_sapflux(dpsi, psi_soil, par_plant, par_env)
    gs = calc_gs_pm(Q, ga, par_env)
    dE_dgs = calc_dE_dgs_pm(gs, ga, par_env)
  end function calc_dE_dgs_pm_from_dpsi

  ! Derivative of E wrt gs
  function calc_dE_dgs_from_gs(gs, par_plant, par_env) result(dE_dgs)
    type(par_plant_type) :: par_plant
    type(par_env_type) :: par_env
    real(kind=dbl8) :: gs, dE_dgs

    if (par_env%et_method == ET_DIFFUSION) then
        dE_dgs = calc_dE_dgs_dif(par_env)
    else if (par_env%et_method == ET_PM) then
        dE_dgs = calc_dE_dgs_pm_from_gs(gs, par_plant, par_env)
    else
        write(*,*) "Unknown et_method:", par_env%et_method
        stop
    end if
  end function calc_dE_dgs_from_gs

  ! Derivative of E wrt gs
  function calc_dE_dgs_from_dpsi(dpsi, psi_soil, par_plant, par_env) result(dE_dgs)
    type(par_plant_type) :: par_plant
    type(par_env_type) :: par_env
    real(kind=dbl8) :: dpsi, psi_soil, dE_dgs

    if (par_env%et_method == ET_DIFFUSION) then
        dE_dgs = calc_dE_dgs_dif(par_env)
    else if (par_env%et_method == ET_PM) then
        dE_dgs = calc_dE_dgs_pm_from_dpsi(dpsi, psi_soil, par_plant, par_env)
    else
        write(*,*) "Unknown et_method:", par_env%et_method
        stop
    end if
  end function calc_dE_dgs_from_dpsi

  ! Derivative of gs wrt dpsi, dgs/ddpsi
  ! This version of the function avoids recomputation of gs when it is already known
  function calc_gsprime(dpsi, gs, psi_soil, par_plant, par_env) result(gsprime)
    type(par_plant_type) :: par_plant
    type(par_env_type) :: par_env
    real(kind=dbl8) :: dpsi, gs, psi_soil, gsprime, Qprime, Eprime

    Qprime = calc_Qprime(dpsi, psi_soil, par_plant, par_env)
    Eprime = calc_dE_dgs_from_gs(gs, par_plant, par_env)
    gsprime = Qprime / Eprime
  end function calc_gsprime

  ! Derivative of gs wrt dpsi, dgs/ddpsi
  ! This version is for use when gs is not known, and needs to be computed anyway
  function calc_gsprime_from_dpsi(dpsi, psi_soil, par_plant, par_env) result(gsprime)
    type(par_plant_type) :: par_plant
    type(par_env_type) :: par_env
    real(kind=dbl8) :: dpsi, psi_soil, gsprime, Qprime, Eprime

    Qprime = calc_Qprime(dpsi, psi_soil, par_plant, par_env)
    Eprime = calc_dE_dgs_from_dpsi(dpsi, psi_soil, par_plant, par_env)
    gsprime = Qprime / Eprime
  end function calc_gsprime_from_dpsi


end module md_phydro_transpiration  

