! module md_photosynth_phydro
!   use md_photosynth
!   use md_precision
!   use md_sofunutils
!   use md_sofunutils_plus
!   use md_phydro_main
!   implicit none

! !   ! list of methods to calculate gs
! !   integer (kind = int4), parameter :: GS_IGF = 0, GS_QNG = 1, GS_APX = 2, GS_APX2 = 3

! !   ! list of methods to calculate ET
! !   integer (kind = int4), parameter :: ET_DIFFUSION = 0, ET_PM = 1

! !   type par_control_type
! !     integer (kind = int4) :: gs_method
! !     integer (kind = int4) :: et_method
! !   end type par_control_type

! !   type par_plant_type
! !     real (kind = dbl8) :: conductivity        ! = ci/ca, leaf-internal to ambient CO2 partial pressure, ci/ca (unitless)
! !     real (kind = dbl8) :: psi50               ! leaf-internal CO2 partial pressure (Pa)
! !     real (kind = dbl8) :: b                   ! ci-limitation factor of light-limited assimilation (unitless)
! !     integer (kind = int4) :: gs_method = GS_IGF
! !   end type par_plant_type

! !   type par_cost_type
! !     real (kind = dbl8) :: alpha               ! = ci/ca, leaf-internal to ambient CO2 partial pressure, ci/ca (unitless)
! !     real (kind = dbl8) :: gamma               ! leaf-internal CO2 partial pressure (Pa)
! !   end type par_cost_type
  
! !   type par_env_type
! !     real (kind = dbl8) :: tc
! !     real (kind = dbl8) :: patm
! !     real (kind = dbl8) :: vpd
! !     real (kind = dbl8) :: viscosity_water
! !     real (kind = dbl8) :: density_water
! !   end type par_env_type

! !   type par_photosynth_type
! !     real (kind = dbl8) ::  kmm
! !     real (kind = dbl8) ::  gammastar
! !     real (kind = dbl8) ::  phi0
! !     real (kind = dbl8) ::  ca
! !     real (kind = dbl8) ::  delta
! !     real (kind = dbl8) ::  Iabs
! !     real (kind = dbl8) ::  patm
! !   end type par_photosynth_type

! !   type dpsi_bounds_type
! !     real (kind = dbl8) ::  exact
! !     real (kind = dbl8) ::  approx_O2
! !     real (kind = dbl8) ::  Iabs_bound
! !   end type dpsi_bounds_type

! !   type dfdx_type
! !     real (kind = dbl8) ::  dPdx
! !     real (kind = dbl8) ::  J
! !     real (kind = dbl8) ::  djmax_dJ
! !     real (kind = dbl8) ::  dJ_dchi
! !   end type dfdx_type

! !   type outtype_phydro
! !     real (kind = dbl8) :: a
! !     real (kind = dbl8) :: e
! !     real (kind = dbl8) :: gs
! !     real (kind = dbl8) :: ci                  ! leaf-internal partial pressure, (Pa)
! !     real (kind = dbl8) :: chi                 ! = ci/ca, leaf-internal to ambient CO2 partial pressure, ci/ca (unitless)
! !     real (kind = dbl8) :: vcmax
! !     real (kind = dbl8) :: jmax
! !     real (kind = dbl8) :: dpsi
! !     real (kind = dbl8) :: psi_l
! !     real (kind = dbl8) :: gammastar           ! temperature-dependent photorespiratory compensation point (Pa)
! !     real (kind = dbl8) :: kmm                 ! Michaelis-Menten coefficient (Pa)
! !     real (kind = dbl8) :: ca                  ! leaf-external (ambient) partial pressure, (Pa)
! !     real (kind = dbl8) :: vcmax25
! !     real (kind = dbl8) :: nfnct;
! !     real (kind = dbl8) :: niter;
! !     real (kind = dbl8) :: mc;
! !     real (kind = dbl8) :: mj;
! !   end type outtype_phydro

! !   interface par_env_type
! !     module procedure :: init_env
! !   end interface

! !   interface par_photosynth_type
! !     module procedure :: init_photo
! !   end interface

! !   contains

! !   function init_env(tc, patm, vpd) result(ppar)
! !     ! use md_photosynth, only: calc_viscosity_h2o, calc_density_h2o
! !     real (kind = dbl8), intent(in) :: tc, patm, vpd        ! = ci/ca, leaf-internal to ambient CO2 partial pressure, ci/ca (unitless)
! !     type(par_env_type) :: ppar
! !     ppar%tc = tc
! !     ppar%vpd = vpd
! !     ppar%patm = patm
! !     ppar%viscosity_water = calc_viscosity_h2o(real(tc), real(patm))
! !     ppar%density_water = calc_density_h2o(real(tc), real(patm))
! !   end

! !   function init_photo(tc, patm, kphio, co2, ppfd, fapar, rdark) result(ppar)
! !     real (kind = dbl8), intent(in) :: tc, patm, kphio, co2, ppfd, fapar, rdark
! !     type(par_photosynth_type) :: ppar
! !     ppar%kmm = calc_kmm(real(tc), real(patm))
! !     ppar%gammastar = calc_gammastar(real(tc), real(patm))
! !     ppar%phi0 = kphio*calc_ftemp_kphio(real(tc), c4 = .false.)
! !     ppar%Iabs = ppfd * fapar
! !     ppar%ca = co2 * patm * 1e-6
! !     ppar%patm = patm
! !     ppar%delta = rdark
! !   end
 


! !   !!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! !   !!! Phydro analytical solver (acclimating)
! !   !!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! !   function calc_J(gs, x, par_photosynth) result(J)
! !     real (kind = dbl8), intent(in) :: gs, x
! !     type(par_photosynth_type), intent(in) :: par_photosynth
! !     real (kind = dbl8) :: g, K, ca, d, J
! !     g = par_photosynth%gammastar / par_photosynth%ca
! !     k = par_photosynth%kmm / par_photosynth%ca
! !     ca = par_photosynth%ca / par_photosynth%patm*1e6
! !     d = par_photosynth%delta
! !     J = 4*gs*ca*(1-x)*(x+2*g)/(x*(1-d)-(g+d*k))
! !   end


! !   function calc_jmax_from_J(J, par_photosynth) result(jmax)
! !     real (kind = dbl8), intent(in) :: J
! !     type(par_photosynth_type), intent(in) :: par_photosynth
! !     real (kind = dbl8) :: pp, pj, jmax
! !     pp = 4*par_photosynth%phi0 * par_photosynth%Iabs;
! !     pj = pp/J;
! !     jmax = pp/sqrt(pj*pj-1);
! !   end


! !   function calc_djmax_dJ(J, par_photosynth) result(djdj)
! !     real (kind = dbl8), intent(in) :: J
! !     type(par_photosynth_type), intent(in) :: par_photosynth
! !     real (kind = dbl8) :: pp, sq, psq, djdj
! !     pp = 4*par_photosynth%phi0 * par_photosynth%Iabs
! !     sq = sqrt(pp*pp-J*J)
! !     psq = pp/sq
! !     djdj = psq*psq*psq
! !   end


! !   function calc_dJ_dchi(gs, x, par_photosynth) result(djdx)
! !     real (kind = dbl8), intent(in) :: gs, x
! !     type(par_photosynth_type), intent(in) :: par_photosynth
! !     real (kind = dbl8) :: g, K, ca, d, djdx, d1
! !     g = par_photosynth%gammastar / par_photosynth%ca
! !     k = par_photosynth%kmm / par_photosynth%ca
! !     ca = par_photosynth%ca / par_photosynth%patm*1e6
! !     d = par_photosynth%delta
! !     ! gs*ca * ((d*(2*g*(k + 1) + k*(2*x - 1) + x^2) + 2*g^2 + g*(2*x - 3) - x^2)/(d*(k + x) + g - x)^2)
! !     d1 = d*(k + x) + g - x;
! !     djdx = 4*gs*ca * ((d*(2*g*(k + 1) + k*(2*x - 1) + x*x) - ((x-g)*(x-g)+3*g*(1-g)))/(d1*d1));
! !     ! gs*ca*(3*(g-1)*g/(g-x)^2 - 1)
! !   end


! !   function calc_dJ_ddpsi(gsprime, x, par_photosynth) result(djdp)
! !     real (kind = dbl8), intent(in) :: gsprime, x
! !     type(par_photosynth_type), intent(in) :: par_photosynth
! !     real (kind = dbl8) :: g, K, ca, d, djdp
! !     g = par_photosynth%gammastar / par_photosynth%ca
! !     k = par_photosynth%kmm / par_photosynth%ca
! !     ca = par_photosynth%ca / par_photosynth%patm*1e6
! !     d = par_photosynth%delta
! !     djdp = 4*gsprime*ca*(1-x)*(x+2*g)/(x*(1-d)-(g+d*k))
! !   end


! !   function calc_x_from_dpsi(dpsi, psi_soil, par_plant, par_env, par_photosynth, par_cost) result(x)
! !     real (kind = dbl8), intent(in) :: dpsi, psi_soil
! !     type(par_plant_type), intent(in) :: par_plant
! !     type(par_env_type), intent(in) :: par_env
! !     type(par_photosynth_type), intent(in) :: par_photosynth
! !     type(par_cost_type), intent(in) :: par_cost

! !     real (kind = dbl8) gstar, Km, ca, br, y, gsprime, ca2, x
! !     gstar = par_photosynth%gammastar/par_photosynth%patm*1e6
! !     Km = par_photosynth%kmm/par_photosynth%patm*1e6
! !     ca = par_photosynth%ca/par_photosynth%patm*1e6
! !     br = par_photosynth%delta
! !     y = par_cost%gamma
    
! !     gsprime = calc_gsprime(dpsi, psi_soil, par_plant, par_env)
  
! !     ca2 = ca*ca;
! !     x = (-2*ca*dpsi*(gstar + br*Km)*y + &
! !       ca2*((3 - 2*br)*gstar + br*Km)*gsprime + &
! !       -sqrt(2.0D+00)*sqrt( &
! !         ca2*dpsi*((-3 + 2*br)*gstar - br*Km)*((-1 + br)*ca + gstar + &
! !                                                   br*Km)*y* &
! !           (-2*dpsi*y + (ca + 2*gstar)* &
! !               gsprime)))/ &
! !       (ca2*(2*(-1 + br)*dpsi*y + ((3 - 2*br)*gstar + br*Km)* &
! !               gsprime))
    
! !     if (x < (gstar + br*Km)/(ca - br*ca)) x = (gstar + br*Km)/(ca - br*ca)+1e-12
! !   end
  

! !   function calc_delta_from_dpsi(dpsi, psi_soil, par_plant, par_env, par_photosynth, par_cost) result(delt)
! !     real (kind = dbl8), intent(in) :: dpsi, psi_soil
! !     type(par_plant_type), intent(in) :: par_plant
! !     type(par_env_type), intent(in) :: par_env
! !     type(par_photosynth_type), intent(in) :: par_photosynth
! !     type(par_cost_type), intent(in) :: par_cost
    
! !     real (kind = dbl8) gstar, Km, ca, br, y
! !     real (kind = dbl8) :: gsprime, delt
    
! !     gstar = par_photosynth%gammastar/par_photosynth%patm*1e6
! !     Km = par_photosynth%kmm/par_photosynth%patm*1e6
! !     ca = par_photosynth%ca/par_photosynth%patm*1e6
! !     br = par_photosynth%delta
! !     y = par_cost%gamma

! !     gsprime = calc_gsprime(dpsi, psi_soil, par_plant, par_env)
! !     delt = (-2*dpsi*y + (ca + 2*gstar)*gsprime);
! !   end


! !   function calc_dpsi_bound(psi_soil, par_plant, par_env, par_photosynth, par_cost) result(bounds)
! !     real (kind = dbl8), intent(in) :: psi_soil
! !     type(par_plant_type), intent(in) :: par_plant
! !     type(par_env_type), intent(in) :: par_env
! !     type(par_photosynth_type), intent(in) :: par_photosynth
! !     type(par_cost_type), intent(in) :: par_cost

! !     type(dpsi_bounds_type) :: bounds 

! !     real (kind = dbl8) :: gstar, ca, y, K, Pox, Ppox, Pppox
! !     real (kind = dbl8) :: a,b,c,del
! !     real (kind = dbl8) :: ex, appo2, iabsb, use_bound

! !     gstar = par_photosynth%gammastar/par_photosynth%patm*1e6
! !     ca = par_photosynth%ca/par_photosynth%patm*1e6
! !     y = par_cost%gamma

! !     K = scale_conductivity(par_plant%conductivity, par_env)/(1.6*par_env%vpd/par_env%patm);
 
! !     Pox = P(psi_soil, par_plant%psi50, par_plant%b);
! !     Ppox = Pprime(psi_soil, par_plant%psi50, par_plant%b);
! !     Pppox = Pprimeprime(psi_soil, par_plant%psi50, par_plant%b);
    
! !     a = (ca + 2*gstar)*K*Pppox*4/8;
! !     b = -(2*y + (ca + 2*gstar)*K*Ppox);
! !     c = (ca + 2*gstar)*K*Pox;
! !     del = b*b-4*a*c;

! !     appo2 = (-b-sqrt(del))/(2*a)
! !     ex = zero(0.0D+00, 10.0D+00, f2, 1D-6)

! !     use_bound = ex
    
! !     !//# cat(psi_soil, ":", exact, " ", approx_O2, " ", use_bound, "\n");
! !     iabsb = zero(use_bound*0.001, use_bound*0.99, f1, 1D-6);
    
! !     !//# dpsi=seq(exact*0.001,exact*0.99, length.out=200);
! !     !//# plot(y=sapply(X = dpsi, FUN = f1), x=dpsi, type="l");
    
! !     bounds = dpsi_bounds_type(ex, appo2, iabsb)

! !     contains

! !     function f2(dpsi) result(gg)
! !       real(kind = dbl8), intent(in) :: dpsi
! !       real(kind = dbl8) :: gg
! !       gg = (-2*dpsi*y + (ca + 2*gstar)*calc_gsprime(dpsi, psi_soil, par_plant, par_env))
! !     end

! !     function f1(dpsi) result(J)
! !       real(kind = dbl8), intent(in) :: dpsi
! !       real(kind = dbl8) :: J, gs, x
! !       gs = calc_gs(dpsi, psi_soil, par_plant, par_env);
! !       x = calc_x_from_dpsi(dpsi, psi_soil, par_plant, par_env, par_photosynth, par_cost)
! !       J = calc_J(gs, x, par_photosynth)-4*par_photosynth%phi0*par_photosynth%Iabs;
! !     end

! !   end


! !   function dFdx(dpsi, psi_soil, par_plant, par_env, par_photosynth, par_cost) result(res)
! !     real (kind = dbl8), intent(in) :: dpsi, psi_soil
! !     type(par_plant_type), intent(in) :: par_plant
! !     type(par_env_type), intent(in) :: par_env
! !     type(par_photosynth_type), intent(in) :: par_photosynth
! !     type(par_cost_type), intent(in) :: par_cost

! !     real (kind = dbl8) :: gs, gsprime, X, J, ca, g, djmax_dJ, dJ_dchi, dP_dx
! !     type(dfdx_type) :: res

! !     gs = calc_gs(dpsi, psi_soil, par_plant, par_env)
! !     gsprime = calc_gsprime(dpsi, psi_soil, par_plant, par_env)

! !     X =  calc_x_from_dpsi(dpsi, psi_soil, par_plant, par_env, par_photosynth, par_cost)

! !     J = calc_J(gs, X, par_photosynth)

! !     ca = par_photosynth%ca / par_photosynth%patm*1e6
! !     g = par_photosynth%gammastar / par_photosynth%ca

! !     djmax_dJ = calc_djmax_dJ(J, par_photosynth)
! !     dJ_dchi  = calc_dJ_dchi(gs, X, par_photosynth)

! !     dP_dx = -gs*ca - par_cost%alpha * djmax_dJ * dJ_dchi

! !     res = dfdx_type(dP_dx, J, djmax_dJ, dJ_dchi)
! !   end


! !   ! Questions for Beni: 
! !   ! 1) no input fapar?
! !   ! 2) How to calc vcmax25?
! !   function phydro(kphio, par_cost, psi_soil, ppfd, fapar, co2, tc, vpd, patm, par_plant, rdark) result(out_phydro)
! !     !//////////////////////////////////////////////////////////////////
! !     ! Implements the P-hydro model, providing predictions for ci, Vcmax, and 
! !     ! Jmax, A, dpsi, and gs. 
! !     !------------------------------------------------------------------
! !     ! arguments
! !     real (kind = dbl8), intent(in) :: kphio        ! apparent quantum yield efficiency       
! !     !real, intent(in) :: beta         ! parameter for the unit cost ratio (corresponding to beta in Prentice et al., 2014)    
! !     ! real, intent(in) :: fapar        ! fraction of absorbed photosynthetically active radiation (unitless) 
! !     real (kind = dbl8), intent(in) :: psi_soil         ! soil water potential (Mpa)
! !     real (kind = dbl8), intent(in) :: ppfd         ! photosynthetic photon flux density (mol m-2 s-1), relevant for acclimated response
! !     real (kind = dbl8), intent(in) :: co2          ! atmospheric CO2 concentration (ppm), relevant for acclimated response
! !     real (kind = dbl8), intent(in) :: tc           ! air temperature (deg C), relevant for acclimated response
! !     real (kind = dbl8), intent(in) :: vpd          ! vapor pressure (Pa), relevant for acclimated response
! !     real (kind = dbl8), intent(in) :: patm         ! atmospheric pressure (Pa), relevant for acclimated response
! !     !logical, intent(in) :: c4        ! whether or not C4 photosynthesis pathway is followed. If .false., it's C3.
! !     type(pxar_plant_type), intent(in) :: par_plant
! !     type(par_cost_type), intent(in) :: par_cost
! !     real (kind = dbl8), intent(in) :: fapar         ! fraction of PAR absorbed
! !     real (kind = dbl8), intent(in) :: rdark         ! dark respiration rate (-), as a fraction of Vcmax
 
! !     ! function return value
! !     type(outtype_phydro) :: out_phydro

! !     type(par_env_type) :: par_env
! !     type(par_photosynth_type) :: par_photosynth

! !     type(dpsi_bounds_type) :: bounds

! !     real (kind = dbl8) :: dpsi_opt, x, gs, J, jmax, vcmax, a

! !     par_photosynth = par_photosynth_type(tc, patm, kphio, co2, ppfd, fapar, rdark)
! !     par_env = par_env_type(tc, patm, vpd);
  
! !     bounds   = calc_dpsi_bound(psi_soil, par_plant, par_env, par_photosynth, par_cost)
! !     dpsi_opt = zero(bounds%Iabs_bound*0.001, bounds%Iabs_bound*0.999, fn, 1.0D-6)
! !     x        = calc_x_from_dpsi(dpsi_opt, psi_soil, par_plant, par_env, par_photosynth, par_cost)
! !     gs       = calc_gs(dpsi_opt, psi_soil, par_plant, par_env)
! !     J        = calc_J(gs, x, par_photosynth)
! !     jmax     = calc_jmax_from_J(J, par_photosynth)
! !     vcmax    = (J/4.0)*(x*par_photosynth%ca + par_photosynth%kmm)/(x*par_photosynth%ca + 2*par_photosynth%gammastar)
! !     a        = gs*(par_photosynth%ca/par_photosynth%patm*1e6)*(1-x)
  
! !     out_phydro%a         = a
! !     out_phydro%e         = 1.6*gs*vpd/par_env%patm
! !     out_phydro%gs        = gs
! !     out_phydro%ci        = x*par_photosynth%ca        ! leaf-internal partial pressure, (Pa)
! !     out_phydro%chi       = x        ! = ci/ca, leaf-internal to ambient CO2 partial pressure, ci/ca (unitless)
! !     out_phydro%vcmax     = vcmax
! !     out_phydro%jmax      = jmax
! !     out_phydro%dpsi      = dpsi_opt
! !     out_phydro%psi_l     = psi_soil - dpsi_opt
! !     out_phydro%gammastar = par_photosynth%gammastar        ! temperature-dependent photorespiratory compensation point (Pa)
! !     out_phydro%kmm       = par_photosynth%kmm         ! Michaelis-Menten coefficient (Pa)
! !     out_phydro%ca        = par_photosynth%ca        ! leaf-external (ambient) partial pressure, (Pa)
! !     out_phydro%vcmax25   = -999
! !     out_phydro%nfnct     = -999
! !     out_phydro%niter     = -999
! !     out_phydro%mc        = (x*par_photosynth%ca - par_photosynth%gammastar) / (x*par_photosynth%ca + par_photosynth%kmm);
! !     out_phydro%mj        = (x*par_photosynth%ca - par_photosynth%gammastar) / (x*par_photosynth%ca + 2*par_photosynth%gammastar);
  
! !     contains 

! !     function fn(dpsi) result(dpdx)
! !       real (kind = dbl8), intent(in) :: dpsi
! !       real (kind = dbl8) :: dpdx
! !       type(dfdx_type) :: dfdx1
! !       dfdx1 = dFdx(dpsi, psi_soil, par_plant, par_env, par_photosynth, par_cost)
! !       dpdx = dfdx1%dPdx
! !     end

! !   end function phydro

! end module md_photosynth_phydro


! ============================
!     Test gs
! ============================ 
subroutine test_gs
  use md_precision
  use md_sofunutils_plus
  use md_sofunutils
  use md_phydro_transpiration

  implicit none

  real (kind = dbl8) :: test
  real(kind = dbl8) :: tc, patm, vpd
  real(kind = dbl8) :: g, psi_s, Q
  type(par_plant_type) :: par_plant
  type(par_env_type) :: par_env
  integer :: N
  integer :: i


  tc   = 25
  patm = calc_patm(0.0)
  vpd  = 1000
	
  call create_par_env(par_env, tc, patm, vpd, 1000.0d0/2, 3.0d0)

  par_plant = par_plant_type(3.0D-17, -2.0D00, 2.0D00)
  par_plant%gs_method = GS_IGF


  N = 10
  do i = 0, N-1, 1
    psi_s = -6.0 + i*(6.0)/(N-1);
    Q = calc_sapflux(1.0d0, psi_s, par_plant, par_env)
    g = calc_gs_from_Q(Q, psi_s, par_plant, par_env);
    print *, psi_s , g;
  end do

  if (abs(g - 0.1116382) < 1e-5) then; print *, "gs_test PASS"
  else; print *, "gs_test FAIL"
  end if

  test = calc_patm(0.0)
  print *, "patm(0) =", test


end



! ============================
!     Test dFdx
! ============================ 
subroutine test_dfdx
  use md_precision
  use md_sofunutils, only: calc_patm
  use md_sofunutils_plus
  use md_phydro_solver

  implicit none

  real(kind = dbl8) :: psi_s, tc, patm, vpd, kphio, co2, ppfd, fapar, rdark
  type(par_plant_type) :: par_plant
  type(par_env_type) :: par_env
  type(par_photosynth_type) :: par_photosynth, ph_test
  type(par_cost_type) :: par_cost
  integer :: N
  integer :: i
  type(dpsi_bounds_type) :: b1
  type(dfdx_type) :: b

  real(dbl8) :: q1, g1, g1p, x1, J1

  tc   = 25
  vpd = 1000
  kphio = 0.087
  co2 = 400
  ppfd = 1000
  fapar = 1
  rdark = 0.02
  patm = calc_patm(0.0)
	
  par_plant = par_plant_type(3.0D-17, -2.0D00, 2.0D00)
  par_plant%gs_method = GS_IGF

  call create_par_env(par_env, tc, patm, vpd, ppfd/2, 3.0d0)
  call create_par_photosynth(par_photosynth, tc, patm, kphio, co2, ppfd, fapar, rdark, &
                             tc, 25.0d0, FV_kumarathunge19, FR_heskel16, FB_atkin15)
  par_cost = par_cost_type(0.1, 1)

  call create_par_photosynth(ph_test, 28.d0, patm, kphio, co2, ppfd, fapar, rdark, &
                             19.d0, 25.0d0, FV_kumarathunge19, FR_heskel16, FB_atkin15)
  call print_par_photosynth(ph_test)

  q1 = calc_sapflux(1.0d0, -1.0d0, par_plant, par_env);
  g1 = calc_gs_from_Q(q1, -1.0d0, par_plant, par_env);
  g1p = calc_gsprime(1.0d0, g1, -1.0d0, par_plant, par_env);
  x1 = calc_x_from_dpsi(1.0d0, g1p, par_photosynth, par_cost);
  J1 = calc_J(g1, x1, par_photosynth)-4.0d0*par_photosynth%phi0*par_photosynth%Iabs;

  print *, "dFdx f2: ", q1, g1, g1p, x1, J1

  N = 10
  do i = 0, N-1, 1
    psi_s = -6.0 + i*(6.0)/(N-1);
    b1 = calc_dpsi_bound(dble(psi_s), par_plant, par_env, par_photosynth, par_cost)
    b = dFdx(b1%Iabs_bound*0.5, dble(psi_s), par_plant, par_env, par_photosynth, par_cost)
    print *, psi_s, b1%exact, b1%Iabs_bound, b%dJ_dchi, b%djmax_dJ, b%dPdx, b%J
  end do

  if ((abs(b1%exact - 3.486955) < 1e-4) .and. &
      (abs(b1%Iabs_bound - 1.846519) < 1e-4) .and. &
      (abs(b%dJ_dchi - -328.536613)/328 < 1e-4) .and. &
      (abs(b%djmax_dJ - 1.182968) < 1e-4) .and. &
      (abs(b%dPdx - -2.701680) < 1e-4) .and. &
      (abs(b%J - 78.110340)/78 < 1e-4)) & 
      then; print *, "dFdx_test PASS"
  else; print *, "dFdx_test FAIL"
  end if 

end


! ============================
!     Test Phydro
! ============================ 
subroutine test_phydro
  use md_precision
  use md_sofunutils_plus
  use md_sofunutils
  use md_phydro_solver
  use md_phydro_main

  implicit none
  real(kind = dbl8) :: psi_soil, tc, patm, vpd, kphio, co2, ppfd, fapar, rdark, elv 
  real(kind = dbl8) :: dpsi1, grad1
  type(par_plant_type) :: par_plant
  type(par_photosynth_type) :: par_photosynth
  type(par_env_type) :: par_env
  type(par_cost_type) :: par_cost
  integer :: N
  integer :: i
  type(phydro_result_type) :: res, res_ref
  type(par_control_type) :: options
  logical :: err = .false.

  elv = 0
  tc   = 25
  vpd = 1000
  kphio = 0.087
  co2 = 400
  ppfd = 300
  fapar = 0.7
  rdark = 0.02
  patm = calc_patm(0.0)

  par_plant = par_plant_type(3.0D-17, -2.0D00, 2.0D00)

  options%gs_method = GS_IGF

  par_cost = par_cost_type(0.1, 1)

  open(unit=1, file='/home/jjoshi/codes/phydro/tests/test_data/psi.tsv')

  print *, "Soil moist response: GS_IGF PM_DIFFUSION"
  N = 20
  do i = 0, N-1, 1
    psi_soil = -6.0 + i*(6.0)/(N-1);
    res = phydro_analytical(tc, tc, ppfd, ppfd/2, vpd, co2, elv, fapar, kphio, psi_soil, rdark, 3.0d0, par_plant, par_cost, options)
    read(1,*) psi_soil, res_ref%jmax, res_ref%dpsi, res_ref%gs, res_ref%a, res_ref%ci, res_ref%chi, res_ref%vcmax
    write (*,10) psi_soil, res%jmax, res%dpsi, res%gs, res%a, res%ci, res%chi, res%vcmax
    if (abs(res%jmax - res_ref%jmax) > 5e-4 .or. &
        abs(res%jmax - res_ref%jmax) > 5e-4 .or. &
        abs(res%dpsi - res_ref%dpsi) > 5e-4 .or. &
        abs(res%gs - res_ref%gs) > 5e-4 .or. &
        abs(res%a - res_ref%a) > 5e-4 .or. &
        abs(res%ci - res_ref%ci) > 5e-4 .or. &
        abs(res%chi - res_ref%chi) > 5e-4 .or. &
        abs(res%vcmax - res_ref%vcmax) > 5e-4) then
          write (*,11) psi_soil, res_ref%jmax, res_ref%dpsi, res_ref%gs, res_ref%a, res_ref%ci, res_ref%chi, res_ref%vcmax
          11  format(8f12.6)
          print *, "Phydro ET_DIFFUSION Test failed at psis = ", psi_soil
       err = .true.
    endif
    10  format(8f12.6)
  end do

  close(unit = 1)


  options%et_method = ET_PM

  open(unit=2, file='/home/jjoshi/codes/phydro/tests/test_data/psi_pm.tsv')
  print *, "Soil moist response: GS_IGF PM_ET"
  N = 20
  do i = 0, N-1, 1
    psi_soil = -6.0 + i*(6.0)/(N-1);
    res = phydro_analytical(tc, tc, ppfd, ppfd/2, vpd, co2, elv, fapar, kphio, psi_soil, rdark, 3.0d0, par_plant, par_cost, options)
    read(2,*) psi_soil, res_ref%jmax, res_ref%dpsi, res_ref%gs, res_ref%a, res_ref%ci, res_ref%chi, res_ref%vcmax
    write (*,20) psi_soil, res%jmax, res%dpsi, res%gs, res%a, res%ci, res%chi, res%vcmax
    if (abs(res%jmax - res_ref%jmax) > 5e-4 .or. &
        abs(res%jmax - res_ref%jmax) > 5e-4 .or. &
        abs(res%dpsi - res_ref%dpsi) > 5e-4 .or. &
        abs(res%gs - res_ref%gs) > 5e-4 .or. &
        abs(res%a - res_ref%a) > 5e-4 .or. &
        abs(res%ci - res_ref%ci) > 5e-4 .or. &
        abs(res%chi - res_ref%chi) > 5e-4 .or. &
        abs(res%vcmax - res_ref%vcmax) > 5e-4) then
        write (*,21) psi_soil, res_ref%jmax, res_ref%dpsi, res_ref%gs, res_ref%a, res_ref%ci, res_ref%chi, res_ref%vcmax
        21  format(8f12.6)
      print *, "Phydro ET_PM Test failed at psis = ", psi_soil
      err = .true.
    endif
   20  format(8f12.6)
  end do

  !---------------------------------------------------------
  !  Instantaneous version
  !---------------------------------------------------------

  tc = 20
  vpd = 810.6
  co2 = 400
  ppfd = 1200
  patm = calc_patm(0.0)
  kphio = 0.087
  fapar = 0.99
  rdark = 0.02
  psi_soil = -0.4137931

  par_cost = par_cost_type(0.118514, 1.227068)
  par_plant = par_plant_type(7.457324e-17, -1.039539, 1)

  options%et_method = ET_DIFFUSION
  options%gs_method = GS_IGF

  call create_par_photosynth(par_photosynth, tc, patm, kphio, co2, ppfd, fapar, rdark, &
                            tc, tc, FV_kumarathunge19, FR_heskel16, FB_atkin15)
  call create_par_env(par_env, tc, patm, vpd, ppfd/2, 3.0d0)

  do i = 0, 49, 1
    dpsi1 = 2.0/49.0*i;
    grad1 = calc_dP_ddpsi(dpsi1, 55.6279401d0, 117.0184518d0, psi_soil, par_plant, par_env, par_photosynth, par_cost);
    print *, dpsi1, grad1
  end do

  open(unit=3, file='/home/jjoshi/codes/phydro/tests/test_data/psi_inst.tsv')
  print *, "Soil moist response: Inst   GS_IGF   ET_DIFFUSION"
  N = 20
  do i = 0, N-1, 1
    psi_soil = -6.0 + i*(6.0)/(N-1);
    res = phydro_instantaneous_analytical(55.6279401d0, 117.0184518d0, tc, tc, ppfd, ppfd/2, vpd, co2, elv, &
                                          fapar, kphio, psi_soil, rdark, 3.0d0, par_plant, par_cost, options)
    read(3,*) psi_soil, res_ref%jmax, res_ref%dpsi, res_ref%gs, res_ref%a, res_ref%ci, res_ref%chi, res_ref%vcmax
    write (*,30) psi_soil, res%jmax, res%dpsi, res%gs, res%a, res%ci, res%chi, res%vcmax
    if (abs(res%jmax - res_ref%jmax) > 5e-4 .or. &
        abs(res%jmax - res_ref%jmax) > 5e-4 .or. &
        abs(res%dpsi - res_ref%dpsi) > 5e-4 .or. &
        abs(res%gs - res_ref%gs) > 5e-4 .or. &
        abs(res%a - res_ref%a) > 5e-4 .or. &
        abs(res%ci - res_ref%ci) > 5e-4 .or. &
        abs(res%chi - res_ref%chi) > 5e-4 .or. &
        abs(res%vcmax - res_ref%vcmax) > 5e-4) then
        write (*,31) psi_soil, res_ref%jmax, res_ref%dpsi, res_ref%gs, res_ref%a, res_ref%ci, res_ref%chi, res_ref%vcmax
        31  format(8f12.6)
      print *, "Phydro Inst Test failed at psis = ", psi_soil
      err = .true.
    endif
   30  format(8f12.6)
  end do


  print *, "--------------------"
  if (err) then 
    print *, "Some tests failed!"
  else 
    print *, "All tests passed!!"
  endif 
end  



! Test program for Environment
subroutine test_par_env
  use md_phydro_env
  implicit none

  type(par_env_type) :: my_env
  type(par_env_type) :: my_env2

  ! Create a ParEnv object with all parameters
  call create_par_env(my_env, 25.0d0, 101325.0d0, 2000.0d0, 300.0d0, 5.0d0)

  ! Create a ParEnv object with default v_wind
  call create_par_env(my_env2, 20.0d0, 101325.0d0, 1000.0d0, 200.0d0, 3.0d0)

  ! Print the ParEnv object
  call print_par_env(my_env)
  call print_par_env(my_env2)

end subroutine test_par_env


! Test program for Environment
subroutine test_dpsi_calc
  use md_phydro_transpiration
  implicit none

  ! real(kind=dbl8) :: tc = 25.0
  ! real(kind=dbl8) :: vpd = 1000.0
  ! real(kind=dbl8) :: kphio = 0.087
  ! real(kind=dbl8) :: co2 = 400.0
  ! real(kind=dbl8) :: ppfd = 1000.0
  ! real(kind=dbl8) :: fapar = 1.0
  ! real(kind=dbl8) :: rdark = 0.02
  real(kind=dbl8) :: psi_s = -1.0
  integer :: N = 10
  integer :: i
  real(kind=dbl8) :: dpsi, trans, dpsi_back, trans_max

  ! Declare ParPlant, ParEnv, and ParPhotosynth as appropriate

  type(par_plant_type) :: P1;
  type(par_env_type) :: E;

  P1%conductivity = 3e-16
  P1%psi50 = -2
  P1%b = 2

  call create_par_env(E, 25.0d0, 101325.0d0, 2000.0d0, 300.0d0, 5.0d0)

  do i = 0, N - 1
      dpsi = 0.0 + i * (10.0) / (N - 1)
      trans = calc_sapflux(dpsi, psi_s, P1, E)
      dpsi_back = calc_dpsi_from_sapflux(trans, psi_s, P1, E)
      write(*, '(F10.6, F20.10, F20.10)') dpsi, trans, dpsi_back

      if (abs(dpsi - dpsi_back) > 2.0e-6) then
          write(*, *) "Error: |dpsi - dpsi_back| > 2.0e-6"
          return
      end if
  end do

  N = 50
  trans_max = calc_max_sapflux(psi_s, P1, E)
  write(*, *) "Max trans = ", trans_max
  do i = 0, N - 1
      trans = trans_max / 2.0 + i * (trans_max) / (N - 1)
      dpsi_back = calc_dpsi_from_sapflux(trans, psi_s, P1, E)
      write(*, '(F20.6, F20.6)') trans, dpsi_back
  end do
end subroutine test_dpsi_calc


subroutine test_jmax_temp_response
  use md_phydro_photosynthesis
  implicit none

  real(kind = dbl8), parameter :: err_tol = 1.0E-6
  integer :: i
  logical :: is_error = .false.
  real(kind = dbl8) :: x, ref, err_result

  ! Sequence of temperature values
  real(kind = dbl8), dimension(46) :: temperature_seq

  ! Expected results for calc_ftemp_inst_jmax with FV_kumarathunge19
  real(kind = dbl8), dimension(46) :: expected_kumarathunge
  data expected_kumarathunge / 0.16224749762063, 0.173636341666431, 0.185731247018254, 0.19856989574356, &
                               0.212191620359412, 0.226637452541204, 0.241950169605546, &
                               0.258174337337494, 0.275356347157117, 0.293544444822553, 0.312788746763373, &
                               0.333141238616272, 0.354655748442466, 0.377387884237216, &
                               0.401394921421314, 0.426735620665162, 0.453469949152501, 0.481658668605213, &
                               0.511362740234979, 0.542642479202111, 0.575556367816386, &
                               0.610159405985626, 0.646500837405835, 0.684621038620159, 0.724547293377837, &
                               0.766288095388894, 0.809825529038732, 0.855105174035919, &
                               0.902022877428985, 0.950407658283492, 1.0, 1.05042491711027, 1.10115957532813, &
                               1.15149606748754, 1.20050142672755, 1.2469793282723, &
                               1.28944133406019, 1.32609979559137, 1.35489880528583, 1.37360195980551, &
                               1.37995303612314, 1.37191431817576, 1.34796519490302, &
                               1.30741376997327, 1.25064729378096, 1.17923975686556 /

  ! Expected results for calc_ftemp_inst_jmax with FV_kattge07
  real(kind = dbl8), dimension(46) :: expected_kattge
  data expected_kattge / 0.108165296210614, 0.1175411330303, 0.127651024795746, 0.138546040580533, &
                           0.150280321556271, 0.162911228097918, 0.176499488675189, &
                           0.19110934874336, 0.20680871702301, 0.2236693053978, 0.241766757029895, &
                           0.261180755011782, 0.281995100691034, 0.30429774637993, &
                           0.328180761031501, 0.353740199003026, 0.381075830401281, 0.410290675603626, &
                           0.441490264927, 0.474781515212076, 0.510271075968114, &
                           0.548062945862371, 0.588255092510554, 0.630934721404157, 0.676171730670059, &
                           0.724009756496315, 0.774454063444923, 0.82745537775353, &
                           0.882888630174217, 0.940525526345156, 1.0, 1.06076609384524, 1.12204890000608, &
                           1.18279119576327, 1.24160166415102, 1.29671527092002, &
                           1.3459821620838, 1.38690701584358, 1.41676324547863, 1.43280134465593, 1.4325533122752, &
                           1.41420402895972, 1.37696211598528, &
                           1.32133373738511, 1.24920421515943, 1.16367598004045 /

  ! Expected results for calc_ftemp_inst_jmax with FV_leuning02
  real(kind = dbl8), dimension(46) :: expected_leuning
  data expected_leuning / 0.120375999968243, 0.130894617257924, 0.142242095547883, 0.154475762525289, &
                             0.167656087918668, 0.181846711992658, 0.197114428717115, &
                             0.213529106725458, 0.231163525951592, 0.250093101138906, 0.270395454897424, &
                             0.292149792245117, 0.315436015172373, 0.340333499283094, &
                             0.366919434633578, 0.395266609349475, 0.425440487691548, 0.457495404927633, &
                             0.491469671831667, 0.527379355975122, 0.565210492197857, &
                             0.604909481860986, 0.646371486311908, 0.68942672777097, 0.733824811005138, &
                             0.779217507449454, 0.825140935113514, 0.870998744771333, &
                             0.916048772557816, 0.959396560997767, 1.0, 1.03668977973091, 1.06820994072775, &
                             1.09328108237436, 1.11068543953201, 1.11936819194594, &
                             1.1185438676821, 1.10779211771609, 1.08712541568215, 1.05701387456219, &
                             1.01835951405898, 0.972422293501319, 0.920710029094153, &
                             0.864850840297628, 0.806468217168184, 0.747075353811709 /

  ! Initialize temperature sequence
  temperature_seq = (/ (-5.0 + dble(i), i = 0, 45) /)  ! -5 to 40 in increments of 1

  ! Test calc_ftemp_inst_jmax with FV_kumarathunge19
  do i = 1, 46
    x = temperature_seq(i)
    ref = expected_kumarathunge(i)
    x = calc_ftemp_inst_jmax(real(x), 25.0, 25.0, 25.0, FV_kumarathunge19)
    print *, x, ref
    err_result = err(x, ref)
    if (err_result > err_tol) then
      print *, "Test failed for FV_kumarathunge19 at temperature ", temperature_seq(i)
      is_error = .true.
    end if
  end do

  ! Test calc_ftemp_inst_jmax with FV_kattge07
  do i = 1, 46
    x = temperature_seq(i)
    ref = expected_kattge(i)
    x = calc_ftemp_inst_jmax(real(x), 25.0, 25.0, 25.0, FV_kattge07)
    print *, x, ref
    err_result = err(x, ref)
    if (err_result > err_tol) then
      print *, "Test failed for FV_kattge07 at temperature ", temperature_seq(i)
      is_error = .true.
    end if
  end do

  ! Test calc_ftemp_inst_jmax with FV_leuning02
  do i = 1, 46
    x = temperature_seq(i)
    ref = expected_leuning(i)
    x = calc_ftemp_inst_jmax(real(x), 25.0, 25.0, 25.0, FV_leuning02)
    print *, x, ref
    err_result = err(x, ref)
    if (err_result > err_tol) then
      print *, "Test failed for FV_leuning02 at temperature ", temperature_seq(i)
      is_error = .true.
    end if
  end do

  if (is_error) then
    print *, "Some tests failed!"
  else
    print *, "All tests passed!"
  end if

  contains

  function err(x, ref) 
    real(kind = dbl8), intent(in) :: x, ref
    real(kind = dbl8) :: err
    err = min(abs(x - ref), abs(x / ref - 1))
  end function err

end subroutine test_jmax_temp_response


program main
  implicit none

  call test_par_env
  
  call test_dpsi_calc

  call test_jmax_temp_response

  call test_gs

  call test_dfdx

  call test_phydro

end program
