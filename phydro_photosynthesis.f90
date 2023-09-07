module md_phydro_photosynthesis
  use md_precision
  use md_photosynth, only: calc_kmm, calc_gammastar, calc_ftemp_kphio

  implicit none

  ! list of methods to model temperature dependencies of Vcmax and Jmax
  integer (kind = int4), parameter :: FV_kattge07 = 0, FV_kumarathunge19 = 1, FV_leuning02 = 2

  ! list of methods to model temperature dependencies of Rd
  integer (kind = int4), parameter :: FR_heskel16 = 0, FR_arrhenius = 1, FR_q10 = 2

  ! list of methods to model temperature dependencies of br
  integer (kind = int4), parameter :: FB_atkin15 = 0, FB_kumarathunge19 = 1

  type par_photosynth_type
    real (kind = dbl8) ::  kmm
    real (kind = dbl8) ::  gammastar
    real (kind = dbl8) ::  phi0
    real (kind = dbl8) ::  ca
    real (kind = dbl8) ::  delta

    integer(kind = int4) :: ftemp_vj_method;
    integer(kind = int4) :: ftemp_rd_method;
    integer(kind = int4) :: ftemp_br_method;
  
    real (kind = dbl8) ::  Iabs
    real (kind = dbl8) ::  patm

    real (kind = dbl8) ::  fT_vcmax;
    real (kind = dbl8) ::  fT_jmax;
    real (kind = dbl8) ::  fT_rd;
  
  end type par_photosynth_type

  type ACi_type
    real(kind=dbl8) :: a
    real(kind=dbl8) :: ci
    logical :: isVcmaxLimited
  end type ACi_type


  contains

  subroutine create_par_photosynth(this, tc, patm, kphio, co2, ppfd, fapar, rdark25, tcgrowth, tchome, &
                                   ftemp_vj_method, ftemp_rd_method, ftemp_br_method)
    
    type(par_photosynth_type), intent(out) :: this
    real(kind = dbl8), intent(in) :: tc, patm, kphio, co2, ppfd, fapar, rdark25, tcgrowth, tchome
    integer(kind = int4), intent(in) :: ftemp_vj_method, ftemp_rd_method, ftemp_br_method

    ! Calculate temperature scaling factors
    this%fT_vcmax = calc_ftemp_inst_vcmax(real(tc), real(tcgrowth), 25.0, ftemp_vj_method)
    this%fT_jmax = calc_ftemp_inst_jmax(real(tc), real(tcgrowth), real(tchome), 25.0, ftemp_vj_method)
    this%fT_rd = calc_ftemp_inst_rd(real(tc), ftemp_rd_method)

    ! Calculate other parameters
    this%kmm = calc_kmm(real(tc), real(patm))
    this%gammastar = calc_gammastar(real(tc), real(patm))
    this%phi0 = kphio * calc_ftemp_kphio(real(tc), .false.)
    this%Iabs = ppfd * fapar
    this%ca = co2 * patm * 1.0d-6
    this%patm = patm
    this%delta = rdark25 * this%fT_rd / this%fT_vcmax

    ! Set the temperature scaling methods
    this%ftemp_vj_method = ftemp_vj_method
    this%ftemp_rd_method = ftemp_rd_method
    this%ftemp_br_method = ftemp_br_method

  end subroutine create_par_photosynth

  subroutine print_par_photosynth(this)
    type(par_photosynth_type), intent(in) :: this

    print *, "ParPhotosynth: "
    print *, "   fT_vcmax", this%fT_vcmax
    print *, "   fT_jmax", this%fT_jmax
    print *, "   fT_rd", this%fT_rd
    print *, "   kmm", this%kmm
    print *, "   gammastar", this%gammastar
    print *, "   phi0", this%phi0
    print *, "   Iabs", this%Iabs
    print *, "   ca", this%ca
    print *, "   patm", this%patm
    print *, "   delta", this%delta
    print *, "   ftemp_vj_method", this%ftemp_vj_method
    print *, "   ftemp_rd_method", this%ftemp_rd_method
    print *, "   ftemp_br_method", this%ftemp_br_method
  end subroutine print_par_photosynth

  function calc_ftemp_arrhenius(tk, dha, tkref) result(ftemp)
    ! Output:   Factor fv to correct for instantaneous temperature response
    !           of Vcmax for:
    !
    !               Vcmax(temp) = fv * Vcmax(25 deg C) 
    !
    ! Input:
    !   tk      - Leaf temperature in Kelvin
    !   dha     - Activation energy (J/mol)
    !   tkref   - Reference temperature in Kelvin (default: 298.15 K)

    real(kind = dbl8), intent(in) :: tk, dha, tkref
    real(kind = dbl8), parameter :: kR = 8.3145 ! Universal gas constant, J/mol/K
    real(kind = dbl8) :: ftemp

    ! Calculate temperature scaling factor using Arrhenius equation
    ftemp = exp(dha * (tk - tkref) / (tkref * kR * tk))

  end function calc_ftemp_arrhenius


  function calc_ftemp_inst_vcmax(tcleaf, tcgrowth, tcref, method_ftemp) result(fv)
    real(kind = flt4), intent(in) :: tcleaf, tcgrowth, tcref
    real(kind = dbl8) :: fv 
    integer(kind = int4), intent(in) :: method_ftemp
    real(kind = dbl8), parameter :: Rgas = 8.3145 ! Universal gas constant (J/mol/K)
    real(kind = dbl8) :: tkref  
    real(kind = dbl8) :: tkleaf 
    real(kind = dbl8) :: Hd, Ha, a_ent, b_ent, dent, fva, fvb
    real(kind = dbl8) :: Sv, term_1, term_2, term_3

    tkref = tcref + 273.15 ! Convert reference temperature to Kelvin
    tkleaf = tcleaf + 273.15 ! Convert leaf temperature to Kelvin

    if (method_ftemp == FV_kattge07 .or. method_ftemp == FV_kumarathunge19) then
        ! Kattge2007 Parametrization
        Hd = 200000.0 ! Deactivation energy (J/mol)
        Ha = 71513.0 ! Activation energy (J/mol)
        a_ent = 668.39 ! Offset of entropy vs. temperature relationship from Kattge & Knorr (2007) (J/mol/K)
        b_ent = 1.07 ! Slope of entropy vs. temperature relationship from Kattge & Knorr (2007) (J/mol/K^2)

        if (method_ftemp == FV_kumarathunge19) then
            ! Kumarathunge2019 Implementation:
            ! local parameters
            a_ent = 645.13 ! Offset of entropy vs. temperature relationship (J/mol/K)
            b_ent = 0.38 ! Slope of entropy vs. temperature relationship (J/mol/K^2)
            
            ! local variables
            Ha = 42600.0 + (1140.0 * tcgrowth) ! Acclimation for vcmax
        end if

        ! Calculate entropy following Kattge & Knorr (2007), negative slope and y-axis intersect is when expressed as a function of temperature in degrees Celsius, not Kelvin!
        dent = a_ent - (b_ent * tcgrowth)  ! 'tcgrowth' corresponds to 'tmean' in Nicks, 'tc25' is 'to' in Nick's

        fva = calc_ftemp_arrhenius(tkleaf, Ha, tkref)
        fvb = (1.0 + exp((tkref * dent - Hd) / (Rgas * tkref))) / (1.0 + exp((tkleaf * dent - Hd) / (Rgas * tkleaf)))
        fv = fva * fvb
    elseif (method_ftemp == FV_leuning02) then
        ! Ref: Leuning, R. (2002). Temperature dependence of two parameters in a photosynthesis model. Plant, Cell & Environment, 25(9), 1205–1210. https://doi.org/10.1046/j.1365-3040.2002.00898.x
        ! Table 2:
        Ha = 73637.0
        Hd = 149252.0
        Sv = 486.0

        term_1 = 1.0 + exp((Sv * tkref - Hd) / (Rgas * tkref))
        term_3 = 1.0 + exp((Sv * tkleaf - Hd) / (Rgas * tkleaf))
        term_2 = exp((Ha / (Rgas * tkref)) * (1.0 - tkref / tkleaf)) ! Careful: In Eq. (1) in Leuning et al. (1992), there is a bracket missing in this term!

        fv = term_1 * term_2 / term_3
    else
        write(*,*) "Invalid method_ftemp:", method_ftemp
        stop
    end if
  end function calc_ftemp_inst_vcmax


  function calc_ftemp_inst_jmax(tcleaf, tcgrowth, tchome, tcref, method_ftemp) result(fv)
    real(kind = flt4), intent(in) :: tcleaf, tcgrowth, tchome, tcref
    integer(kind = int4), intent(in) :: method_ftemp

    real(kind = dbl8), parameter :: Rgas = 8.3145 ! Universal gas constant (J/mol/K)
    real(kind = dbl8) :: tkref ! Convert reference temperature to Kelvin
    real(kind = dbl8) :: tkleaf ! Convert leaf temperature to Kelvin
    real(kind = dbl8) :: fv

    real(kind = dbl8) :: Hd ! Deactivation energy (J/mol)
    real(kind = dbl8) :: Ha ! Activation energy (J/mol)
    real(kind = dbl8) :: a_ent ! Offset of entropy vs. temperature relationship from Kattge & Knorr (2007) (J/mol/K)
    real(kind = dbl8) :: b_ent ! Slope of entropy vs. temperature relationship from Kattge & Knorr (2007) (J/mol/K^2)
    real(kind = dbl8) :: c_ent
    real(kind = dbl8) :: dent ! Entropy calculation, equations given in Celsius, not in Kelvin
    real(kind = dbl8) :: fva
    real(kind = dbl8) :: fvb

    real(kind = dbl8) :: Sv, term_1, term_2, term_3

    tkref = tcref + 273.15
    tkleaf = tcleaf + 273.15

    if (method_ftemp == FV_kattge07 .or. method_ftemp == FV_kumarathunge19) then
      Hd = 200000.0
      Ha = 49884.0
      a_ent = 659.70
      b_ent = 0.75

      dent = a_ent - b_ent * tcgrowth

      if (method_ftemp == FV_kumarathunge19) then
        Ha = 40710.0
        a_ent = 658.77
        b_ent = 0.84
        c_ent = 0.52

        dent = a_ent - (b_ent * tchome) - c_ent * (tcgrowth - tchome)
      end if

      fva = calc_ftemp_arrhenius(tkleaf, Ha, tkref)
      fvb = (1.0 + exp((tkref * dent - Hd) / (Rgas * tkref))) / (1.0 + exp((tkleaf * dent - Hd) / (Rgas * tkleaf)))
      fv = fva * fvb

    elseif (method_ftemp == FV_leuning02) then
      Ha = 50300.0
      Hd = 152044.0
      Sv = 495.0

      term_1 = 1.0 + exp((Sv * tkref - Hd) / (Rgas * tkref))
      term_3 = 1.0 + exp((Sv * tkleaf - Hd) / (Rgas * tkleaf))
      term_2 = exp((Ha / (Rgas * tkref)) * (1.0 - tkref / tkleaf))

      fv = term_1 * term_2 / term_3

    else
        write(*,*) "Invalid method_ftemp:", method_ftemp
        stop
    end if

  end function calc_ftemp_inst_jmax


  function calc_ftemp_inst_rd(tc_leaf, method_rd_scale) result(f)
    real(kind = dbl8) :: f
    real(kind = flt4), intent(in) :: tc_leaf
    integer(kind=int4), intent(in) :: method_rd_scale
    real(kind = dbl8) :: apar, bpar, dha

    if (method_rd_scale == FR_heskel16) then
      ! Heskel et al. (2016) temperature scaling
      apar = 0.1012
      bpar = 0.0005
      f = exp(apar * (tc_leaf - 25.0) - bpar * (tc_leaf*tc_leaf - 25.0*25.0))
    elseif (method_rd_scale == FR_arrhenius) then
      ! Arrhenius temperature scaling
      dha = 20700.0 ! Activation energy taken from Kumarathunge et al. (2019), Table 1, Mature Natural Environment
      f = calc_ftemp_arrhenius(dble(tc_leaf) + 273.15, dha, 298.15_dbl8) ! Convert temperature to Kelvin and call calc_ftemp_arrh function
    elseif (method_rd_scale == FR_q10) then
      ! Q10 temperature scaling according to Tjoelker et al. (2001)
      f = (3.22 - 0.046 * tc_leaf)**(tc_leaf - 25.0) / 10.0
    else
      write(*,*) "Invalid method_rd_scale:", method_rd_scale
      stop
    end if

  end function calc_ftemp_inst_rd


  function calc_brd25(method_rd25, tc_growth) result(rd_to_vcmax)
    real(kind = dbl8) :: rd_to_vcmax
    real(kind = dbl8), intent(in) :: tc_growth
    integer(kind = int4), intent(in) :: method_rd25

    if (method_rd25 == FB_atkin15) then
      rd_to_vcmax = 0.015 ! Ratio of Rdark to Vcmax25, Atkin et al., 2015 for C3 herbaceous
    elseif (method_rd25 == FB_kumarathunge19) then
      rd_to_vcmax = 0.0360 - 0.0010 * tc_growth ! Acclimated rd_to_vcmax taken from Kumarathunge et al. (2019), Table 1, Mature Natural Environment
    else
      write(*,*) "Invalid method_rd25:", method_rd25
      stop
    end if

  end function calc_brd25


  !-------------------------------------------------------
  !   Ac / Aj calculations
  !-------------------------------------------------------
  function QUADM(A, B, C)
    real(kind=dbl8) :: QUADM
    real(kind=dbl8), intent(in) :: A, B, C
    QUADM = (-B - sqrt(B*B - 4.0d0*A*C)) / (2.0d0*A)
  end function QUADM

  function QUADP(A, B, C)
    real(kind=dbl8) :: QUADP
    real(kind=dbl8), intent(in) :: A, B, C
    QUADP = (-B + sqrt(B*B - 4.0d0*A*C)) / (2.0d0*A)
  end function QUADP


  function calc_assim_rubisco_limited(gs_in, vcmax, par_photosynth) result(res)
    real(kind=dbl8), intent(in) :: gs_in
    real(kind=dbl8), intent(in) :: vcmax
    type(ACi_type) :: res
    type(par_photosynth_type) :: par_photosynth
    real(kind=dbl8) :: ca, d, A, B, C, gs

    gs = gs_in

    ca = par_photosynth%ca
    gs = gs * 1.0d6 / par_photosynth%patm
    d = par_photosynth%delta

    A = -gs
    B = gs * ca - gs * par_photosynth%kmm - vcmax*(1.0d0-d)
    C = gs * ca * par_photosynth%kmm + vcmax * (par_photosynth%gammastar + par_photosynth%kmm*d)

    res%ci = QUADM(A, B, C)
    res%a = gs * (ca - res%ci)
    res%isVcmaxLimited = .true.

  end function calc_assim_rubisco_limited


  function calc_assim_light_limited(gs_in, jmax, par_photosynth) result(res)
    real(kind=dbl8), intent(in) :: gs_in
    real(kind=dbl8), intent(in) :: jmax
    type(ACi_type) :: res
    type(par_photosynth_type) :: par_photosynth
    real(kind=dbl8) :: ca, d, phi0iabs, jj, jlim, A, B, C, gs

    gs = gs_in

    ca = par_photosynth%ca
    gs = gs * 1.0d6 / par_photosynth%patm
    gs = gs + 1.0d-12
    d = par_photosynth%delta

    phi0iabs = par_photosynth%phi0 * par_photosynth%Iabs
    jj = 4.0d0 * phi0iabs / jmax
    jlim = phi0iabs / sqrt(1.0d0 + jj*jj)

    A = -gs
    B = gs * ca - gs * 2.0d0 * par_photosynth%gammastar - jlim * (1.0d0-d)
    C = gs * ca * 2.0d0 * par_photosynth%gammastar + jlim * (par_photosynth%gammastar + d*par_photosynth%kmm)

    res%ci = QUADM(A, B, C)
    res%a = gs * (ca - res%ci)
    res%isVcmaxLimited = .false.

  end function calc_assim_light_limited


  function calc_assimilation_limiting(vcmax, jmax, gs, par_photosynth) result(Aout)
    real(kind=dbl8), intent(in) :: vcmax, jmax
    real(kind=dbl8), intent(in) :: gs
    type(ACi_type) :: Ac, Aj, Aout
    type(par_photosynth_type) :: par_photosynth

    Ac = calc_assim_rubisco_limited(gs, vcmax, par_photosynth)
    Aj = calc_assim_light_limited(gs, jmax, par_photosynth)

    if (Ac%ci > Aj%ci) then
      Aout = Ac
    else
      Aout = Aj
    end if
  end function calc_assimilation_limiting  

end module md_phydro_photosynthesis

