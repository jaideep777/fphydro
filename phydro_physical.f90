module md_phydro_physical
  use md_precision
  implicit none
  
  contains
  
  function calc_esat(TdegC, patm) result(esatval)
    real(kind = dbl8), intent(in) :: TdegC, patm
    real(kind = dbl8) :: esatval
    real(kind = dbl8) :: a, b, c, f
  
    a = 611.21
    b = 17.502
    c = 240.97
    f = 1.0007 + 3.46e-8 * patm
  
    esatval = f * a * exp(b * TdegC / (c + TdegC))
  end function calc_esat
  
  function calc_density_air(tc_air, patm, vpd, moist) result(rho)
    real(kind = dbl8), intent(in) :: tc_air, patm, vpd
    logical, intent(in) :: moist
    real(kind = dbl8) :: rho, tk, R
    real(kind = dbl8) :: vp, rv, tv
  
    tk = tc_air + 273.16
    R = 287.052874
  
    if (.not. moist) then
    rho = patm / R / tk
    else
    vp = calc_esat(tc_air, patm) - vpd
    rv = 0.622 * vp / (patm - vp)
    tv = tk * (1.0 + rv / 0.622) / (1.0 + rv)
  
    rho = patm / R / tv
    end if
  end function calc_density_air
  
  function calc_enthalpy_vap(tc) result(enthalpy)
    real(kind = dbl8), intent(in) :: tc
    real(kind = dbl8) :: enthalpy, tk, a
  
    tk = tc + 273.15
    a = tk / (tk - 33.91)
  
    enthalpy = 1.91846e6 * a**2
  end function calc_enthalpy_vap
  
  function calc_cp_moist_air(tc) result(cp)
    real(kind = dbl8), intent(in) :: tc
    real(kind = dbl8) :: cp, my_tc
  
    my_tc = max(min(tc, 100.0), 0.0)
    
    cp =         (1.0045714270 +     &
       my_tc * (2.050632750e-3 +   &
       my_tc * (-1.631537093e-4 +  &
       my_tc * (6.212300300e-6 -   &
       my_tc * (8.830478888e-8 -   &
       my_tc * 5.071307038e-10))))) * 1e3
       
  end function calc_cp_moist_air
  
  function calc_psychro(tc, patm) result(psychro)
    real(kind = dbl8), intent(in) :: tc, patm
    real(kind = dbl8) :: psychro, Ma, Mv, cp, lv
  
    Ma = 0.02896
    Mv = 0.018016
  
    cp = calc_cp_moist_air(tc)
    lv = calc_enthalpy_vap(tc)
  
    psychro = cp * patm / ((Mv / Ma) * lv)
  end function calc_psychro
  
  function calc_sat_slope(tc) result(slope)
    real(kind = dbl8), intent(in) :: tc
    real(kind = dbl8) :: slope
  
    slope = 17.269 * 237.3 * 610.78 * exp(tc * 17.269 / (tc + 237.3)) / ((tc + 237.3)**2)
  end function calc_sat_slope
  
  
end module md_phydro_physical