module md_precision
  implicit none

  integer, parameter :: int4=SELECTED_INT_KIND(4)
  integer, parameter :: flt4=SELECTED_REAL_KIND(6,37)
  integer, parameter :: dbl8=SELECTED_REAL_KIND(15,307)
end module md_precision