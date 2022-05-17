module missingMod
  integer *1 :: missing_i1
  integer *2 :: missing_i2
  integer *4 :: missing_i4
  real *4 :: missing_r4
  real *8 :: missing_r8
contains
  subroutine init_missing_flags
    missing_i1 = -99
    missing_i2 = -9999
    missing_i4 = -9999
    missing_r4 = -9999.9
    missing_r8 = -9999.9
  end subroutine init_missing_flags
end module missingMod

