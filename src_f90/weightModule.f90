module weight
  real, dimension(:,:), pointer :: weightTable10
  real, dimension(:,:), pointer :: weightTable19
  real, dimension(:,:), pointer :: weightTable21
  real, dimension(:,:), pointer :: weightTable37
  real, dimension(:,:), pointer :: weightTable85
  integer, dimension(:), pointer  :: iflag
contains
  subroutine initWFlag(nfreq)
    integer :: nfreq
    if(associated(iflag)) return
    allocate(iflag(nfreq))
    iflag=0
  end subroutine initWFlag
end module weight
