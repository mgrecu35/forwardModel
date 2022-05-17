module nbinMod
  integer :: nbin
  integer :: imemb
  real    :: dpia13srt(100)
  integer :: n9(9)
  integer :: iTcount
!begin  WSO 8/8/13
  integer :: ntransition
!end    WSO 8/8/13
contains
  subroutine init_nbin
    use ran_mod
    nbin=88
    iTcount=0
    do i=1,100
       dpia13srt(i)=normal2(0.,1.)
    enddo
!begin  WSO 8/14/13
    ntransition = 10
!end    WSO 8/14/13
    
  end subroutine init_nbin
end module nbinMod
