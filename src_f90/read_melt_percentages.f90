module Tables_frac

   ! define table variables for storing melt liquid water content
   ! and melt rate fractions based upon Liang Liao and Bob Meneghinis
   ! melting model calculations

   ! Bill Olson  July, 2013

      integer :: d0_nx, height_ny, mu_nz  ! number of d0, height, and mu values represented in tables

      real :: mu_init, mu_interval ! first mu [] and mu interval [] values
      real :: d0_init, d0_interval ! first d0 [mm] and d0 interval [mm] values
      real :: height_init, height_interval ! first height [m] and height interval [m] values

      integer, allocatable :: max_j(:, :) !level of maximum reflectivity
      real, allocatable :: mlwc_fraction(:, :, :), mrate_fraction(:, :, :) ! mass fraction of melt lwc and rate

end module Tables_frac

      subroutine read_melt_percentages

   ! read tables of melt liquid water content
   ! and melt rate fractions based upon Liang Liao and Bob Meneghinis
   ! melting model calculations

   ! Bill Olson  July, 2013

      use Tables_frac

      implicit none
 
      logical :: printflag
      integer :: i, j, k

      printflag = .true.

!begin  WSO 12/12/13 updated melting fraction file to Dm indexing
      open(11, file = 'MeltFrac/melt_lwc_rate_percentages_Dm.rho100.dat', &
       form = 'formatted', status='old', action='read')
!end    WSO 12/12/13
      read(11, '(i10, f7.2, f10.4)') mu_nz, mu_init, mu_interval
      read(11, '(i10, f7.4, f12.8)') d0_nx, d0_init, d0_interval
      read(11, '(i10, f7.1, f10.6)') height_ny, height_init, height_interval

      if(printflag) write(6, '("mu_nz:  ", i5, "  mu_init: ", f7.2, "  mu_interval: ", f10.4)'), mu_nz, mu_init, mu_interval
      if(printflag) write(6, '("d0_nx:  ", i5, "  d0_init: ", f7.4, "  d0_interval: ", f12.8)'), d0_nx, d0_init, d0_interval
      if(printflag) write(6, '("height_ny:  ", i5, "  height_init: ", f7.1, "  height_interval: ", f10.6)'), &
       height_ny, height_init, height_interval
 
      allocate(max_j(d0_nx, mu_nz), mlwc_fraction(d0_nx, height_ny, mu_nz), mrate_fraction(d0_nx, height_ny, mu_nz))


      do k = 1, mu_nz

    ! read lwc melt fractions
        do i = 1, d0_nx
          read(11, '(i5, 2x, 51f5.2)') max_j(i, k), (mlwc_fraction(i, j, k), j = 1, height_ny)
        end do

    ! read rate melt fractions
        do i = 1, d0_nx
          read(11, '(i5, 2x, 51f5.2)') max_j(i, k), (mrate_fraction(i, j, k), j = 1, height_ny)
        end do

      end do

      if(printflag) write(6, '("read the percentages in")') 

      return
      end subroutine read_melt_percentages
