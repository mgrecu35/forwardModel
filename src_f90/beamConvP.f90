! ###########################################################
!            Satellite Data Simulation Unit v2
!                     - BMCVparameters -
!
! This module defines variables for use by the beam convolution
! routine (BEAM_CONV.f90).
!
!        Updated from BMCVparameters.inc by H. M. (May 2009)
!
 Module BMCVparameters
!
!  ######## microwave radiometer #######
!
! fov_ct_microw : cross-track FOV size at each channel
! fov_dt_microw : down-track FOV size at each channel
!                 (both defined in 'param_set_bmcv.f90')
!
   real, dimension(:), allocatable :: fov_ct_microw
   real, dimension(:), allocatable :: fov_dt_microw
!
!  ############## radar ###############
!
! fov_ct_radar : cross-track FOV size at each channel
! fov_dt_radar : down-track FOV size at each channel
!                 (both defined in 'param_set_bmcv.f90')
!
   real, dimension(:), allocatable :: fov_ct_radar
   real, dimension(:), allocatable :: fov_dt_radar
!
!  ######## visible/IR imager ##########
!
! fov_ct_visir : cross-track FOV size
! fov_dt_visir : down-track FOV size
!
   real, dimension(:), allocatable :: fov_ct_visir
   real, dimension(:), allocatable :: fov_dt_visir
!
! ##########      Work       ##########
!
   Real, dimension(:,:), allocatable :: radinp
!
 End Module BMCVparameters
