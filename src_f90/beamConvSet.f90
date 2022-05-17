!  SFM  04/06/2013  Code changes from M.Grecu
!
! ############################################################
!            Satellite Data Simulation Unit v2
!                     -  param_set_bmcv -
!
! User-defined parameters for the beam convolution routine
! are set up in this subroutine.
!
!              Updated from param_set-BMCV.f by H. M. (May, 2009)
!
 Subroutine param_set_BMCV(mxfreq_microw)
   Use BMCVparameters
   Implicit NONE
!
! 1. Microwave radiometer -------------------------------------------
!
! 1.2 Sensor field of views in [km]*
!    (*: the unit of FOV can be arbitrarily defined but must be shared
!        by "delta_x" and "delta_y".)
!
   integer :: mxfreq_microw
   if(allocated(fov_ct_microw)) return
   Allocate ( fov_ct_microw(mxfreq_microw), fov_dt_microw(mxfreq_microw) )
!
!  Cross-track (x-direction) FOV sizes
!
   fov_dt_microw(1) = 36./2.
   fov_ct_microw(1) = 60./2.
   if(mxfreq_microw==1) return
   fov_dt_microw(2) = 18./2.
   fov_ct_microw(2) = 30./2.
   if(mxfreq_microw==2) return
   fov_dt_microw(3) = 17./2.
   fov_ct_microw(3) = 27./2.
   if(mxfreq_microw==3) return
   fov_dt_microw(4) = 9.7/2.
   fov_ct_microw(4) = 16./2.
   if(mxfreq_microw==4) return
   fov_dt_microw(5) = 4.2
   fov_ct_microw(5) = 6.8
!
!  Down-track (y-direction) FOV sizes
!

   
 End Subroutine param_set_BMCV
