! This is include file "diags_vrt3d.h"
!  ==== == ======= ==== ==========
!

#ifdef DIAGNOSTICS_VRT3D

      real vrt3dXadv(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE vrt3dXadv(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrt3dstretch(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE vrt3dstretch(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrt3dadvecti(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE vrt3dadvecti(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrt3dYadv(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE vrt3dYadv(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrt3dVadv(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE vrt3dVadv(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrt3dPrsgrd(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE vrt3dPrsgrd(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrt3dCor(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE vrt3dCor(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrt3dVmix(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE vrt3dVmix(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrt3dHmix(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE vrt3dHmix(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrt3dHdiff(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE vrt3dHdiff(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrt3dVmix2(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE vrt3dVmix2(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrt3drate(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE vrt3drate(BLOCK_PATTERN,*) BLOCK_CLAUSE
# if defined DIAGNOSTICS_BARO
      real vrtBaro(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE vrtBaro(BLOCK_PATTERN,*) BLOCK_CLAUSE
# endif
# if defined M3FAST
      real vrt3dfast(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE vrt3dfast(BLOCK_PATTERN,*) BLOCK_CLAUSE
# endif
# ifdef AVERAGES
      real timediags_vrt_avg
      real vrt3dXadv_avg(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE vrt3dXadv_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrt3dstretch_avg(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE vrt3dstretch_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrt3dadvecti_avg(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE vrt3dadvecti_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrt3dYadv_avg(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE vrt3dYadv_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrt3dVadv_avg(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE vrt3dVadv_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrt3dPrsgrd_avg(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE vrt3dPrsgrd_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrt3dCor_avg(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE vrt3dCor_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrt3dVmix_avg(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE vrt3dVmix_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrt3dHmix_avg(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE vrt3dHmix_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrt3dHdiff_avg(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE vrt3dHdiff_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrt3drate_avg(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE vrt3drate_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vrt3dVmix2_avg(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE vrt3dVmix2_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
#  if defined DIAGNOSTICS_BARO
      real vrtBaro_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE vrtBaro_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
#  endif
#  if defined M3FAST
      real vrt3dfast_avg(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE vrt3dfast_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
#  endif
# endif



      common /diag_vrt3dXadv/vrt3dXadv
     &       /diag_vrt3dstretch/vrt3dstretch
     &       /diag_vrt3dadvecti/vrt3dadvecti
     &       /diag_vrt3dYadv/vrt3dYadv
     &       /diag_vrt3dVadv/vrt3dVadv
     &       /diag_vrt3dPrsgrd/vrt3dPrsgrd
     &       /diag_vrt3dCor/vrt3dCor
     &       /diag_vrt3dVmix/vrt3dVmix
     &       /diag_vrt3dHmix/vrt3dHmix
     &       /diag_vrt3dHdiff/vrt3dHdiff
     &       /diag_vrt3drate/vrt3drate
     &       /diag_vrt3dVmix2/vrt3dVmix2
# if defined DIAGNOSTICS_BARO
     &       /diag_vrtBaro/vrtBaro
# endif
# if defined M3FAST
     &       /diag_vrt3dfast/vrt3dfast
# endif

# ifdef AVERAGES
      common /diag_timediags_vrt_avg/timediags_vrt_avg
      common /diag_vrt3dXadv_avg/vrt3dXadv_avg
     &       /diag_vrt3dstretch_avg/vrt3dstretch_avg
     &       /diag_vrt3dadvecti_avg/vrt3dadvecti_avg
     &       /diag_vrt3dYadv_avg/vrt3dYadv_avg
     &       /diag_vrt3dVadv_avg/vrt3dVadv_avg
     &       /diag_vrt3dPrsgrd_avg/vrt3dPrsgrd_avg
     &       /diag_vrt3dCor_avg/vrt3dCor_avg
     &       /diag_vrt3dVmix_avg/vrt3dVmix_avg
     &       /diag_vrt3dHmix_avg/vrt3dHmix_avg
     &       /diag_vrt3dHdiff_avg/vrt3dHdiff_avg
     &       /diag_vrt3drate_avg/vrt3drate_avg
     &       /diag_vrt3dVmix2_avg/vrt3dVmix2_avg
#  if defined DIAGNOSTICS_BARO
     &       /diag_vrtBaro_avg/vrtBaro_avg
#  endif
#  if defined M3FAST
     &       /diag_vrt3dfast_avg/vrt3dfast_avg
#  endif
# endif

#endif /* DIAGNOSTICS_VRT3D */




