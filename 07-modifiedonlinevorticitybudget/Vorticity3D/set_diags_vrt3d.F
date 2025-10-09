#include "cppdefs.h"
#if defined DIAGNOSTICS_VRT3D


      subroutine set_diags_vrt3d (tile)
      implicit none
      integer tile
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      call set_diags_vrt3d_tile (istr,iend,jstr,jend,tile)
      return
      end

      subroutine set_diags_vrt3d_tile (istr,iend,jstr,jend,tile)
!
      implicit none
# include "param.h"
      integer istr,iend,jstr,jend, i,j, ilc, iflux,
     &     imin,imax,jmin,jmax,tile
# ifdef SOLVE3D
     &                       , itrc, k
# endif
      real cff,cff1,cff2,
     &     dH(N), jstri(2), istri(2),
     &     jendi(2), iendi(2)
# include "scalars.h"
# include "ncscrum.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# if defined DIAGNOSTICS_UV || defined DIAGNOSTICS_VRT3D
#  include "diagnostics.h"
# endif
# include "diags_vrt3d.h"
# include "compute_auxiliary_bounds.h"

! ------- ------------- ------ ------ - -----

        jstri(1) = jstr-1
        jstri(2) = jstr

        istri(1) = istr
        istri(2) = istr-1

        jendi(1) = jend
        jendi(2) = jend

        iendi(1) = iend
        iendi(2) = iend

        !jstri(1) = jstrR
        !jstri(2) = jstrR

        !istri(1) = istrR
        !istri(2) = istrR

        !jendi(1) = jendR
        !jendi(2) = jendR

        !iendi(1) = iendR
        !iendi(2) = iendR

      if (WESTERN_EDGE) then     ! Restrict extended ranges one
        imin=istr+1              ! point inward near the physical
      else                       ! boundary. Note that this version
        imin=istr                ! of code is suitable for MPI
      endif                      ! configuration too.

      if (SOUTHERN_EDGE) then
        jmin=jstr+1
      else
        jmin=jstr
      endif

! ------- --------- -- -- -------- --------
! Compute Vorticity of 3D momentum equation
! ------- --------- -- -- -------- --------
        do j=jmin,jend
          do i=imin,iend
            cff = 0.25*(pm(i,j) + pm(i-1,j) + pm(i,j-1) + pm(i-1,j-1))
     &         * 0.25*(pn(i,j) + pn(i-1,j) + pn(i,j-1) + pn(i-1,j-1))
            cff1 = cff
            cff2 = cff
            cff = 0.25*(pm(i,j) + pm(i-1,j) + pm(i,j-1) + pm(i-1,j-1))

            do k=1,N
! Rate
              vrt3drate(i,j,k) = (Mrate(i,j,k,2)-Mrate(i-1,j,k,2)) * cff1/cff
     &           - (Mrate(i,j,k,1)-Mrate(i,j-1,k,1))*cff2/cff
! x-advection          
              vrt3dXadv(i,j,k) = (MXadv(i,j,k,2)-MXadv(i-1,j,k,2)) * cff1/ cff
     &           - (MXadv(i,j,k,1) - MXadv(i,j-1,k,1)) * cff2/cff
! y-advection          
              vrt3dYadv(i,j,k) = (MYadv(i,j,k,2)-MYadv(i-1,j,k,2)) * cff1/cff
     &           - (MYadv(i,j,k,1) - MYadv(i,j-1,k,1)) * cff2/cff
! vertical advection
              vrt3dVadv(i,j,k) = (MVadv(i,j,k,2)-MVadv(i-1,j,k,2)) * cff1/cff
     &           - (MVadv(i,j,k,1) - MVadv(i,j-1,k,1)) * cff2/cff
! Coriolis 
              vrt3dCor(i,j,k) = (MCor(i,j,k,2)-MCor(i-1,j,k,2))*cff1/cff
     &           - (MCor(i,j,k,1) - MCor(i,j-1,k,1)) * cff2/cff
! Vmix
              vrt3dVmix(i,j,k) = (MVmix(i,j,k,2)-MVmix(i-1,j,k,2))*cff1/cff
     &           - (MVmix(i,j,k,1) - MVmix(i,j-1,k,1)) * cff2/cff
! Pressure gradient
              vrt3dPrsgrd(i,j,k) = (MPrsgrd(i,j,k,2)-MPrsgrd(i-1,j,k,2))*cff1/cff
     &           - (MPrsgrd(i,j,k,1) - MPrsgrd(i,j-1,k,1)) * cff2/cff
! Horizontal mixing
              vrt3dHmix(i,j,k) = (MHmix(i,j,k,2,nstp)-MHmix(i-1,j,k,2,nstp))*cff1/cff
     &           - (MHmix(i,j,k,1,nstp) - MHmix(i,j-1,k,1,nstp)) * cff2/cff
! Horizontal diffusion
              vrt3dHdiff(i,j,k) = (MHdiff(i,j,k,2)-MHdiff(i-1,j,k,2))*cff1/cff
     &           - (MHdiff(i,j,k,1) - MHdiff(i,j-1,k,1)) * cff2/cff
! Nudging term
              vrt3dVmix2(i,j,k) = (MVmix2(i,j,k,2)-MVmix2(i-1,j,k,2))*cff1/cff
     &           - (MVmix2(i,j,k,1) - MVmix2(i,j-1,k,1)) * cff2/cff

#  if defined M3FAST   
! Fast
              vrt3dfast(i,j,k)=(Mfast(i,j,k,2)-Mfast(i-1,j,k,2))*cff1/cff
     &           - (Mfast(i,j,k,1) - Mfast(i,j-1,k,1)) * cff2/cff
#  endif                
            enddo
          
# if defined DIAGNOSTICS_BARO
            vrtBaro(i,j)=(wrkBaro(i,j,2)-wrkBaro(i-1,j,2))*cff1
     &           - (wrkBaro(i,j,1) - wrkBaro(i,j-1,1)) * cff2
# endif
          enddo
        enddo

! ------- ------------- ------ ------ - -----
      return
      end
#else /* DIAGNOSTICS_VRT3D*/
      subroutine set_diags_vrt3d_empty
      end
#endif /* DIAGNOSTICS_VRT3D */
