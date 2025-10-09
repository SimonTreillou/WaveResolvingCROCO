!======================================================================
! CROCO is a branch of ROMS developped at IRD, INRIA, 
! Ifremer, CNRS and Univ. Toulouse III  in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
#include "cppdefs.h"
#if defined DIAGNOSTICS_VRT3D
!
!---------------------------------------------------------------
! Write diagnostics fields at requested levels into diagnostics
! netCDF file.

!---------------------------------------------------------------
!
      subroutine wrt_diags_vrt3d

      implicit none
      integer ierr,  record, lstr, lvar, lenstr
     &   , start(2), count(2), ibuff(4), nf_fwrite, type
# if defined MPI & !defined PARALLEL_FILES
      include 'mpif.h'
      integer status(MPI_STATUS_SIZE), blank
# endif
# include "param.h"
# include "scalars.h"
# include "ncscrum.h"
# include "forces.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "mixing.h"
# include "diags_vrt3d.h"
# include "mpi_cpl.h"
# ifdef SEDIMENT
#  include "sediment.h"
# endif
# ifdef BBL
#  include "bbl.h"
# endif
# ifdef SOLVE3D
      integer tile,itrc,i,j,k,ivar
#  ifdef SEDIMENT
     &      , indxWrk
#  endif
#  include "work.h"
# endif
# include "netcdf.inc"
# if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
      if (mynode.gt.0) then
        call MPI_Recv (blank, 1, MPI_INTEGER, mynode-1,
     &                 1, MPI_COMM_WORLD, status, ierr)
      endif
# endif
!
! Create/open diagnostic file; write grid arrays, if so needed.
!
      call def_diags_vrt3d (nciddiags_vrt3d, nrecdiags_vrt3d, ierr)
      if (ierr .ne. nf_noerr) goto 99
      lstr=lenstr(diags_vrt3dname)
!                                          !!! WARNING: Once time
! Set record within the file.              !!! stepping has been
!
                                           !!! started, it is assumed
      if (iic.eq.0) nrecdiags_vrt3d=nrecdiags_vrt3d+1  !!! that global history
      if (nrpfdiags_vrt3d.eq.0) then            !!! record index "nrecdiags_vrt3d"
        record=nrecdiags_vrt3d                  !!! is advanced by main.
      else
        record=1+mod(nrecdiags_vrt3d-1, nrpfdiags_vrt3d)
      endif
!
!---------------------------------------------------------------
! Write out evolving model variables:
!---------------------------------------------------------------
!
! Time step number and record numbers.
!
      type=filetype_diags_vrt
!
      ibuff(1)=iic
      ibuff(2)=nrecrst
      ibuff(3)=nrechis
      ibuff(4)=nrecdiags_vrt3d
      start(1)=1
      start(2)=record
      count(1)=4
      count(2)=1
      ierr=nf_put_vara_int (nciddiags_vrt3d, diags_vrt3dTstep, start, count, ibuff)
      if (ierr .ne. nf_noerr) then
        write(stdout,1) 'time_step', record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! Time
!
      ierr=nf_put_var1_FTYPE (nciddiags_vrt3d, diags_vrt3dTime, record, time)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxTime))
        write(stdout,1) vname(1,indxTime)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! Time2
!
      ierr=nf_put_var1_FTYPE (nciddiags_vrt3d, diags_vrt3dTime2, record, time)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxTime2))
        write(stdout,1) vname(1,indxTime2)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif

!---------------------------------------------------------------
! Momentum diagnostic variables.
!---------------------------------------------------------------
!
!
        itrc = 1
          if (wrtdiags_vrt3d(itrc)) then
            ivar=p2dvar
!
! indxvrt3drate
!
      !     work2d=vrtrateS
           call fillvalue3d(vrt3drate,nciddiags_vrt3d,diags_vrt3drate(itrc),
     &        vname(1,indxvrt3drate+itrc-1),record,r3dvar,type)
!
! indxvrt3dXadv
!
           call fillvalue3d(vrt3dXadv,nciddiags_vrt3d,diags_vrt3dXadv(itrc),
     &                      vname(1,indxvrt3dXadv+itrc-1),
     &                      record,r3dvar,type)
!
! indxvrt3dYadv
!
           call fillvalue3d(vrt3dYadv,nciddiags_vrt3d,diags_vrt3dYadv(itrc),
     &                      vname(1,indxvrt3dYadv+itrc-1),
     &                      record,r3dvar,type)
!
! indxvrt3dVadv
!
           call fillvalue3d(vrt3dVadv,nciddiags_vrt3d,diags_vrt3dVadv(itrc),
     &                      vname(1,indxvrt3dVadv+itrc-1),
     &                      record,r3dvar,type)
!
! indxvrt3dPrsgrd
!
           call fillvalue3d(vrt3dPrsgrd,nciddiags_vrt3d,diags_vrt3dPrsgrd(itrc),
     &                      vname(1,indxvrt3dPrsgrd+itrc-1),
     &                      record,r3dvar,type) 
!
! indxvrt3dCor
!
           call fillvalue3d(vrt3dCor,nciddiags_vrt3d,diags_vrt3dCor(itrc),
     &                      vname(1,indxvrt3dCor+itrc-1),
     &                      record,r3dvar,type)  
!
! indxvrt3dVmix
!
           call fillvalue3d(vrt3dVmix,nciddiags_vrt3d,diags_vrt3dVmix(itrc),
     &                      vname(1,indxvrt3dVmix+itrc-1),
     &                      record,r3dvar,type)
!
! indxvrt3dHmix
!
           call fillvalue3d(vrt3dHmix,nciddiags_vrt3d,diags_vrt3dHmix(itrc),
     &                      vname(1,indxvrt3dHmix+itrc-1),
     &                      record,r3dvar,type)
!
! indxvrt3dHdiff
!
           call fillvalue3d(vrt3dHdiff,nciddiags_vrt3d,diags_vrt3dHdiff(itrc),
     &                      vname(1,indxvrt3dHdiff+itrc-1),
     &                      record,r3dvar,type) 

!
! indxvrt3dstretch
!
           call fillvalue3d(vrt3dstretch,nciddiags_vrt3d,diags_vrt3dstretch(itrc),
     &                      vname(1,indxvrt3dstretch+itrc-1),
     &                      record,r3dvar,type)
     
!
! indxvrt3dadvecti
!
!           call fillvalue3d(vrt3dadvecti,nciddiags_vrt3d,diags_vrt3dadvecti(itrc),
!     &                      vname(1,indxvrt3dadvecti+itrc-1),
!     &                      record,r3dvar,type)  

!
! indxvrt3dVmix2
!
           call fillvalue3d(vrt3dVmix2,nciddiags_vrt3d,diags_vrt3dVmix2(itrc),
     &                      vname(1,indxvrt3dVmix2+itrc-1),
     &                      record,r3dvar,type)

# if defined DIAGNOSTICS_BARO
!
! indxvrtBaro
!
           work2d=vrtBaro(:,:)
           call fillvalue2d(work2d,nciddiags_vrt3d,diags_vrtBaro(itrc),
     &                      vname(1,indxvrtBaro+itrc-1),
     &                      record,ivar,type)
# endif

# if defined M3FAST
!
! indxvrt3dfast
!
           call fillvalue3d(vrt3dfast,nciddiags_vrt3d,diags_vrt3dfast(itrc),
     &        vname(1,indxvrt3dfast+itrc-1),record,r3dvar,type)
# endif

        endif

  1   format(/1x,'WRT_DIAGS_VRT3D ERROR while writing variable ''', A,
     &          ''' into diag file.', /11x, 'Time record:',
     &                   I6,3x,'netCDF error code',i4,3x,a,i4)
      goto 100
  99  may_day_flag=3
 100  continue

!
!---------------------------------------------------------------
! Synchronize netCDF file to disk to allow other processes
! to access data immediately after it is written.
!---------------------------------------------------------------
!
# if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
      ierr=nf_close (nciddiags_vrt3d)
      if (nrpfdiags_vrt3d.gt.0 .and. record.ge.nrpfdiags_vrt3d) nciddiags_vrt3d=-1
# else
      if (nrpfdiags_vrt3d.gt.0 .and. record.ge.nrpfdiags_vrt3d) then
        ierr=nf_close (nciddiags_vrt3d)
        nciddiags_vrt3d=-1
      else
        ierr=nf_sync(nciddiags_vrt3d)
      endif
# endif


      if (ierr .eq. nf_noerr) then
        MPI_master_only write(stdout,'(6x,A,2(A,I4,1x),A,I3)')
     &            'WRT_DIAGS_VRT3D -- wrote',
     &            ' diag fields into time record =', record, '/',
     &             nrecdiags_vrt3d  MYID
      else
        MPI_master_only write(stdout,'(/1x,2A/)')
     &             'WRT_DIAGS_VRT3D ERROR: Cannot ',
     &             'synchronize/close diag netCDF file.'
        may_day_flag=3
      endif

# if defined MPI & !defined PARALLEL_FILES & !defined NC4PAR
      if (mynode .lt. NNODES-1) then
        call MPI_Send (blank, 1, MPI_INTEGER, mynode+1,
     &                        1, MPI_COMM_WORLD,  ierr)
      endif
# endif
      return
      end

#else
      subroutine wrt_diags_vrt3d_empty
      end
#endif /* (DIAGNOSTICS_VRT) */
