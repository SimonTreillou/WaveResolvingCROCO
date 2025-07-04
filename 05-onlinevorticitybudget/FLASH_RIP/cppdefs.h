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
/*
   This is "cppdefs.h": MODEL CONFIGURATION FILE
   ==== == ============ ===== ============= ====
*/
/*
        SELECT ACADEMIC TEST CASES
*/
#undef  BASIN           /* Basin Example */
#undef  CANYON          /* Canyon Example */
#undef  EQUATOR         /* Equator Example  */
#undef  INNERSHELF      /* Inner Shelf Example */
#undef  SINGLE_COLUMN   /* 1DV vertical mixing Example */
#undef  RIVER           /* River run-off Example */
#undef  OVERFLOW        /* Gravitational/Overflow Example */
#undef  SEAMOUNT        /* Seamount Example */
#undef  SHELFRONT       /* Shelf Front Example */
#undef  SOLITON         /* Equatorial Rossby Wave Example */
#undef  THACKER         /* Thacker wetting-drying Example */
#undef  UPWELLING       /* Upwelling Example */
#undef  VORTEX          /* Baroclinic Vortex Example */
#undef  INTERNAL        /* Internal Tide Example */
#undef  IGW             /* COMODO Internal Tide Example */
#undef  JET             /* Baroclinic Jet Example */
#undef  SHOREFACE       /* Shoreface Test Case on a Planar Beach */
#undef  RIP             /* Rip Current Test Case */
#define  FLASH_RIP       /* Flash Rip Test Case */
#undef  SANDBAR         /* Bar-generating Flume Example */
#undef  SWASH           /* Swash Test Case on a Planar Beach */
#undef  TANK            /* Tank Example */
#undef  MOVING_BATHY    /* Moving Bathymetry Example */
#undef  ACOUSTIC        /* Acoustic wave Example */
#undef  GRAV_ADJ        /* Graviational Adjustment Example */
#undef  ISOLITON        /* Internal Soliton Example */
#undef  KH_INST         /* Kelvin-Helmholtz Instability Example */
#undef  TS_HADV_TEST    /* Horizontal tracer advection Example */
#undef  DUNE            /* Dune migration Example */
#undef  SED_TOY         /* 1DV sediment toy Example */
#undef  TIDAL_FLAT      /* 2DV tidal flat Example */
#undef  ESTUARY         /* 3D tidal estuary Example */
#undef  KILPATRICK      /* 2D sst front*/
#undef  SEAGRASS        /* 2DV over seagrass using OBSTRUCTION module*/
/*
        ... OR REALISTIC CONFIGURATIONS
*/
#undef  COASTAL         /* COASTAL Applications */
#undef REGIONAL        /* REGIONAL Applications */



#if defined REGIONAL
/*
!====================================================================
!               REGIONAL (realistic) Configurations
!====================================================================
!
!----------------------
! BASIC OPTIONS
!----------------------
!
*/
                      /* Configuration Name */
# define BENGUELA_LR
                      /* Parallelization */
# undef  OPENMP
# undef  MPI
                      /* Non-hydrostatic option */
# undef  NBQ
# undef  CROCO_QH
                      /* Nesting */
# undef  AGRIF
# undef  AGRIF_2WAY
                      /* OA and OW Coupling via OASIS (MPI) */
# undef  OA_COUPLING
# undef  OW_COUPLING
# ifdef OW_COUPLING
#  undef OW_COUPLING_FULL
#  undef WAVE_SMFLUX
# endif
                      /* Wave-current interactions */
# undef  MRL_WCI
                      /* Open Boundary Conditions */
# undef  TIDES
# define OBC_EAST
# define OBC_WEST
# define OBC_NORTH
# define OBC_SOUTH
                      /* Applications */
# undef  BIOLOGY
# undef  FLOATS
# undef  STATIONS
# undef  PASSIVE_TRACER
# undef  SEDIMENT
# undef  MUSTANG
# undef  BBL
                      /* I/O server */
# undef  XIOS
                      /* Calendar */
# undef  USE_CALENDAR
                      /* dedicated croco.log file */
# undef  LOGFILE
/*
!-------------------------------------------------
! PRE-SELECTED OPTIONS
!
! ADVANCED OPTIONS ARE IN CPPDEFS_DEV.H
!-------------------------------------------------
*/
                      /* Parallelization */
# ifdef MPI
#  undef  PARALLEL_FILES
#  undef  NC4PAR
#  undef  MPI_NOLAND
#  undef  MPI_TIME
# endif
# undef  AUTOTILING
                      /* Grid configuration */
# define CURVGRID
# define SPHERICAL
# define MASKING
# undef  WET_DRY
# define NEW_S_COORD
                      /* Model dynamics */
# define SOLVE3D
# define UV_COR
# define UV_ADV
                      /* Equation of State */
# define SALINITY
# define NONLIN_EOS
                      /* Lateral Forcing */
# undef CLIMATOLOGY
# ifdef CLIMATOLOGY
#  define ZCLIMATOLOGY
#  define M2CLIMATOLOGY
#  define M3CLIMATOLOGY
#  define TCLIMATOLOGY

#  define ZNUDGING
#  define M2NUDGING
#  define M3NUDGING
#  define TNUDGING
#  undef  ROBUST_DIAG
# endif

# define FRC_BRY
# ifdef FRC_BRY
#  define Z_FRC_BRY
#  define M2_FRC_BRY
#  define M3_FRC_BRY
#  define T_FRC_BRY
# endif
                      /* Surface Forcing */
/*
! Bulk flux algorithms (options)
! by default : COARE3p0 paramet with GUSTINESS effects
!
! To change bulk param, define one the following keys (exclusive) :
! - define BULK_ECUMEV0 : ECUME_v0 param
! - define BULK_ECUMEV6 : ECUME_v6 param
! - define BULK_WASP    : WASP param
! Note : gustiness effects can be added for all params
!        by defining BULK_GUSTINESS
*/
# define BULK_FLUX
# ifdef BULK_FLUX
#  undef  BULK_ECUMEV0
#  undef  BULK_ECUMEV6
#  undef  BULK_WASP
#  define BULK_GUSTINESS
#  define BULK_LW
#  undef  SST_SKIN
#  undef  ANA_DIURNAL_SW
#  undef  ONLINE
#  ifdef ONLINE
#   undef  AROME
#   undef  ERA_ECMWF
#  endif
#  undef READ_PATM
#  ifdef READ_PATM
#   define OBC_PATM
#  endif
#  undef  ABL1D
#  ifdef  ABL1D
#   undef  ANA_ABL_LSDATA
#   undef  ANA_ABL_VGRID
#   define STRESS_AT_RHO_POINTS
#   define ABL_NUDGING
#   define ABL_NUDGING_DYN
#   define ABL_NUDGING_TRA
#   undef  ABL_DYN_RESTORE_EQ
#   undef  SFLUX_CFB
#  endif
# else
#  define QCORRECTION
#  define SFLX_CORR
#  undef  SFLX_CORR_COEF
#  define ANA_DIURNAL_SW
# endif
# undef  SFLUX_CFB
# undef  SEA_ICE_NOFLUX
                      /* Lateral Momentum Advection (default UP3) */
# define UV_HADV_UP3
# undef  UV_HADV_UP5
# undef  UV_HADV_WENO5
# undef  UV_HADV_TVD
                      /* Lateral Explicit Momentum Mixing */
# undef  UV_VIS2
# ifdef UV_VIS2
#  define UV_VIS_SMAGO
# endif
                      /* Vertical Momentum Advection */
# define UV_VADV_SPLINES
# undef  UV_VADV_WENO5
# undef  UV_VADV_TVD
                      /* Lateral Tracer Advection (default UP3) */
# undef  TS_HADV_UP3
# define TS_HADV_RSUP3
# undef  TS_HADV_UP5
# undef  TS_HADV_WENO5
                      /* Lateral Explicit Tracer Mixing */
# undef  TS_DIF2
# undef  TS_DIF4
# undef  TS_MIX_S
                      /* Vertical Tracer Advection  */
# define TS_VADV_SPLINES
# undef  TS_VADV_AKIMA
# undef  TS_VADV_WENO5
                      /* Sponge layers for UV and TS */
# define SPONGE
                      /* Semi-implicit Vertical Tracer/Mom Advection */
# undef  VADV_ADAPT_IMP
                      /* Bottom friction in fast 3D step */
# define LIMIT_BSTRESS
# undef  BSTRESS_FAST
                      /* Vertical Mixing */
# undef  BODYFORCE
# undef  BVF_MIXING
# define LMD_MIXING
# undef  GLS_MIXING
# ifdef LMD_MIXING
#  define LMD_SKPP
#  define LMD_BKPP
#  define LMD_RIMIX
#  define LMD_CONVEC
#  define LMD_NONLOCAL
#  undef  LMD_DDMIX
#  undef  LMD_LANGMUIR
# endif
                      /* Wave-current interactions */
# ifdef OW_COUPLING
#  define MRL_WCI
#  define BBL
# endif
# ifdef MRL_WCI
#  ifndef OW_COUPLING
#   undef  WAVE_OFFLINE
#   define ANA_WWAVE
#   undef  WKB_WWAVE
#  endif
#  undef  WAVE_ROLLER
#  define WAVE_STREAMING
#  define WAVE_FRICTION
#  define WAVE_RAMP
#  ifdef WKB_WWAVE
#   undef  WKB_OBC_NORTH
#   undef  WKB_OBC_SOUTH
#   define WKB_OBC_WEST
#   undef  WKB_OBC_EAST
#  endif
# endif
                      /* Bottom Forcing */
# define ANA_BSFLUX
# define ANA_BTFLUX
                      /* Point Sources - Rivers */
# undef PSOURCE
# undef PSOURCE_NCFILE
# ifdef PSOURCE_NCFILE
#  undef PSOURCE_NCFILE_TS
# endif
                      /* Open Boundary Conditions */
# ifdef TIDES
#  define SSH_TIDES
#  define UV_TIDES
#  define POT_TIDES
#  undef  TIDES_MAS
#  define TIDERAMP
# endif
# define OBC_M2CHARACT
# undef  OBC_M2ORLANSKI
# define OBC_M3ORLANSKI
# define OBC_TORLANSKI
# undef  OBC_M2SPECIFIED
# undef  OBC_M3SPECIFIED
# undef  OBC_TSPECIFIED
                      /* Input/Output */
# define AVERAGES
# define AVERAGES_K
# undef  OUTPUTS_SURFACE
# undef  HOURLY_VELOCITIES
                     /* Exact restart */
# undef EXACT_RESTART
/*
!                        Diagnostics
!--------------------------------------------
! 3D Tracer & momentum balance
! 2D Mixing layer balance
! Depth-mean vorticity and energy balance
! Eddy terms
!--------------------------------------------
!
*/
# undef DO_NOT_OVERWRITE
# undef RESTART_DIAGS

# undef DIAGNOSTICS_TS
# undef DIAGNOSTICS_UV
# ifdef DIAGNOSTICS_TS
#  undef  DIAGNOSTICS_TS_ADV
#  undef  DIAGNOSTICS_TS_MLD
#  ifdef DIAGNOSTICS_TS_MLD
#   define DIAGNOSTICS_TS_MLD_CRIT
#  endif

# endif

# undef DIAGNOSTICS_TSVAR
# ifdef DIAGNOSTICS_TSVAR
#  define DIAGNOSTICS_TS
#  define DIAGNOSTICS_TS_ADV
# endif

# undef  DIAGNOSTICS_VRT
# undef  DIAGNOSTICS_EK
# ifdef DIAGNOSTICS_EK
#  undef DIAGNOSTICS_EK_FULL
#  undef DIAGNOSTICS_EK_MLD
# endif

# undef DIAGNOSTICS_BARO
# undef DIAGNOSTICS_PV
# undef DIAGNOSTICS_DISS
# ifdef DIAGNOSTICS_DISS
#  define DIAGNOSTICS_PV
# endif

# undef DIAGNOSTICS_EDDY

# undef TENDENCY
# ifdef TENDENCY
#  define DIAGNOSTICS_UV
# endif
/*
!           Applications:
!---------------------------------
! Biology, floats, Stations,
! Passive tracer, Sediments, BBL
!---------------------------------
!
   Quasi-monotone lateral advection scheme (WENO5)
   for passive/biology/sediment tracers
*/
# if defined PASSIVE_TRACER || defined BIOLOGY || defined SEDIMENT \
                                               || defined MUSTANG
#  define BIO_HADV_WENO5
# endif
                      /*   Choice of Biology models   */
# ifdef BIOLOGY
#  undef  PISCES
#  undef  BIO_NChlPZD
#  undef  BIO_N2ChlPZD2
#  define BIO_BioEBUS
                      /*   Biology options    */
#  ifdef PISCES
#   undef  DIURNAL_INPUT_SRFLX
#   define key_pisces
#   define key_ligand
#   undef key_pisces_quota
#   undef key_pisces_npzd
#   undef key_sediment
#  endif
#  ifdef BIO_NChlPZD
#   define OXYGEN
#  endif
#  ifdef BIO_BioEBUS
#   define NITROUS_OXIDE
#  endif
                      /*   Biology diagnostics    */
#  define DIAGNOSTICS_BIO
#  if defined DIAGNOSTICS_BIO && defined PISCES
#   define key_trc_diaadd
#  endif
# endif
                      /*   Lagrangian floats model    */
# ifdef FLOATS
#  undef  FLOATS_GLOBAL_ATTRIBUTES
#  undef  IBM
#  undef  RANDOM_WALK
#  ifdef RANDOM_WALK
#   define DIEL_MIGRATION
#   define RANDOM_VERTICAL
#   define RANDOM_HORIZONTAL
#  endif
# endif
                      /*   Stations recording    */
# ifdef STATIONS
#  define ALL_SIGMA
# endif
                      /*   USGS Sediment model     */
# ifdef SEDIMENT
#  define SUSPLOAD
#  define BEDLOAD
#  define MORPHODYN
# endif
                      /*   MUSTANG Sediment model     */
# ifdef MUSTANG
#  undef  key_MUSTANG_V2
#  undef  key_MUSTANG_bedload
#  undef  MORPHODYN
#  define key_sand2D
#  define MUSTANG_CORFLUX
#  undef  key_tauskin_c_upwind
#  undef  WAVE_OFFLINE
# endif

#elif defined COASTAL
/*
!====================================================================
!               COASTAL (realistic) Configurations
!====================================================================
!
!----------------------
! BASIC OPTIONS
!----------------------
!
*/
                      /* Configuration Name */
# define VILAINE
                      /* Parallelization */
# undef  OPENMP
# undef  MPI
                      /* Open Boundary Conditions */
# define TIDES
# undef  OBC_EAST
# define OBC_WEST
# undef  OBC_NORTH
# define OBC_SOUTH
                      /* Applications */
# undef  BIOLOGY
# undef  FLOATS
# undef  STATIONS
# undef  PASSIVE_TRACER
# undef  SEDIMENT
# undef  BBL
# define MUSTANG
                      /* I/O server */
# undef  XIOS
                     /* Custion IO */
# define ZETA_DRY_IO
# define FILLVAL
                      /* Calendar */

# define USE_CALENDAR
                      /* dedicated croco.log file */
# undef  LOGFILE
/*!
!-------------------------------------------------
! PRE-SELECTED OPTIONS
!
! ADVANCED OPTIONS ARE IN CPPDEFS_DEV.H
!-------------------------------------------------
*/
                      /* Parallelization */
# ifdef MPI
#  define NC4PAR
#  undef  MPI_NOLAND
#  undef  MPI_TIME
# endif
# undef  AUTOTILING
                      /* Non-hydrostatic options */
# ifdef NBQ
#  define W_HADV_WENO5
#  define W_VADV_WENO5
# endif
                      /* Grid configuration */
# define ANA_INITIAL
# define CURVGRID
# define SPHERICAL
# define MASKING
# define WET_DRY
# define NEW_S_COORD
                      /* Model dynamics */
# define SOLVE3D
# define UV_COR
# define UV_ADV
                      /* Equation of State */
# define SALINITY
# define NONLIN_EOS
                      /* Lateral Momentum Advection (default UP3) */
# undef  UV_HADV_UP3
# define UV_HADV_WENO5
                      /* Lateral Explicit Momentum Mixing */
# define UV_VIS2
# ifdef UV_VIS2
#  define UV_VIS_SMAGO
# endif
                      /* Vertical Momentum Advection  */
# undef  UV_VADV_SPLINES
# define UV_VADV_WENO5
                      /* Lateral Tracer Advection (default UP3) */
# undef  TS_HADV_UP3
# define TS_HADV_WENO5
                      /* Lateral Explicit Tracer Mixing */
# define TS_DIF2
# define TS_MIX_S
                      /* Vertical Tracer Advection  */
# undef  TS_VADV_SPLINES
# define TS_VADV_WENO5
                      /* Sponge layers for UV and TS */
# define SPONGE
                      /* Semi-implicit Vertical Tracer/Mom Advection */
# undef  VADV_ADAPT_IMP
                      /* Bottom friction in fast 3D step */
# define LIMIT_BSTRESS
# undef  BSTRESS_FAST
                      /* Vertical Mixing */
# define GLS_MIXING
                      /* Surface Forcing */
# define BULK_FLUX
# ifdef BULK_FLUX
#  undef  ECUMEv0
#  undef  ECUMEv6
#  undef  WASP
#  define GUSTINESS
#  undef  BULK_LW
#  undef  SST_SKIN
#  undef  ANA_DIURNAL_SW
#  define ONLINE
#  ifdef ONLINE
#   define AROME
#   undef  ERA_ECMWF
#  endif
#  define READ_PATM
#  ifdef READ_PATM
#   define OBC_PATM
#  endif
# else
#  undef  QCORRECTION
#  undef  SFLX_CORR
#  undef  SFLX_CORR_COEF
#  undef  ANA_DIURNAL_SW
# endif
# undef  ANA_SSFLUX
# define ANA_STFLUX
                      /* Lateral Forcing */
# undef  ANA_BRY
# define FRC_BRY
# ifdef FRC_BRY
#  define Z_FRC_BRY
#  define M2_FRC_BRY
#  undef  M3_FRC_BRY
#  define T_FRC_BRY
# endif
                      /* Bottom Forcing */
# define ANA_BSFLUX
# define ANA_BTFLUX
                      /* Point Sources - Rivers */
# define PSOURCE
# undef  PSOURCE_MASS
# define PSOURCE_NCFILE
# ifdef PSOURCE_NCFILE
#  define PSOURCE_NCFILE_TS
# endif
                      /* Open Boundary Conditions */
# ifdef TIDES
#  define M2FILTER_NONE
#  define SSH_TIDES
#  define UV_TIDES
#  undef  POT_TIDES
#  define TIDES_MAS
#  ifndef UV_TIDES
#   define OBC_REDUCED_PHYSICS
#  endif
#  define TIDERAMP
# endif
# define OBC_M2CHARACT
# define OBC_M3ORLANSKI
# define OBC_TORLANSKI
                      /* Input/Output */
# undef  AVERAGES
# undef  AVERAGES_K
# undef  OUTPUTS_SURFACE
# undef  HOURLY_VELOCITIES
/*
!           Applications:
!---------------------------------
! Biology, floats, Stations,
! Passive tracer, Sediments, BBL
!---------------------------------
!
   Quasi-monotone lateral advection scheme (WENO5)
   for passive/biology/sediment tracers
*/
# if defined PASSIVE_TRACER || defined BIOLOGY || defined SEDIMENT \
                            || defined SUBSTANCE || defined MUSTANG
#  define BIO_HADV_WENO5
# endif
                      /*     USGS Sediment model     */
# ifdef SEDIMENT
#  define SUSPLOAD
#  define BEDLOAD
#  define MORPHODYN
# endif
                      /*   MUSTANG Sediment model     */
# ifdef MUSTANG
#  undef  key_MUSTANG_V2
#  undef  key_MUSTANG_bedload
#  undef  MORPHODYN
#  define key_sand2D
#  define MUSTANG_CORFLUX
#  undef  key_tauskin_c_upwind
#  define WAVE_OFFLINE
# endif
/*
!
!==========================================================
!              IDEALIZED CONFIGURATIONS
!==========================================================
!
*/
#elif defined BASIN
/*
!                       Basin Example
!                       ===== =======
*/
# undef  OPENMP
# undef  MPI
# define UV_ADV
# define UV_COR
# define UV_VIS2
# define SOLVE3D
# define TS_DIF2
# define BODYFORCE
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_BTFLUX
# define NO_FRCFILE

#elif defined CANYON
/*
!                       Canyon Example
!                       ====== =======
*/
# undef  OPENMP
# undef  MPI
# define CANYON_STRAT
# define UV_ADV
# define UV_COR
# define SOLVE3D
# define EW_PERIODIC
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_BTFLUX
# define NO_FRCFILE

#elif defined EQUATOR
/*
!                       Equator Example
!                       ======= =======
! Boccaletti, G., R.C. Pacanowski, G.H. Philander and A.V. Fedorov, 2004,
! The Thermal Structure of the Upper Ocean, J.Phys.Oceanogr., 34, 888-902.
*/
# undef  OPENMP
# undef  MPI
# define UV_ADV
# define UV_COR
# define UV_VIS2
# define SOLVE3D
# define TS_DIF2
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SRFLUX
# define ANA_SSFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define QCORRECTION
# define ANA_SST
# define LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define NO_FRCFILE

#elif defined INNERSHELF
/*
!                       Inner Shelf Example
!                       ===== ===== =======
*/
# undef  OPENMP
# undef  MPI
# undef  NBQ
# define INNERSHELF_EKMAN
# define INNERSHELF_APG
# define SOLVE3D
# define UV_COR
# define ANA_GRID
# define ANA_INITIAL
# define AVERAGES
# define ANA_SSFLUX
# define ANA_SRFLUX
# define ANA_STFLUX
# define ANA_BSFLUX
# define ANA_BTFLUX
# define ANA_SMFLUX
# define NS_PERIODIC
# define OBC_WEST
# define SPONGE
# ifndef INNERSHELF_EKMAN
#  define UV_ADV
#  define SALINITY
#  define NONLIN_EOS
#  define LMD_MIXING
#  undef  GLS_MIXING
#  ifdef LMD_MIXING
#   define LMD_SKPP
#   define LMD_BKPP
#   define LMD_RIMIX
#   define LMD_CONVEC
#  endif
#  undef WAVE_MAKER_INTERNAL
#  ifdef WAVE_MAKER_INTERNAL
#   define ANA_BRY
#   define Z_FRC_BRY
#   define M2_FRC_BRY
#   define M3_FRC_BRY
#   define T_FRC_BRY
#  endif
# endif
# define NO_FRCFILE

#elif defined SINGLE_COLUMN
/*
!                       Single Column Example
!                       ====== ====== =======
!
!                              Seven  sets up are encompassed :
*/
/* erosion of linear strat by constant wind stress */
# define KATO_PHILIPS
/* erosion of linear strat by constant surf buoyancy loss */
# undef  WILLIS_DEARDORFF
/* erosion of linear strat by constant surf buoyancy loss */
# undef  DIURNAL_CYCLE
/* forced Ekman bottom boundary layer */
# undef  FORCED_EKBBL
/* forced Ekman bottom and surface boundary layers */
# undef  FORCED_DBLEEK
/* non rotating forced bottom boundary layer : Prandt layer */
# undef  FORCED_NONROTBBL
/* non rotating oscillatory forced bottom boundary layer */
# undef  FORCED_OSCNONROTBBL
# undef  OPENMP
# undef  MPI
# define UV_ADV
# define NEW_S_COORD
# define UV_COR
# define SOLVE3D
# undef  LMD_MIXING
# define GLS_MIXING
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_SRFLUX
# define ANA_STFLUX
# define ANA_BTFLUX
# define EW_PERIODIC
# define NS_PERIODIC

#elif defined INTERNAL
/*
!                       Internal Tide Example
!                       ======== ==== =======
!
! Di Lorenzo, E, W.R. Young and S.L. Smith, 2006, Numerical and anlytical estimates of M2
! tidal conversion at steep oceanic ridges, J. Phys. Oceanogr., 36, 1072-1084.
*/
# undef  OPENMP
# undef  MPI
# define SOLVE3D
# define UV_COR
# define UV_ADV
# define BODYTIDE
# define ANA_GRID
# define ANA_INITIAL
# define ANA_BTFLUX
# define ANA_SMFLUX
# define ANA_SRFLUX
# define ANA_STFLUX
# define ANA_VMIX
# define EW_PERIODIC
# define NS_PERIODIC
# ifdef INTERNALSHELF
#  undef   EW_PERIODIC
#  define OBC_EAST
#  define OBC_WEST
#  define SPONGE
#  define ANA_SSH
#  define ANA_M2CLIMA
#  define ANA_M3CLIMA
#  define ANA_TCLIMA
#  define ZCLIMATOLOGY
#  define M2CLIMATOLOGY
#  define M3CLIMATOLOGY
#  define TCLIMATOLOGY
#  define M2NUDGING
#  define M3NUDGING
#  define TNUDGING
# endif
# define NO_FRCFILE

#elif defined IGW
/*
!                  COMODO Internal Tide Example
!                  ====== ======== ==== =======
!
! Pichon, A., 2007: Tests academiques de maree, Rapport interne n 21 du 19 octobre 2007,
! Service Hydrographique et Oceanographique de la Marine.
*/

# define EXPERIMENT3
# undef  OPENMP
# undef  MPI
# undef  NBQ
# define NEW_S_COORD
# define TIDES
# define TIDERAMP
# define SSH_TIDES
# define UV_TIDES
# define SOLVE3D
# define UV_ADV
# define UV_COR
# define UV_VIS2
# undef  VADV_ADAPT_IMP
# define SPHERICAL
# define CURVGRID
# define ANA_INITIAL
# define ANA_VMIX
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SRFLUX
# define ANA_SSFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define NS_PERIODIC
# define OBC_EAST
# define OBC_WEST
# undef  SPONGE
# define ANA_SSH
# define ANA_M2CLIMA
# define ANA_M3CLIMA
# define ANA_TCLIMA
# define ZCLIMATOLOGY
# define M2CLIMATOLOGY
# define M3CLIMATOLOGY
# define TCLIMATOLOGY
# define M2NUDGING
# define M3NUDGING
# define TNUDGING
# undef  ONLINE_ANALYSIS

#elif defined RIVER
/*
!                       River run-off test problem
!                       ==========================
*/
# undef  OPENMP
# undef  MPI
# define SOLVE3D
# define UV_ADV
# define UV_COR
# define NONLIN_EOS
# define SALINITY
# define ANA_GRID
# define MASKING
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_SRFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define LMD_MIXING
# define LMD_SKPP
# define LMD_BKPP
# define LMD_RIMIX
# define LMD_CONVEC
# define PSOURCE
# undef  PSOURCE_MASS
# define ANA_PSOURCE
# define NS_PERIODIC
# undef  FLOATS
# ifdef FLOATS
#  define RANDOM_WALK
#  ifdef RANDOM_WALK
#   define DIEL_MIGRATION
#   define RANDOM_VERTICAL
#   define RANDOM_HORIZONTAL
#  endif
# endif
# define NO_FRCFILE

#elif defined SEAMOUNT
/*
!                       Seamount Example
!                       ======== =======
*/
# undef  OPENMP
# undef  MPI
# define UV_ADV
# define UV_COR
# define SOLVE3D
# define SALINITY
# define NONLIN_EOS
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_SRFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define NO_FRCFILE

#elif defined SHELFRONT
/*
!                       Shelf Front Example
!                       ===== ===== =======
*/
# undef  OPENMP
# undef  MPI
# define UV_ADV
# define UV_COR
# define SOLVE3D
# define SALINITY
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_SRFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define EW_PERIODIC
# define NO_FRCFILE

#elif defined SOLITON
/*
!                       Equatorial Rossby Wave Example
!                       ========== ====== ==== =======
*/
# undef  OPENMP
# undef  MPI
# define UV_COR
# define UV_ADV
# define ANA_GRID
# define ANA_INITIAL
# define AVERAGES
# define EW_PERIODIC
# define ANA_SMFLUX
# define ANA_BTFLUX
# define NO_FRCFILE

#elif defined THACKER
/*
!                       Thacker Example
!                       ======= =======
!
! Thacker, W., (1981), Some exact solutions to the nonlinear
! shallow-water wave equations. J. Fluid Mech., 107, 499-508.
*/
# undef  OPENMP
# undef  MPI
# define THACKER_2DV
# define SOLVE3D
# define UV_COR
# define UV_ADV
# undef  UV_VIS2
# define WET_DRY
# define NEW_S_COORD
# define ANA_GRID
# define ANA_INITIAL
# define ANA_BTFLUX
# define ANA_SMFLUX
# define ANA_SRFLUX
# define ANA_STFLUX
# define NO_FRCFILE

#elif defined OVERFLOW
/*
!                       Gravitational/Overflow Example
!                       ====================== =======
*/
# undef  OPENMP
# undef  MPI
# define UV_ADV
# define UV_COR
# define UV_VIS2
# define TS_DIF2
# define TS_MIX_GEO
# define SOLVE3D
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_BTFLUX
# define NO_FRCFILE

#elif defined UPWELLING
/*
!                       Upwelling Example
!                       ========= =======
*/
# undef  OPENMP
# undef  MPI
# define SOLVE3D
# define UV_COR
# define UV_ADV
# define ANA_GRID
# define ANA_INITIAL
# define AVERAGES
# define SALINITY
# define NONLIN_EOS
# define ANA_SSFLUX
# define ANA_SRFLUX
# define ANA_STFLUX
# define ANA_BSFLUX
# define ANA_BTFLUX
# define ANA_SMFLUX
# define LMD_MIXING
# define LMD_SKPP
# define LMD_BKPP
# define LMD_RIMIX
# define LMD_CONVEC
# define EW_PERIODIC
# define NO_FRCFILE

#elif defined VORTEX
/*
!                       Baroclinic Vortex Example (TEST AGRIF)
!                       ========== ====== ======= ===== ======
*/
# undef  OPENMP
# undef  MPI
# undef  AGRIF
# undef  AGRIF_2WAY
# undef  NBQ
# define SOLVE3D
# define UV_COR
# define UV_ADV
# define ANA_STFLUX
# define ANA_SMFLUX
# define ANA_BSFLUX
# define ANA_BTFLUX
# define ANA_VMIX
# define OBC_EAST
# define OBC_WEST
# define OBC_NORTH
# define OBC_SOUTH
# define SPONGE
# define ZCLIMATOLOGY
# define M2CLIMATOLOGY
# define M3CLIMATOLOGY
# define TCLIMATOLOGY
# define ZNUDGING
# define M2NUDGING
# define M3NUDGING
# define TNUDGING
# define NO_FRCFILE

#elif defined JET
/*
!                       Baroclinic JET Example
!                       ========== === =======
*/
# define ANA_JET
# undef  MPI
# undef  NBQ
# define SOLVE3D
# define UV_COR
# define UV_ADV
# define UV_VIS2
# ifdef ANA_JET
#  define ANA_GRID
#  define ANA_INITIAL
# endif
# define ANA_STFLUX
# define ANA_SMFLUX
# define ANA_BSFLUX
# define ANA_BTFLUX
# define ANA_VMIX
# define EW_PERIODIC
# define CLIMATOLOGY
# ifdef CLIMATOLOGY
#  define ZCLIMATOLOGY
#  define M2CLIMATOLOGY
#  define M3CLIMATOLOGY
#  define TCLIMATOLOGY
#  define ZNUDGING
#  define M2NUDGING
#  define M3NUDGING
#  define TNUDGING
#  define ROBUST_DIAG
#  define ZONAL_NUDGING
#  ifdef ANA_JET
#   define ANA_SSH
#   define ANA_M2CLIMA
#   define ANA_M3CLIMA
#   define ANA_TCLIMA
#  endif
# endif
# define LMD_MIXING
# ifdef  LMD_MIXING
#  undef  ANA_VMIX
#  define ANA_SRFLUX
#  undef  LMD_KPP
#  define LMD_RIMIX
#  define LMD_CONVEC
# endif
# define NO_FRCFILE

#elif defined SHOREFACE
/*
!                       PLANAR BEACH Example
!                       ====== ===== =======
!
!   Uchiyama, Y., McWilliams, J.C. and Shchepetkin, A.F. (2010):
!      Wave-current interaction in an oceanic circulation model with a
!      vortex force formalism: Application to the surf zone.
!      Ocean Modelling Vol. 34:1-2, pp.16-35.
*/
# undef  OPENMP
# undef  MPI
# define SOLVE3D
# define UV_ADV
# undef  MASKING
# define WET_DRY
# define NEW_S_COORD
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_SRFLUX
# define ANA_SST
# define ANA_BTFLUX
# define NS_PERIODIC
# define OBC_WEST
# define SPONGE
# define MRL_WCI
# ifdef MRL_WCI
#  undef  WAVE_OFFLINE
#  ifndef WAVE_OFFLINE
#   define WKB_WWAVE
#   define WKB_OBC_WEST
#   define WAVE_FRICTION
#   undef  WAVE_ROLLER
#   undef  MRL_CEW
#  endif
# endif
# define LMD_MIXING
# define LMD_SKPP
# define LMD_BKPP
# undef  BBL
# undef  SEDIMENT
# ifdef SEDIMENT
#  define SUSPLOAD
#  define BEDLOAD
#  define TCLIMATOLOGY
#  define TNUDGING
#  define ANA_TCLIMA
# endif

#elif defined SANDBAR
/*
!                       SANDBAR Example
!                       ======= =======
!
!   Roelvink, J. A. and Reniers, A. (1995). Lip 11d delta flume experiments
!   – data report. Technical report, Delft, The Netherlands, Delft Hydraulics
*/
/* LIP-1B */
# define SANDBAR_OFFSHORE
/* LIP-1C */
# undef  SANDBAR_ONSHORE
# undef  OPENMP
# undef  MPI
# undef  NBQ
# define SOLVE3D
# define UV_ADV
# define NEW_S_COORD
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_SRFLUX
# define ANA_SST
# define ANA_BTFLUX
# define OBC_WEST
# define SPONGE
# define WET_DRY
/* ! NBQ */
# ifndef NBQ 
#  define MRL_WCI
#  ifdef MRL_WCI
#   define WKB_WWAVE
#   define MRL_CEW
#   define WKB_OBC_WEST
#   define WAVE_ROLLER
#   define WAVE_FRICTION
#   define WAVE_BREAK_TG86
#   define WAVE_STREAMING
#   undef  WAVE_RAMP
#  endif
#  define GLS_MIXING
#  define GLS_KOMEGA
#  undef  LMD_MIXING
#  ifdef LMD_MIXING
#   define LMD_SKPP
#   define LMD_BKPP
#  endif
#  define BBL
#  define BBL_BREAKING_STIR
/* NBQ */
# else
#  define MPI
#  define NBQ_PRECISE
#  define WAVE_MAKER
#  define UV_ADV
#  define UV_HADV_WENO5
#  define UV_VADV_WENO5
#  define W_HADV_WENO5
#  define W_VADV_WENO5
#  define GLS_MIXING_3D
#  define GLS_KOMEGA
#  define ANA_BRY
#  define Z_FRC_BRY
#  define M2_FRC_BRY
#  define M3_FRC_BRY
#  define T_FRC_BRY
#  define AVERAGES
#  define AVERAGES_K
#  define DIAGNOSTICS_EDDY
/* NBQ */
# endif 
# define SEDIMENT
# ifdef SEDIMENT
#  define SUSPLOAD
#  ifndef NBQ
#   define BEDLOAD
#  endif
#  define MORPHODYN
#  define TCLIMATOLOGY
#  define TNUDGING
#  define ANA_TCLIMA
# endif
# undef  STATIONS
# ifdef STATIONS
#  define ALL_SIGMA
# endif
# undef  DIAGNOSTICS_TS
# ifdef DIAGNOSTICS_TS
#  define DIAGNOSTICS_TS_ADV
# endif
# define NO_FRCFILE

#elif defined RIP
/*
!                       Rip Current Example
!                       === ======= =======
!
!   Weir, B., Uchiyama, Y. (2010):
!      A vortex force analysis of the interaction of rip
!      currents and surface gravity wave
!      JGR Vol. 116
!
!  Default is an idealized Duck Beach wave-averaged simulation
!  with 3D topography
!
!  RIP_TOPO_2D: Logshore uniform topography
!  BISCA: realistic case with Input files
!  GRANDPOPO: idealized Grand Popo Beach in Benin
!  WAVE_MAKER & NBQ : wave-resolving (#undef MRL_WCI)
*/
# undef  BISCA
# undef  RIP_TOPO_2D
# undef  GRANDPOPO
# ifdef GRANDPOPO
#  define RIP_TOPO_2D
# endif
# undef  ANA_TIDES
# undef  OPENMP
# undef  MPI
# define SOLVE3D
# define NEW_S_COORD
# define UV_ADV
# undef  NBQ
# ifdef NBQ
#  define NBQ_PRECISE
#  define WAVE_MAKER
#  define WAVE_MAKER_SPECTRUM
#  define WAVE_MAKER_DSPREAD
#  define UV_HADV_WENO5
#  define UV_VADV_WENO5
#  define W_HADV_WENO5
#  define W_VADV_WENO5
#  define GLS_MIXING_3D
#  define GLS_KOMEGA
#  undef  ANA_TIDES
#  undef  MRL_WCI
#  define OBC_SPECIFIED_WEST
#  define FRC_BRY
#  define ANA_BRY
#  define Z_FRC_BRY
#  define M2_FRC_BRY
#  define M3_FRC_BRY
#  define T_FRC_BRY
#  define AVERAGES
#  define AVERAGES_K
# else
#  define UV_VIS2
#  define UV_VIS_SMAGO
#  define LMD_MIXING
#  define LMD_SKPP
#  define LMD_BKPP
#  define MRL_WCI
# endif
# ifdef MRL_WCI
#  define WKB_WWAVE
#  define WKB_OBC_WEST
#  define WAVE_ROLLER
#  define WAVE_FRICTION
#  define WAVE_STREAMING
#  define MRL_CEW
#  ifdef RIP_TOPO_2D
#   define WAVE_RAMP
#  endif
# endif
# ifndef BISCA
#  define ANA_GRID
# endif
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_SRFLUX
# define ANA_SST
# define ANA_BTFLUX
# define OBC_WEST
# define SPONGE
# define WET_DRY
# if !defined BISCA && !defined ANA_TIDES
#  define NS_PERIODIC
# else
#  define OBC_NORTH
#  define OBC_SOUTH
# endif
# ifdef ANA_TIDES
#  define ANA_SSH
#  define ANA_M2CLIMA
#  define ANA_M3CLIMA
#  define ZCLIMATOLOGY
#  define M2CLIMATOLOGY
#  define M3CLIMATOLOGY
#  define M2NUDGING
#  define M3NUDGING
# endif
# ifdef BISCA
#  define BBL
# endif
# undef SEDIMENT
# ifdef SEDIMENT
#  define SUSPLOAD
#  define BEDLOAD
#  undef  MORPHODYN
# endif
# undef  DIAGNOSTICS_UV

#elif defined FLASH_RIP
/*
!                       Flash Rip Example
!                       ===== === =======
!
!   Semi-idealized case based on the IB09 experiment.
!
!   Hally-Rosendahl, K., Feddersen, F. (2016):
!      Modeling surfzone to inner-shelf tracer exchange
!      JGR Vol. 121
!   Treillou, S. et al. (in review for JPO, 2025):
!      Tracer dispersion by surfzone eddies: assessing 
!      the impact of undertow vertical shear
*/
# define  MPI
# undef  NC4PAR
# define SOLVE3D
# define NEW_S_COORD
# define UV_ADV
# define NO_TRACER
# define NO_TEMPERATURE
# undef  PASSIVE_TRACER
# undef  ANA_TIDES
# define NBQ
# define NBQ_PRECISE
# define LIMIT_BSTRESS
# define WAVE_MAKER
# define WAVE_MAKER_SPECTRUM
# define WAVE_MAKER_DSPREAD
# define UV_HADV_WENO5
# define UV_VADV_WENO5
# define W_HADV_WENO5
# define W_VADV_WENO5
# define GLS_MIXING_3D
# define GLS_KOMEGA
# define NS_PERIODIC
# define OBC_WEST
# define OBC_SPECIFIED_WEST
# define SPONGE
# define FRC_BRY
# define ANA_BRY
# define Z_FRC_BRY
# define M2_FRC_BRY
# define M3_FRC_BRY
# define T_FRC_BRY
# define WET_DRY
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_SRFLUX
# define ANA_SST
# define ANA_BTFLUX
# define AVERAGES
# define AVERAGES_K
# define DIAGNOSTICS_EDDY
# define DIAGNOSTICS_VRT

#elif defined SWASH
/*
!                       SWASH PLANAR BEACH Example
!                       ===== ====== ===== =======
!
!  Simulation of GLOBEX experiments.
!  Marchesiello et al. (Ocean Modelling 2021):
!  
*/
# define SWASH_GLOBEX_B3
# undef  SWASH_GLOBEX_B2
# undef  SWASH_GLOBEX_A3
# undef  OPENMP
# undef  MPI
# define SOLVE3D
# define AVERAGES
# define NBQ
# define NBQ_PRECISE
# define WAVE_MAKER
# define UV_ADV
# define UV_HADV_WENO5
# define UV_VADV_WENO5
# define W_HADV_WENO5
# define W_VADV_WENO5
# define GLS_MIXING_3D
# define GLS_KOMEGA
# define NEW_S_COORD
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_SRFLUX
# define ANA_SST
# define ANA_BTFLUX
# define OBC_WEST
# undef  OBC_SPECIFIED_WEST
# define ANA_BRY
# define Z_FRC_BRY
# define M2_FRC_BRY
# define M3_FRC_BRY
# define T_FRC_BRY
# define WET_DRY
# define NO_FRCFILE

#elif defined TANK
/*
!                       Tank Example
!                       ==== =======
!
! Chen, X.J., 2003. A fully hydrodynamic model for three-dimensional,
! free-surface flows.
! Int. J. Numer. Methods Fluids 42, 929â~@~S952.
*/
# undef  MPI
# define NBQ
# define NBQ_PRECISE
# define SOLVE3D
# undef  UV_ADV
# define NEW_S_COORD
# define ANA_GRID
# define ANA_INITIAL
# define ANA_BTFLUX
# define ANA_SMFLUX
# define ANA_SRFLUX
# define ANA_STFLUX
# define NO_FRCFILE

#elif defined MOVING_BATHY
/*
!                       Moving Bathy Example
!                       ====== ===== =======
!
!  Auclair et al., Ocean Mod. 2014: Implementation of a time-dependent
!     bathymetry in a free-surface ocean model: Application to internal
!     wave generation
!
*/
# undef  MPI
# define ANA_MORPHODYN
# define NBQ
# define NBQ_PRECISE
# define SOLVE3D
# define NEW_S_COORD
# undef  PASSIVE_TRACER
# define UV_ADV
# define TS_HADV_WENO5
# define TS_VADV_WENO5
# define UV_HADV_WENO5
# define UV_VADV_WENO5
# define W_HADV_WENO5
# define W_VADV_WENO5
# define ANA_GRID
# define ANA_INITIAL
# define ANA_VMIX
# define ANA_BTFLUX
# define ANA_SMFLUX
# define ANA_SRFLUX
# define ANA_STFLUX
# define NO_FRCFILE

#elif defined ACOUSTIC
/*
!                       ACOUSTIC WAVE TESTCASE
!                       ======== ==== ========
*/
# undef  MPI
# define NBQ
# ifdef NBQ
#  undef  NBQ_PRECISE
# endif
# undef  UV_VIS2
# define SOLVE3D
# define NEW_S_COORD
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SRFLUX
# define ANA_BTFLUX
# define NO_FRCFILE

#elif defined GRAV_ADJ
/*
!                       Gravitational Adjustment Example
!                       ============= ========== =======
*/
# undef  OPENMP
# undef  MPI
# undef  NBQ
# undef  XIOS
# define SOLVE3D
# define NEW_S_COORD
# define UV_ADV
# define TS_HADV_WENO5
# define TS_VADV_WENO5
# define UV_HADV_WENO5
# define UV_VADV_WENO5
# ifdef NBQ
#  define W_HADV_WENO5
#  define W_VADV_WENO5
# endif
# undef  UV_VIS2
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_BTFLUX
# undef  PASSIVE_TRACER
# define NO_FRCFILE

#elif defined ISOLITON
/*
!                       Gravitational Adjustment Example
!                       ============= ========== =======
!
!  Internal soliton case ISOLITON (non-hydrostatic) is setup from:
!  Horn, D.A., J. Imberger, & G.N. Ivey, (2001).
!  The degeneration of large-scale interfacial gravity waves in lakes.
!  J. Fluid Mech., 434:181-207.
!
*/
# undef  MPI
# define NBQ
# undef  XIOS
# define SOLVE3D
# define NEW_S_COORD
# define UV_ADV
# define TS_HADV_WENO5
# define TS_VADV_WENO5
# define UV_HADV_WENO5
# define UV_VADV_WENO5
# define W_HADV_WENO5
# define W_VADV_WENO5
# undef  UV_VIS2
# undef  TS_DIF2
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_BTFLUX
# undef  PASSIVE_TRACER
# define NO_FRCFILE

#elif defined KH_INST
/*
!                       Kelvin-Helmholtz Instability Example
!                       ================ =========== =======
!
*/
# undef  KH_INSTY
# undef  KH_INST3D
# undef  MPI
# define NBQ
# undef  NBQ_PRECISE
# undef  XIOS
# define SOLVE3D
# define NEW_S_COORD
# define UV_ADV
# define TS_HADV_WENO5
# define TS_VADV_WENO5
# define UV_HADV_WENO5
# define UV_VADV_WENO5
# define W_HADV_WENO5
# define W_VADV_WENO5
# undef  SALINITY
# undef  PASSIVE_TRACER
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# undef  ANA_SRFLUX
# define ANA_BTFLUX
# define ANA_SSFLUX
# define ANA_BSFLUX
# ifndef KH_INSTY
#  define EW_PERIODIC
# else
#  define NS_PERIODIC
# endif
# define NO_FRCFILE

#elif defined TS_HADV_TEST
/*
!                       Horizontal TRACER ADVECTION EXAMPLE
!                       ========== ====== ========= =======
!
*/
/* Example with spatially varying velocity */
# undef  SOLID_BODY_ROT
/* Constant advection in the diagonal   */
# undef  DIAGONAL_ADV
/* Example with a space and time-varying velocity */
# define SOLID_BODY_PER

# undef  OPENMP
# undef  MPI
# undef  UV_ADV
# define NEW_S_COORD
# undef  UV_COR
# define SOLVE3D
# define M2FILTER_NONE
# define ANA_VMIX
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_SRFLUX
# define ANA_STFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define ANA_SSFLUX
# define NO_FRCFILE
# define SALINITY
# define EW_PERIODIC
# define NS_PERIODIC

/* Choose specific advection scheme */
# define TS_HADV_UP3
# undef  TS_HADV_C4
# undef  TS_HADV_UP5
# undef  TS_HADV_WENO5
# undef  TS_HADV_C6

#elif defined DUNE
/*
!                       Dune test case example
!                       ==== ==== ==== =======
!
*/
/* Analytical test case (Marieu) */
# undef  ANA_DUNE
/* 3D Dune example */
# undef  DUNE3D

# undef  OPENMP
# undef  MPI
# define M2FILTER_NONE
# define UV_ADV
# define NEW_S_COORD
# undef  UV_COR
# define SOLVE3D
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SSFLUX
# define ANA_SRFLUX
# define ANA_STFLUX
# define ANA_BSFLUX
# define ANA_BTFLUX
# define ANA_SMFLUX
# define OBC_WEST
# define OBC_EAST
# define ANA_SSH
# define ZCLIMATOLOGY
# define ANA_M2CLIMA
# define M2CLIMATOLOGY
# define SEDIMENT
# undef  MUSTANG
# define MORPHODYN
# ifdef SEDIMENT
#  undef  SUSPLOAD
#  define BEDLOAD
#  undef  BEDLOAD_WENO5
#  ifdef ANA_DUNE
#   define BEDLOAD_MARIEU
#  else
#   define BEDLOAD_WULIN
#   define TAU_CRIT_WULIN
#  endif
# endif
# ifdef MUSTANG
#  define key_MUSTANG_V2
#  define key_MUSTANG_bedload
#  define key_tauskin_c_upwind
# endif
# define GLS_MIXING
# define NO_FRCFILE

#elif defined SED_TOY
/*
!                       SED TOY (1D Single Column example)
!                       === === === ====== ====== ========
!
*/ 
/* Choose an experiment :               */
/*   Rouse                              */
# define SED_TOY_ROUSE
/*   Consolidation                      */ 
# undef  SED_TOY_CONSOLID
/*   Erosion and sediment resuspension  */
# undef  SED_TOY_RESUSP
/*   Flocculation                       */
# undef  SED_TOY_FLOC_0D
/*   Flocculation                       */
# undef  SED_TOY_FLOC_1D

# undef  OPENMP
# undef  MPI
# define NEW_S_COORD
# define SOLVE3D
# undef  NONLIN_EOS
# define SALINITY
# undef  UV_VIS2
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_SRFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define EW_PERIODIC
# define NS_PERIODIC

# ifdef SED_TOY_ROUSE
#  define ANA_VMIX
#  define BODYFORCE
# endif

# ifdef SED_TOY_FLOC_1D
#  define ANA_VMIX
#  define BODYFORCE
# endif

# ifdef SED_TOY_FLOC_0D
#  define ANA_VMIX
#  define BODYFORCE
# endif

# define SEDIMENT
# undef  MUSTANG

# ifdef MUSTANG
#  if defined SED_TOY_FLOC_0D || defined SED_TOY_FLOC_1D
#   define key_MUSTANG_flocmod
#   define GLS_MIXING
#   define GLS_KOMEGA
#  endif
# endif


# ifdef SEDIMENT
#  define SUSPLOAD
#  undef  BEDLOAD

#  ifdef SED_TOY_ROUSE
#   define SED_TAU_CD_CONST
#  endif

#  if defined SED_TOY_FLOC_1D || defined SED_TOY_CONSOLID || \
	defined SED_TOY_RESUSP
#   undef  BBL
#   define GLS_MIXING
#   define GLS_KOMEGA
#   define MIXED_BED
#   undef  COHESIVE_BED
#  endif

#  if defined SED_TOY_FLOC_0D || defined SED_TOY_FLOC_1D
#   define FLOC_TURB_DISS
#   undef FLOC_BBL_DISS
#   define SED_FLOCS
#   undef SED_DEFLOC
#  endif

# endif
# undef  MORPHODYN
# define NO_FRCFILE

#elif defined KILPATRICK
/*
!                       KILPATRICK  Example
!                       ==========  =======
*/
# define MPI
# define AVERAGES
# define NONLIN_EOS
# define SOLVE3D
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_BTFLUX
# define NO_FRCFILE
# define ABL1D
# ifdef  ABL1D
#  define BULK_FLUX
#  undef  BULK_ECUMEV0
#  undef  BULK_ECUMEV6
#  define BULK_GUSTINESS
#  define ANA_ABL_LSDATA
#  define ANA_ABL_VGRID
#  define STRESS_AT_RHO_POINTS
#  undef  ABL_NUDGING
#  undef  ABL_NUDGING_DYN
#  undef  ABL_NUDGING_TRA
#  undef  ABL_DYN_RESTORE_EQ
#  undef  SFLUX_CFB
# else
#  undef BULK_FLUX
# endif

#elif defined TIDAL_FLAT
/*
!                       TIDAL_FLAT  Example
!                       ==========  =======
*/
# undef  OPENMP
# undef  MPI
# undef  NONLIN_EOS
# define NEW_S_COORD
# define SALINITY
# define UV_ADV
# define TS_HADV_WENO5
# define TS_VADV_WENO5
# define UV_HADV_WENO5
# define UV_VADV_WENO5
# define UV_COR
# define SOLVE3D
# define UV_VIS2
# define GLS_MIXING
# define ANA_INITIAL
# define WET_DRY
# define TS_DIF2
# define SPONGE
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_SRFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define OBC_WEST
# define FRC_BRY
# ifdef FRC_BRY
#  define ANA_BRY
#  define Z_FRC_BRY
#  define OBC_M2CHARACT
#  define OBC_REDUCED_PHYSICS
#  define M2_FRC_BRY
#  undef  M3_FRC_BRY
#  define T_FRC_BRY
# endif
# undef  SEDIMENT
# define MUSTANG
# ifdef SEDIMENT
#  define SUSPLOAD
#  undef  BEDLOAD
# endif
# ifdef MUSTANG
#  define key_sand2D
#  undef  key_MUSTANG_V2
# endif
# define NO_FRCFILE
# undef  ZETA_DRY_IO

#elif defined ESTUARY
/*
!                       ESTUARY  Example
!                       =======  =======
*/
# undef  OPENMP
# undef  MPI
# undef  NONLIN_EOS
# define NEW_S_COORD
# define SALINITY
# define UV_ADV
# define TS_HADV_WENO5
# define TS_VADV_WENO5
# define UV_HADV_WENO5
# define UV_VADV_WENO5
# define UV_COR
# define SOLVE3D
# define UV_VIS2
# define GLS_MIXING
# define ANA_INITIAL
# define WET_DRY
# define TS_DIF2
# define SPONGE
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_SRFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define OBC_WEST
# define FRC_BRY
# ifdef FRC_BRY
#  define ANA_BRY
#  define Z_FRC_BRY
#  define OBC_M2CHARACT
#  define OBC_REDUCED_PHYSICS
#  define M2_FRC_BRY
#  undef  M3_FRC_BRY
#  define T_FRC_BRY
# endif
# undef  SEDIMENT
# define MUSTANG
# ifdef SEDIMENT
#  define SUSPLOAD
#  undef  BEDLOAD
# endif
# ifdef MUSTANG
#  define key_sand2D
#  undef  key_MUSTANG_V2
#  define MUSTANG_CORFLUX
#  define SUBSTANCE_SUBMASSBALANCE 
# endif
# define PSOURCE
# define ANA_PSOURCE
# define MASKING
# define NO_FRCFILE
# undef  ZETA_DRY_IO

#elif defined SEAGRASS
/*
!                       Seagrass example
!                       ================
*/
# define OBSTRUCTION

# undef  OPENMP
# undef  MPI
# define SOLVE3D
# define UV_ADV
# define UV_COR
# define NONLIN_EOS
# define SALINITY
# define ANA_GRID
# define MASKING
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_SRFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define GLS_MIXING
# define PSOURCE
# undef  PSOURCE_MASS
# define ANA_PSOURCE
# define NS_PERIODIC
# undef  FLOATS
# define NO_FRCFILE
# define USE_CALENDAR
# define NEW_S_COORD
# define OBC_EAST
# define FRC_BRY
# ifdef FRC_BRY
#  define ANA_BRY
#  define Z_FRC_BRY
#  define OBC_M2CHARACT
#  define OBC_REDUCED_PHYSICS
#  define M2_FRC_BRY
#  undef  M3_FRC_BRY
#  define T_FRC_BRY
# endif

#endif /* END OF CONFIGURATION CHOICE */

#include "cppdefs_dev.h"
#include "set_global_definitions.h"
