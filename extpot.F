#include "symbol.inc"
!***********************************************************************
! These are stubs for the Yu/Libisch/Carter PAW DFET embedding routines
!***********************************************************************

   subroutine EXTPT_EXTERNAL_POT_ADD(GRIDC, LATT_CUR, CVTOT)
      USE mextpot
      USE prec
      USE base
      USE lattice
      USE mpimy
      USE mgrid
      USE poscar
      USE constant
      implicit none

      TYPE (grid_3d)     GRIDC
      TYPE (latt)        LATT_CUR
      COMPLEX(q) ::      CVTOT(GRIDC%MPLWV)
      call EXTERNAL_POT_ADD(GRIDC, LATT_CUR, CVTOT)
   end subroutine EXTPT_EXTERNAL_POT_ADD

   subroutine EXTPT_EXTERNAL_POT_ADD_PAW(POTAE, POT, NDIM, NCDIJ, LMMAX, NIP)
      USE mextpot
      USE prec
      USE base
      USE lattice
      USE mpimy
      USE mgrid
      USE poscar
      USE constant
      implicit none

      INTEGER             NDIM, NCDIJ, LMMAX, NIP
      REAL(q) POT(NDIM,LMMAX,NCDIJ), POTAE(NDIM,LMMAX,NCDIJ)
      call ADD_EXTERNAL_POTENTIAL(POTAE, POT, NDIM, NCDIJ, LMMAX, NIP)
   end subroutine EXTPT_EXTERNAL_POT_ADD_PAW

   subroutine EXTPT_CALC_VLM(LATT_CUR, GRIDC, T_INFO, P, LMDIM, WDES, IU0)
      use mextpot
      use prec
      use base
      use lattice
      use mpimy
      use constant
      use mgrid
      use poscar
      use pseudo
      use paw
      use wave
      use us
      use asa
      implicit none

      type (grid_3d)    GRIDC
      type (latt)       LATT_CUR
      type (type_info)  T_INFO
      type (potcar),TARGET::  P(T_INFO%NTYP)
      type (wavedes)    WDES
      integer :: LMDIM, IU0
      call CALC_VLM(LATT_CUR, GRIDC, T_INFO, P, LMDIM, WDES, IU0)
   end subroutine EXTPT_CALC_VLM

   subroutine EXTPT_DE_DVRN(WDES, GRID_SOFT, GRIDC, SOFT_TO_C, LATT_CUR, &
        P, T_INFO, LMDIM, CRHODE, CHTOT_, CHDEN, IRDMAX )
      use mextpot
      use prec
      use base
      use mgrid
      use mpimy
      use lattice
      use poscar
      use pseudo
      use paw
      use wave
      use asa
      use constant
      implicit none

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d)     GRIDC
      TYPE (grid_3d)     GRID_SOFT
      TYPE (latt)        LATT_CUR
      TYPE (transit)     SOFT_TO_C
      TYPE (wavedes)     WDES
      integer :: LMDIM, IRDMAX
      COMPLEX(q) CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)
      OVERLAP   CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      COMPLEX(q), TARGET :: CHTOT_(DIMREAL(GRIDC%MPLWV),WDES%NCDIJ)
      call DE_DVRN(WDES, GRID_SOFT, GRIDC, SOFT_TO_C, LATT_CUR, &
        P, T_INFO, LMDIM, CRHODE, CHTOT_, CHDEN, IRDMAX)
   end subroutine EXTPT_DE_DVRN
