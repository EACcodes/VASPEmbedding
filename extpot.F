#include "symbol.inc"
!****************** module use external potential  *********************
! Humbly written by Kuang... Apr-2014
! This whole code is only tested with PAW
!***********************************************************************

  module mextpot
    use prec
    implicit none

    type extpot
      logical :: lextpot, dextpot
      integer :: nx, ny, nz ! should be equal to the dimension of potential grid
      real(q) :: FCORR
    endtype extpot

    type (extpot) EXTPT
    RGRID, allocatable :: EXTPT_GRID(:,:,:) ! the external potential grid in normal layout
    RGRID, allocatable :: VLM(:,:,:)        ! external potential for each channel of each ion
    RGRID, allocatable :: SCX(:), SCY(:), SCZ(:)
    RGRID, allocatable :: SGRID(:,:,:)      ! the B-spline interpolation coefficients
    ! parameters for real space VLM projection
    integer, parameter :: KMAX=20, di_smear=5, NGRID=50
    real(q), parameter :: Delta = 0.1_q

  contains

    function EXTPT_LEXTPOT()
      logical :: EXTPT_LEXTPOT
      EXTPT_LEXTPOT=EXTPT%lextpot
    end function EXTPT_LEXTPOT

    function EXTPT_DEXTPOT()
      logical :: EXTPT_DEXTPOT
      EXTPT_DEXTPOT=EXTPT%dextpot
    end function EXTPT_DEXTPOT

    ! read INCAR FILE for external potential information
    subroutine EXTPT_READER(NGXF,NGYF,NGZF,IU0,IU5,IU6)
      use base
      use reader_tags
      use vaspxml
      implicit none

      integer :: NGXF, NGYF, NGZF, IU0, IU5, IU6, IU_POT=99
      ! local variables
      integer :: ix, iy, iz, ntot
      logical :: lopen, LDUM
      INTEGER IDUM, IERR, N
      REAL(q) RDUM
      COMPLEX(q) CDUM
      CHARACTER(1) :: CHARAC

      lopen = .false.
      OPEN(UNIT=IU5,FILE=INCAR,STATUS='OLD')
      EXTPT%lextpot = .false.
      EXTPT%dextpot = .false.                 ! initialiting logical dextpot
      EXTPT%FCORR = 0.67_q
      EXTPT%nx = -1
      EXTPT%ny = -1
      EXTPT%nz = -1
      CALL PROCESS_INCAR(INCAR_F, 'LEXTPOT', EXTPT%lextpot)
      CALL PROCESS_INCAR(INCAR_F, 'DEXTPOT', EXTPT%dextpot)
      if( EXTPT%lextpot ) then
         CALL PROCESS_INCAR(INCAR_F, 'FCORR', EXTPT%FCORR)
      endif

      if( EXTPT%lextpot ) then
         open(unit=IU_POT,file='EXTPOT',status='old')
         read(IU_POT,*) EXTPT%nx, EXTPT%ny, EXTPT%nz
         ! check grid size consistency
         if ( EXTPT%nx /= NGXF .or. EXTPT%ny /= NGYF .or. EXTPT%nz /= NGZF ) then
            write(IU6,*) 'ERROR: external potential grid size mismatch:'
            write(IU6,*) '  NGXF, NGYF, NGZF:', NGXF, NGYF, NGZF
            write(IU6,*) '    NX,   NY,   NZ:', EXTPT%nx, EXTPT%ny, EXTPT%nx
            write(IU0,*) 'ERROR: external potential grid size mismatch:'
            write(IU0,*) '  NGXF, NGYF, NGZF:', NGXF, NGYF, NGZF
            write(IU0,*) '    NX,   NY,   NZ:', EXTPT%nx, EXTPT%ny, EXTPT%nx
            stop
         endif
         allocate(EXTPT_GRID(EXTPT%nx, EXTPT%ny, EXTPT%nz))
         ! follow the tradition of CHGCAR format, ix be the fastest index
         read(IU_POT,*) (((EXTPT_GRID(ix,iy,iz), ix=1,EXTPT%nx), &
                                            iy=1,EXTPT%ny), iz=1,EXTPT%nz)
         close(IU_POT)
      endif

      return
    end subroutine EXTPT_READER

!******************** SUBROUTINE EXTERNAL_POT_ADD **********************
! A wrapper deal with the wierd variable type compatibility issue...
! In realmode, CVTOT allocated as complex array, but used as real grid
!***********************************************************************

    SUBROUTINE EXTPT_EXTERNAL_POT_ADD(GRIDC, LATT_CUR, CVTOT)

      USE prec
      USE base
      USE lattice
      USE mpimy
      USE mgrid
      USE poscar
      USE constant

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRIDC
      TYPE (latt)        LATT_CUR
      COMPLEX(q) ::      CVTOT(GRIDC%MPLWV)

      call EXTERNAL_POT_ADD_(GRIDC, LATT_CUR, CVTOT(1))

    end subroutine EXTPT_EXTERNAL_POT_ADD

    subroutine EXTPT_EXTERNAL_POT_ADD_PAW(POTAE, POT, NDIM, NCDIJ, LMMAX, NIP)

      USE prec
      USE base
      USE lattice
      USE mpimy
      USE mgrid
      USE poscar
      USE constant
      implicit none

      INTEGER             NDIM, NCDIJ, LMMAX, NIP
      REAL(q) POT(:,:,:), POTAE(:,:,:)

      integer ISP, LM, II

      do ISP = 1, NCDIJ
         do LM = 1, LMMAX
            do II = 1, NDIM
               POTAE(II,LM,ISP) = POTAE(II,LM,ISP) + VLM(II,LM,NIP) !+1.0_q
               POT(II,LM,ISP) = POT(II,LM,ISP) + VLM(II,LM,NIP) !+1.0_q
            enddo
         enddo
      enddo

    end subroutine EXTPT_EXTERNAL_POT_ADD_PAW



!*********************** SUBROUTINE CALC_VLM ***************************
! This subroutine calculates the LM component of the external potential
! This code was written assuming realmode and NGZhalf = .True.
!***********************************************************************

    SUBROUTINE EXTPT_CALC_VLM(LATT_CUR, GRIDC, T_INFO, P, LMDIM, WDES, IU0)

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
      ! local variables
      type (potcar),POINTER:: PP
      complex(q) :: PHASE
      real(q), allocatable :: CVEXT(:)
      real(q) :: RCUT, R, RI(3),RVEC(3)
      real(q) :: CTMP, RMAX, RSWITCH, val, RPREV, CTMPPREV
      integer, parameter :: nprint=8
      real(q) :: VLMK(KMAX), dV
      integer :: NXMAX, NYMAX, NZMAX, NG, ofile=666
      integer :: NI, NDIM, NT, NIP, LMMAX, LYMAX, L, M, LM, K, RNMAX, NR
      integer :: iprint, IRDMAX, INDMAX, IND
      integer, external :: MAXL_AUG, MAXL1
      RGRID, allocatable :: VLMFULL(:,:,:)
      integer, allocatable :: NLI(:)
      real(q),allocatable :: YLM(:,:), DIST(:), XS(:), YS(:), ZS(:)

      allocate(CVEXT(DIMREAL(GRIDC%MPLWV)))
      NXMAX = EXTPT%nx
      NYMAX = EXTPT%ny
      NZMAX = EXTPT%nz

      ! copy the external potential to CVEXT, change normal layout to VASP layout
      CVEXT = 0.0_q
      call EXTERNAL_POT_ADD_(GRIDC, LATT_CUR, CVEXT(1))

      ! set up the maximum radial grid dimension, copied from SET_DD_PAW
      NDIM=0
      DO NT=1,T_INFO%NTYP
         IF (ASSOCIATED(P(NT)%QPAW)) THEN
            NDIM=MAX(NDIM, P(NT)%R%NMAX)
         END IF
      ENDDO
      IF (NDIM == 0) RETURN
      LYMAX =MAXL_AUG(T_INFO%NTYP,P)
      LMMAX=(LYMAX+1)**2
      ! allocate the VLM array, distributed over COMM_INB
      if (.not. allocated(VLM)) allocate(VLM(NDIM, LMMAX, WDES%NIONS))
      VLM = 0.0_q
      ! allocate the full matrix
      allocate(VLMFULL(NDIM,LMMAX,T_INFO%NIONS))

      ! Initialize B-spline interpolation
      allocate( SGRID(NXMAX,NYMAX,NZMAX) )
      allocate(SCX(-di_smear:di_smear), SCY(-di_smear:di_smear), SCZ(-di_smear:di_smear))
      ! 1D interpolate for delta grid, used later as basis for 3D B-spline
      call CALC_Sdelta_COEFF(NXMAX, SCX, di_smear)
      call CALC_Sdelta_COEFF(NYMAX, SCY, di_smear)
      call CALC_Sdelta_COEFF(NZMAX, sCZ, di_smear)
      call BSPLINE_INTERPOLATE_COEFF( EXTPT_GRID, SGRID, NXMAX, NYMAX, NZMAX, di_smear, &
                                        SCX, SCY, SCZ, GRIDC )

      ! allocate temporary arrays
#ifdef MPI
      IRDMAX = NGRID*NGRID*NGRID/GRIDC%COMM%NCPU+1
#else
      IRDMAX = NGRID*NGRID*NGRID+1
#endif
      allocate( YLM(IRDMAX,LMMAX), DIST(IRDMAX), XS(IRDMAX), &
     &        YS(IRDMAX), ZS(IRDMAX) )
      ! set up YLM for the fine integration grid, this is shared by all atom
      ! centers
      call SETYLM_FINE_GRID( GRIDC, NGRID, YLM, LYMAX, LMMAX, XS, YS, ZS, DIST &
                          , IRDMAX, INDMAX, .True.)

      ! start projection
      VLMFULL = 0.0_q
      ! do ion loops
      ion1: do NI=1,T_INFO%NIONS
         RI(:) = T_INFO%POSION(1:3,NI) ! ion position in direct coordinate
         NT=T_INFO%ITYP(NI)
         PP=> PP_POINTER(P, NI, NT)
         RNMAX =PP%R%NMAX
         RMAX = PP%R%R(RNMAX)
         RCUT = RMAX + Delta ! Want to use RMAX+D as real cutoff
         dV = (RCUT*2/NGRID)**3
         LYMAX =MAXL1(PP)*2
         ! for each VLM
         do L = 0, LYMAX
           do M = 0, L*2
             LM=L*L+M+1
             VLMK(:) = 0.0_q
             do K = 1, KMAX ! radial expansion
               do IND = 1, INDMAX
                 ! \sum YLM *RK *V(r) *dV
                 ! Note this summation is distributed over GRIDC%COMM
                 R = DIST(IND)*RCUT ! in Angstrom
                 if ( R <= RCUT ) then
                   RVEC(1) = XS(IND)*R
                   RVEC(2) = YS(IND)*R
                   RVEC(3) = ZS(IND)*R
                   RVEC = matmul( transpose(LATT_CUR%B(1:3,1:3)), RVEC(1:3) )
                   RVEC(:) = RVEC(:) + RI(:) ! grid coordinate in direct
                   val = bspline( RVEC, SGRID, NXMAX, NYMAX, NZMAX )
                   VLMK(K) = VLMK(K) + dV* FSWITCH(R,RMAX,Delta)* &
                        &  YLM(IND,LM)*RK(R,K,RCUT)*val
                 endif
               enddo
             enddo
             do NR = RNMAX, 1, -1
               R = PP%R%R(NR)
               if ( R>RMAX/KMAX .or. L==0 ) then
                 CTMP = 0.0_q
                 do K = 1, KMAX
                   CTMP = CTMP + RK(R,K,RCUT) * VLMK(K)
                 enddo
                 CTMPPREV = CTMP
                 RPREV = R
               else
                 CTMP = CTMPPREV * (R/RPREV)**L
               endif
               VLMFULL(NR,LM,NI) = CTMP
             enddo
           enddo
         enddo
      enddo ion1
      ! Reduce VLM over GRIDC%COMM
      CALLMPI( M_sum_d(GRIDC%COMM, VLMFULL, NDIM*LMMAX*T_INFO%NIONS) )

      ! distribute the full matrix over COMM_INB to save memory space
      VLM = 0.0_q
      ion2: do NI=1,T_INFO%NIONS
         ! local storage index
         ! if this element is not treated locally CYCLE
         IF (.NOT. DO_LOCAL(NI)) cycle ion2
         NIP=NI_LOCAL(NI, WDES%COMM_INB)
         do L = 0, LYMAX
            do M = 0, L*2
               LM=L*L+M+1
               VLM(:,LM,NIP) = VLMFULL(:,LM,NI)
            enddo
         enddo
      enddo ion2
      CALLMPI( M_sum_d(WDES%COMM_INTER, VLM, NDIM*LMMAX*WDES%NIONS) )


      ! test print
#ifdef MPI
      if (GRIDC%COMM%NODE_ME == 1) then
#endif
         open(unit=ofile, file='VLM.data', action='write')
         do NR = 1, RNMAX
            write(ofile, '(5E15.6)') PP%R%R(NR), VLM(NR,1,1), VLM(NR,2,1),VLM(NR,3,1), VLM(NR,4,1)
         enddo
#ifdef MPI
      endif
#endif
      close(ofile)

#ifdef MPI
      if (GRIDC%COMM%NODE_ME == 1) then
#endif
         write(IU0,'(A)') ' External potential tranformed to VLM.'
#ifdef MPI
      endif
#endif
      deallocate(VLMFULL, CVEXT)

      return
    END SUBROUTINE EXTPT_CALC_VLM

!*************************** subroutine SETYLM_FINE_GRID **************************
! this subroutine setup the fine grid for the real space projection integeration
! I made up this layout to parallel the real space projection code
! Use LPAR to choose whether to distritute the local fine grid points or not
! Note:
! In CALC_VLM, we can not parallelize over uniform grid (GRIDC), so set LPAR = .Ture.
! to parallelize over the local fine grid
! In DE_DVRN, parallelize over GRIDC, so set LPAR = .False., making the local
! fine grid layout serial.
!**********************************************************************************
    subroutine SETYLM_FINE_GRID( GRIDC, NGRID, YLM, LYMAX, LMMAX, &
                                 XS, YS, ZS, DIST, IRDMAX, INDMAX, &
                                 LPAR )
    use prec
    use constant
    use mgrid
    use mpimy
    use asa
    implicit none
      type (grid_3d)    GRIDC
      integer :: NG, LYMAX, LMMAX, IRDMAX, NGRID
      integer :: INDMAX
      real(q) :: XS(IRDMAX), YS(IRDMAX), ZS(IRDMAX)
      real(q) :: YLM(IRDMAX,LMMAX), DIST(IRDMAX)
      integer :: ix, iy, iz, NODE_ME, NCPU, IND
      logical :: LPAR

#ifdef MPI
      NCPU = GRIDC%COMM%NCPU
      NODE_ME = GRIDC%COMM%NODE_ME
#else
      NCPU = 1
      NODE_ME = 1
#endif
      NG = 0
      IND = 0
      do ix = 1, NGRID
        do iy = 1, NGRID
          do iz = 1, NGRID
            NG = NG + 1
            if ( (.not. LPAR) .or. (mod(NG,NCPU)+1 == NODE_ME) ) then
              IND = IND + 1
              XS(IND) = (ix-0.5_q)*2.0_q/NGRID - 1.0_q
              YS(IND) = (iy-0.5_q)*2.0_q/NGRID - 1.0_q
              ZS(IND) = (iz-0.5_q)*2.0_q/NGRID - 1.0_q
              DIST(IND) = sqrt(XS(IND)**2 + YS(IND)**2 + ZS(IND)**2)
              if ( DIST(IND) > 1.0e-8 ) then
                XS(IND) = XS(IND) / DIST(IND)
                YS(IND) = YS(IND) / DIST(IND)
                ZS(IND) = ZS(IND) / DIST(IND)
              endif
            endif
          enddo
        enddo
      enddo
      INDMAX = IND
      call SETYLM(LYMAX,INDMAX,YLM,XS,YS,ZS)

      return
    end subroutine SETYLM_FINE_GRID

!*************************** function bspline **************************
! this function computes the interpolate value, using the interpolation
! coefficeint grid
! r should be the direct coordinate
!***********************************************************************
    function bspline( r, s_grid, nx, ny, nz )
    use prec
    use constant
    implicit none
      integer :: nx, ny, nz
      real(q) :: r(3), s_grid(nx, ny, nz), bspline
      real(q) :: Bx, By, Bz, x, y, z, dx, dy, dz, h
      integer :: dix, diy, diz, ix0, iy0, iz0, ix, iy, iz

      bspline = 0.0_q
      h = 1.0_q
      ! r ~ 0.0-1.0
      ix0 = floor(r(1)*nx)+1
      iy0 = floor(r(2)*ny)+1
      iz0 = floor(r(3)*nz)+1
      do diz = -1, 2
        iz = mod(iz0+diz-1 + 10*nz, nz) + 1
        dz = r(3)*nz - (iz0-1) - diz
        Bz = Bfunc(dz, h)
        do diy = -1, 2
          iy = mod(iy0+diy-1 + 10*ny, ny) + 1
          dy = r(2)*ny - (iy0-1) - diy
          By = Bfunc(dy, h)
          do dix = -1, 2
            ix = mod(ix0+dix-1 + 10*nx, nx) + 1
            dx = r(1)*nx - (ix0-1)  - dix
            Bx = Bfunc(dx, h)
            bspline = bspline + s_grid(ix,iy,iz)*Bx*By*Bz
          enddo
        enddo
      enddo

      return
    end function bspline

!*********************** SUBROUTINE BSPLINE_INTERPOLATE_COEFF ******************
! This subroutine conducts a 3-d B-spline interpolation
! f_grid: the actual input grid
! s_grid: the output spline coefficient grid
! sx, sy, sz : the 1-d coefficients for Delta-grid
!*******************************************************************************
    subroutine BSPLINE_INTERPOLATE_COEFF( f_grid, s_grid, nx, ny, nz, di_smear, &
                                        sx, sy, sz, GRIDC )
    use prec
    use constant
    use mgrid
    use mpimy
    implicit none
      integer :: nx, ny, nz, di_smear, ic, ir
      real(q) :: f_grid(nx,ny,nz), s_grid(nx,ny,nz)
      integer :: ix, iy, iz, ixp, iyp, izp, dix, diy, diz
      real(q) :: sx(-di_smear:di_smear),sy(-di_smear:di_smear),sz(-di_smear:di_smear)
      type (grid_3d)    GRIDC

      s_grid = 0.0_q
      do ic = 1, GRIDC%RL%NCOL
        if (GRIDC%RL%NFAST==3) then ! mpi version: x-> N2, y-> N3, z-> N1
          ix = GRIDC%RL%I2(ic)
          iy = GRIDC%RL%I3(ic)
        else ! conventional version: x-> N1, y-> N2, z-> N3
          iy = GRIDC%RL%I2(ic)
          iz = GRIDC%RL%I3(ic)
        endif
        do ir = 1, GRIDC%RL%NROW
          if (GRIDC%RL%NFAST==3) then
            iz = ir
          else
            ix = ir
          endif
          do diz = -di_smear, di_smear
            izp = mod(iz+diz + 10*nz-1, nz) + 1
            do diy = -di_smear, di_smear
              iyp = mod(iy+diy + 10*ny-1, ny) + 1
              do dix = -di_smear, di_smear
                ixp = mod(ix+dix + 10*nx-1, nx) + 1
                s_grid(ixp,iyp,izp) = s_grid(ixp,iyp,izp) + &
                         sx(dix)*sy(diy)*sz(diz)*f_grid(ix,iy,iz)
              enddo
            enddo
          enddo
        enddo
      enddo
      CALLMPI( M_sum_d(GRIDC%COMM,s_grid,nx*ny*nz))

      return
    end subroutine BSPLINE_INTERPOLATE_COEFF

!*********************** SUBROUTINE CALC_Sdelta_COEFF ******************
! Calculate the 1-D B-spline coefficients for a "delta" grid
!***********************************************************************
    subroutine CALC_Sdelta_COEFF(N, s, di_smear)
    use prec
    use constant
    implicit none
      integer :: N, di_smear
      real(q) :: s(-di_smear:di_smear)
      real(q) :: A(N,N), B(N)
      integer :: i, j, ip, info

      A = 0.0_q
      B = 0.0_q
      do i = 1, N
        do j = -1, 1
          ip = i+j
          ip = mod( (ip-1 + 10*N), N ) + 1
          if ( j==0 ) then
            A(i,ip) = 4.0_q
          else
            A(i,ip) = 1.0_q
          endif
        enddo
      enddo
      B(1) = 1.0_q

      ! use mkl function to find solution
      call DPOTRF( 'U', N, A, N, info ) ! Cholesky decomposition
      call DPOTRS( 'U', N, 1, A, N, B, N, info ) ! solve Ax=B
      do i = 1, di_smear
        s(i) = B(i+1)
        s(-i) = B(n-i+1)
      enddo
      s(0) = B(1)
      s = s*6.0_q

      return
    end subroutine CALC_Sdelta_COEFF

!*********************** SUBROUTINE Bfunc ******************************
! This is the basis functions used in the B-spline interpolation
!***********************************************************************
    function Bfunc(x, h)
    use prec
    use constant
    implicit none
      real(q) :: x, h, Bfunc
      if ( x<-2*h .or. x>2*h) then
        Bfunc = 0.0_q
      else if ( x<-h .or. x>h) then
        Bfunc = 1._q/6*(2*h-abs(x))**3
      else
        Bfunc = 2._q/3*h**3 - 0.5_q*x**2*(2*h-abs(x))
      endif
      return
    end function Bfunc


!*********************** SUBROUTINE READ_ELEM **************************
! Wrappers to read elements from CVTOT
!***********************************************************************
    FUNCTION READ_ELEM(GRIDC, CVTOT, I)
      USE prec
      USE base
      USE lattice
      USE mpimy
      USE mgrid
      USE poscar
      USE constant
      implicit none

      TYPE (grid_3d)     GRIDC
      COMPLEX(q)         CVTOT(DIMREAL(GRIDC%MPLWV))
      INTEGER :: I
      REAL(q) :: READ_ELEM

      CALL READ_ELEM_(GRIDC, CVTOT, I, READ_ELEM)
      return
    END FUNCTION READ_ELEM

!*************************** FUNCTION RK *******************************
!
! This function calculates the radial projection functions
! Right now use the 1/r*sqrt(2/RCUT)*sin(pi*K/RCUT*r)
!
!***********************************************************************
    FUNCTION RK( R, K, L )
      use prec
      use constant
      implicit none

      real(q) :: R, L, RK
      integer :: K

      if ( R < 0.0_q .or. R > L ) then
        write(6,*) 'Internal Error: wrong R value for RK,', R, L
      endif
      if ( R > 1e-5 ) then
        RK = sqrt(2._q/L)*sin(PI*K*R/L) / R
      else
        RK = sqrt(2._q/L)*PI*K/L
      endif

      return

    END FUNCTION RK

!************************ FUNCTION FSWITCH *****************************
!
! This function calculates the switch function f(r) for V(r)
! In principle, can be anything as long as:
! f(RMAX) = 1
! f'(RMAX) = 0
! f(RMAX+Delta) = 0
!
!***********************************************************************
    FUNCTION FSWITCH( R, RMAX, Delta )
    use prec
    use constant
    implicit none
      real*8 :: FSWITCH, R, RMAX, Delta
      real*8 :: dr
      if ( R<RMAX ) then
        FSWITCH = 1.0
      else if ( R>RMAX+Delta ) then
        FSWITCH = 0.0
      else
        dr = (R - RMAX)/Delta
        FSWITCH = cos(dr*PI/2.0)
      endif
      return
    end function FSWITCH

!************************* ENERGY_DERIVATIVE ***************************
! This subroutine augments the electron density to give dE/dV(r)
! Theoretically it should be the AE density
! Pratically it should be softer than AECCAR2 but harder than CHGCAR
! Need:
! CRHODE: PAW occupancy matrix
! CHDEN: soft tilde(n) without compensation charge hat(n)
!        defined in GRID_SOFT
! CHTOT: tilde(n) + hat(n) as input, correct dE/dV as output
!***********************************************************************
    SUBROUTINE EXTPT_DE_DVRN(WDES, GRID_SOFT, GRIDC, SOFT_TO_C, LATT_CUR, &
        P, T_INFO, LMDIM, CRHODE, CHTOT_, CHDEN, IRDMAX )
      use prec
      use base
      use mgrid
      use mpimy
      use lattice
      use poscar
      use pseudo
      use us
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
      COMPLEX(q),TARGET  :: CHTOT_(DIMREAL(GRIDC%MPLWV),WDES%NCDIJ)
      ! local variable
      integer, external :: MAXL_AUG, MAXL1
      integer :: ISP, NT, NI, NIP, LYMAX, LMMAX, INDMAX, IND, NDIM, I
      integer :: LT, LPT, MT, MPT, IGRID
      integer :: L, LP, LL, LLP, M, MP, MPLOW, FAKLLP, LM, LMP, RNMAX, FAKQ
      integer :: LMINDX, ISTART, IEND, IC, NR, NC, NCPU, ICPU
      integer :: NXMAX, NYMAX, NZMAX, NALLOC, NGF, JLL, JLM, ierror
      integer, allocatable :: NLI(:), LLOCAL(:)
      real(q), allocatable :: CWORK(:)
      COMPLEX(q),POINTER :: CHTOT(:)
      integer, parameter :: nprint=8
      integer :: ofile
      integer :: INDFMAX, INDF, K, ix, iy, iz, ixp, iyp, izp, dix, diy, diz
      real(q) :: VLMK(KMAX), dV, RMAX, RCUT
      real(q) :: RI(3), RVEC(3), RVECF(3), DRVEC(3), val
      real(q) :: STMP1, STMP2, STMP3, STMP4, R, STMP1PREV, RPREV
      real(q), allocatable :: DIST(:), QIJ(:)
      real(q), allocatable :: XS(:),YS(:), ZS(:), SUM(:)
      real(q), allocatable :: YLMF(:,:), XSF(:), YSF(:), ZSF(:), DISTF(:)
      OVERLAP, allocatable :: CTMP(:,:,:,:)
      character(100) :: fn
      character(50) :: buffer
      logical :: LINFO

      ! temparary storage matrix for rho_ij
      ALLOCATE(CTMP(LMDIM,LMDIM,T_INFO%NIONS,WDES%NCDIJ))
      allocate(CWORK(DIMREAL(GRIDC%MPLWV)))

      LYMAX=MAXL_AUG(T_INFO%NTYP,P)
      LMMAX=(LYMAX+1)**2          ! number of lm pairs
#ifdef MPI
      NCPU = GRIDC%COMM%NCPU
#else
      NCPU = 1
#endif

      ! Construct the full rho_ij  matrix
      CTMP=0
      DO ISP=1,WDES%NCDIJ
         DO NI=1,T_INFO%NIONS
            NIP=NI_LOCAL(NI, WDES%COMM_INB)
            IF (NIP/=0) THEN
               CTMP(:,:,NI,ISP)=CRHODE(:,:,NIP,ISP)
            ENDIF
         ENDDO
      ENDDO
      ! Reduce rho_ij for all nodes, assume realmode
      CALLMPI( M_sum_d(WDES%COMM_INB,CTMP,LMDIM*LMDIM*T_INFO%NIONS*WDES%NCDIJ))

      ALLOCATE( XS(IRDMAX), YS(IRDMAX), ZS(IRDMAX), &
               DIST(IRDMAX),SUM(IRDMAX), &
               NLI(IRDMAX), LLOCAL(IRDMAX))
      NDIM=0
      DO NT=1,T_INFO%NTYP
         IF (ASSOCIATED(P(NT)%QPAW)) THEN
            NDIM=MAX(NDIM, P(NT)%R%NMAX)
         END IF
      ENDDO
      allocate(QIJ(NDIM))

#ifdef MPI
      ICPU = WDES%COMM%NODE_ME
#else
      ICPU = 1
#endif
      NXMAX = EXTPT%nx
      NYMAX = EXTPT%ny
      NZMAX = EXTPT%nz
      NALLOC = NGRID*NGRID*NGRID
      allocate( YLMF(NALLOC,LMMAX), DISTF(NALLOC), XSF(NALLOC), &
     &        YSF(NALLOC), ZSF(NALLOC) )
      CALLMPI( MPI_barrier( GRIDC%COMM, ierror ))
      inquire(DIRECTORY="DERINFO", EXIST=LINFO) ! does external information exists?
      CALLMPI( MPI_barrier( GRIDC%COMM, ierror ))
      if ( .not. LINFO ) call system('mkdir DERINFO')
      CALLMPI( MPI_barrier( GRIDC%COMM, ierror ))
      write( buffer, "(I)" ) ICPU
      buffer = trim(adjustl(buffer))
      write( fn, "('DERINFO/',A)" ) buffer
      ofile = 666+ICPU-1
      ! if DERINFO exists, open files to read
      if ( LINFO ) then
        if ( ICPU == 1 ) print *, 'read DERINFO'
        open( unit=ofile, file=fn, action='read' )
      ! else, open files to write and prepare for the real calculations
      else
        if ( ICPU == 1 ) print *, 'write DERINFO'
        open( unit=ofile, file=fn, action='write' )
        ! set up YLM for the fine integration grid, this is shared by all atom
        ! centers
        call SETYLM_FINE_GRID( GRIDC, NGRID, YLMF, LYMAX, LMMAX, XSF &
                              , YSF, ZSF, DISTF, NALLOC, INDFMAX, .false.)
      endif

      ! Do the actual work, only work with the total (alpha+beta) component
      spin: do ISP = 1, 1 !WDES%NCDIJ
      CHTOT => CHTOT_(:,ISP)
      ion: do NI = 1, T_INFO%NIONS
         if ( ICPU == 1 ) write(*,*) 'Ion:', NI
         CWORK = 0.0_q
         NT = T_INFO%ITYP(NI)
         RNMAX = P(NT)%R%NMAX
         RMAX = P(NT)%R%R(RNMAX)
         RCUT = RMAX + Delta
         dV = (RCUT*2/NGRID)**3
         RI = T_INFO%POSION(1:3,NI)
         ! Find the relevant points on GRIDC, in my equation, there is no
         ! YLM(hat(r)), hence no need to use SETYLM_AUG
         ! In this version, we parallelize over local fine grid, so for each node,
         ! should have a complete grid point list, rather than local list
         CALL SETRN(GRIDC,LATT_CUR,RI,P(NT)%PSDMAX*EXTPT%FCORR,NPSRNL,&
               IRDMAX,INDMAX, 0.0,0.0,0.0,DIST(1),NLI(1),XS,YS,ZS,LLOCAL)
         rn_loop: do IND = 1, INDMAX
            SUM(IND) = 0.0_q
#ifdef MPI
            if ( mod(IND,GRIDC%COMM%NCPU)+1 /= ICPU ) goto 321 ! if not local grid point, skip
#endif
            SGRID = 0.0_q
            if ( .not. LINFO ) then
              RVEC(1) = XS(IND)*DIST(IND) ! in angstrom
              RVEC(2) = YS(IND)*DIST(IND)
              RVEC(3) = ZS(IND)*DIST(IND)
              RVEC = matmul( transpose(LATT_CUR%B(1:3,1:3)), RVEC(1:3)) ! in direct
              ix = nint((RVEC(1)+RI(1))*NXMAX)
              iy = nint((RVEC(2)+RI(2))*NYMAX)
              iz = nint((RVEC(3)+RI(3))*NZMAX)
              ix = mod(ix + 10*NXMAX, NXMAX) + 1
              iy = mod(iy + 10*NYMAX, NYMAX) + 1
              iz = mod(iz + 10*NZMAX, NZMAX) + 1
              ! construct SGRID for a delta grid centered at (XS,YS,ZS)
              do diz = -di_smear, di_smear
                izp = mod(iz+diz + 10*NZMAX-1, NZMAX) + 1
                do diy = -di_smear, di_smear
                  iyp = mod(iy+diy + 10*NYMAX-1, NYMAX) + 1
                  do dix = -di_smear, di_smear
                    ixp = mod(ix+dix + 10*NXMAX-1, NXMAX) + 1
                    SGRID(ixp,iyp,izp) = SGRID(ixp,iyp,izp) + &
                                    SCX(dix)*SCY(diy)*SCZ(diz)
                  enddo
                enddo
              enddo
            endif

            ! four DO loops over all the channels: l, m, lp, mp
            LM = 1
            l_loop:  DO L =1,P(NT)%LMAX
            LMP=LM
            lp_loop: DO LP=L,P(NT)%LMAX

            if ( .not. LINFO ) then
              ! setup QIJ
              call SETQIJ(GRIDC, L,LP,P(NT),QIJ, NDIM)
            endif
            ! Electron density scaled by occupation number
            LL =P(NT)%LPS(L )
            LLP=P(NT)%LPS(LP)
            if ( .not. LINFO ) then
              ! setup the YLM3 lookup index
              CALL YLM3LOOKUP(LL,LLP,LMINDX)
            endif

            m_loop:  DO M=1,2*LL+1
            MPLOW=1
            IF (L==LP) MPLOW=M
            mp_loop: DO MP=1,2*LLP+1 !MPLOW,2*LLP+1
               FAKLLP=1
               !IF (LMP+MP/=LM+M) FAKLLP=2
               IF (L/=LP) FAKLLP=2
               if ( .not. LINFO ) then
                 ! setup YLM3 for each (l, m, lp, mp)
                 LMINDX=LMINDX+1
                 ISTART=INDCG(LMINDX)
                 IEND  =INDCG(LMINDX+1)
               endif

               ! Loop over L, M
               ! JL(IC) would point to L
               ! JS(IC) would point to L*(L+1) + M + 1, for YLM3(IC)
               STMP3 = 0.0_q

               if (LINFO) then
               !*****************************
               ! read STMP3 from DERINFO
               !*****************************
               read(ofile,'(I4,I8,4I4,E25.10)') NIP,IGRID,LT,MT,LPT,MPT,STMP3
               ! check consistency
               if (NIP/=NI .or. IGRID/=ABS(LLOCAL(IND)) .or. LT/=L .or. MT/=M &
                  .or. LPT/=LP .or. MPT/=MP ) then
                  write(6,*) 'Internal error, mismatch with DERINFO'
               endif

               else
               !*****************************
               ! compute and write DERINFO
               !*****************************
               DO  IC=ISTART,IEND-1
                  JLL = JL(IC)
                  JLM = JS(IC)

                  ! find out VLMK for Delta(XS, YS, ZS)
                  VLMK(:) = 0.0_q
                  do K = 1, KMAX
                    ! \sum YLM *RK *V(r) *dV, in fine grid
                    do INDF = 1, INDFMAX
                      R = DISTF(INDF)*RCUT ! in Angstrom
                      if ( R <= RCUT ) then
                        RVECF(1) = XSF(INDF)*R
                        RVECF(2) = YSF(INDF)*R
                        RVECF(3) = ZSF(INDF)*R
                        RVECF = matmul( transpose(LATT_CUR%B(1:3,1:3)), RVECF(1:3) )
                        DRVEC = RVECF - RVEC ! in direct
                        if(abs(DRVEC(1)*NXMAX)<di_smear+1 .and. &
                           abs(DRVEC(2)*NYMAX)<di_smear+1 .and. &
                           abs(DRVEC(3)*NZMAX)<di_smear+1) then
                          RVECF(:) = RVECF(:) + RI(:) ! grid coordinate in direct
                          val = bspline( RVECF, SGRID, NXMAX, NYMAX, NZMAX )
                          VLMK(K) = VLMK(K) + dV* FSWITCH(R,RMAX,Delta)* &
                            &  YLMF(INDF,JLM)*RK(R,K,RCUT)*val
                        endif
                      endif
                    enddo
                  enddo
                  ! loop over r
                  STMP2 = 0.0_q
                  do NR = RNMAX, 1, -1
                    R = P(NT)%R%R(NR)
                    if ( R>RMAX/KMAX .or. JLL==0 ) then
                      STMP1 = 0.0_q
                      do K = 1, KMAX
                        STMP1 = STMP1 + RK(R,K,RCUT) * VLMK(K)
                      enddo
                      STMP1PREV = STMP1
                      RPREV = R
                    else ! in the r->0 region
                      STMP1 = STMP1PREV * ((R/RPREV)**JLL)
                    endif
                    ! \sum_r W(r)*QIJ(r)
                    STMP2 = STMP2 + STMP1*P(NT)%R%SI(NR)* &
                           (QIJ(NR) - P(NT)%QPAW(L,LP,JLL)*P(NT)%AUG(NR,JLL))
                  enddo
                  ! \sum_LM C(LM, l, m, lp, mp)
                  STMP3 = STMP3 + STMP2*YLM3(IC)
               ENDDO
               ! write the external info
               write(ofile,'(I4,I8,4I4,E25.10)') NI,abs(LLOCAL(IND)),L,M,LP,MP,STMP3
               endif
               ! \sum_llpmmp rho_ij
               SUM(IND) = SUM(IND) + STMP3*CTMP(LM+M-1,LMP+MP-1,NI,ISP)*FAKLLP

            ENDDO mp_loop
            ENDDO m_loop
            LMP=LMP+2*LLP+1
            ENDDO lp_loop
            LM =LM +2*LL +1
            ENDDO l_loop

 321        continue
         enddo rn_loop

         ! synchronize the rn_loop
         CALLMPI( M_sum_d(GRIDC%COMM, SUM, INDMAX) )

         do IND = 1, INDMAX
            ! when adding, only deal with local grid points
            if (LLOCAL(IND) > 0) then
              CWORK(NLI(IND)) = CWORK(NLI(IND)) + SUM(IND)*NXMAX*NYMAX*NZMAX
            endif
         enddo
         ! FFT to reciprocal space
         CALL FFT_RC_SCALE(CWORK(1),CWORK(1),GRIDC)
         CHTOT = CWORK + CHTOT
      enddo ion
      enddo spin

      deallocate(XS,YS,ZS,DIST,SUM,NLI,LLOCAL,QIJ,YLMF,DISTF,XSF,YSF,ZSF,CTMP,CWORK)
      close( ofile )

      return
    END SUBROUTINE EXTPT_DE_DVRN

end module mextpot

!************************ SUBROUTINE SETRN  ****************************
!
! this subroutine performes the following tasks
! ) finds the points, which are within a certain cutoff around one ion
! ) calculates the distance of each of this points from the ion
! DISX,Y,Z are additional displacements of the ions
!
! mind that the cutoff-sphere extends up to PSDMAX*(NPSRNL-1)/NPSRNL
! This code is a direct copy from SETYLM_AUG, with the YLM part removed
!
!***********************************************************************

    SUBROUTINE SETRN(GRID,LATT_CUR,POSION,PSDMAX,NPSRNL, &
             IRMAX,INDMAX,DISX,DISY,DISZ,DIST,NLI, &
             XS,YS,ZS, LLOCAL)
      USE prec

      USE mpimy
      USE mgrid
      USE lattice
      USE asa
      USE constant
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)  GRID
      TYPE (latt)     LATT_CUR

      DIMENSION POSION(3)

      DIMENSION DIST(IRMAX)
! work-arrays
      DIMENSION NLI(IRMAX)
      REAL(q) XS(IRMAX),YS(IRMAX),ZS(IRMAX)
      INTEGER :: LLOCAL(IRMAX)
      XS=0;YS=0;ZS=0; LLOCAL=0
!=======================================================================
! find lattice points contained within the cutoff-sphere
! mind that the cutoff-sphere extends up to PSDMAX*(NPSRNL-1)/NPSRNL
! which is a somewhat strange convention
!=======================================================================
      F1=1._q/GRID%NGX
      F2=1._q/GRID%NGY
      F3=1._q/GRID%NGZ

      ARGSC=NPSRNL/PSDMAX
!-----------------------------------------------------------------------
! restrict loop to points contained within a cubus around the ion
!-----------------------------------------------------------------------
      D1= PSDMAX*LATT_CUR%BNORM(1)*GRID%NGX
      D2= PSDMAX*LATT_CUR%BNORM(2)*GRID%NGY
      D3= PSDMAX*LATT_CUR%BNORM(3)*GRID%NGZ

      N3LOW= INT(POSION(3)*GRID%NGZ-D3+10*GRID%NGZ+.99_q)-10*GRID%NGZ
      N2LOW= INT(POSION(2)*GRID%NGY-D2+10*GRID%NGY+.99_q)-10*GRID%NGY
      N1LOW= INT(POSION(1)*GRID%NGX-D1+10*GRID%NGX+.99_q)-10*GRID%NGX

      N3HI = INT(POSION(3)*GRID%NGZ+D3+10*GRID%NGZ)-10*GRID%NGZ
      N2HI = INT(POSION(2)*GRID%NGY+D2+10*GRID%NGY)-10*GRID%NGY
      N1HI = INT(POSION(1)*GRID%NGX+D1+10*GRID%NGX)-10*GRID%NGX

!-----------------------------------------------------------------------
! MPI version z ist the fast index
!-----------------------------------------------------------------------
      IND=1

      DO N2=N2LOW,N2HI
      X2=(N2*F2-POSION(2))
      N2P=MOD(N2+10*GRID%NGY,GRID%NGY)

      DO N1=N1LOW,N1HI
      X1=(N1*F1-POSION(1))
      N1P=MOD(N1+10*GRID%NGX,GRID%NGX)

      NCOL=GRID%RL%INDEX(N1P,N2P)
      ! Nolonger cycle, but set a tag(LLOCAL) later
      ! IF (NCOL==0) CYCLE
      if (NCOL>0) then
        IF (GRID%RL%I2(NCOL) /= N1P+1 .OR. GRID%RL%I3(NCOL) /= N2P+1) THEN
          WRITE(*,*)'internal ERROR SETYLM_AUG:',N1P+1,N2P+1,NCOL, &
             GRID%RL%I2(NCOL),N1P , GRID%RL%I3(NCOL),N2P
          STOP
        ENDIF
     endif

!OCL SCALAR
      DO N3=N3LOW,N3HI
        X3=(N3*F3-POSION(3))

        X= X1*LATT_CUR%A(1,1)+X2*LATT_CUR%A(1,2)+X3*LATT_CUR%A(1,3)
        Y= X1*LATT_CUR%A(2,1)+X2*LATT_CUR%A(2,2)+X3*LATT_CUR%A(2,3)
        Z= X1*LATT_CUR%A(3,1)+X2*LATT_CUR%A(3,2)+X3*LATT_CUR%A(3,3)

        D=SQRT(X*X+Y*Y+Z*Z)
        ARG=(D*ARGSC)+1
        NADDR=INT(ARG)

        IF (NADDR<NPSRNL) THEN
          N3P=MOD(N3+10*GRID%NGZ,GRID%NGZ)
          if (NCOL == 0) then
            LLOCAL(IND) = -(1+N3P+N2P*GRID%NGZ+N1P*GRID%NGZ*GRID%NGY)
          else
            LLOCAL(IND) = 1+N3P+N2P*GRID%NGZ+N1P*GRID%NGZ*GRID%NGY
          endif
          NLI (IND)=1+N3P+ GRID%NGZ*(NCOL-1)

          ZZ=Z-DISZ
          YY=Y-DISY
          XX=X-DISX
! the calculation of the | R(ion)-R(mesh)+d | for displaced ions
! was done using  | R+d | = | R | + d . R/|R|
!          IF (D<1E-4_q) THEN
!            DIST(IND)=1E-4_q
!          ELSE
!            DIST(IND)=MAX(D-(DISX*X+DISY*Y+DISZ*Z)/D,1E-10_q)
!          ENDIF
! new version correct to second order
          DIST(IND)=MAX(SQRT(XX*XX+YY*YY+ZZ*ZZ),1E-15_q)

          XS(IND)  =XX/DIST(IND)
          YS(IND)  =YY/DIST(IND)
          ZS(IND)  =ZZ/DIST(IND)

          IND=IND+1
        ENDIF
      ENDDO
      ENDDO; ENDDO
!-----------------------------------------------------------------------
!  compare maximum index with INDMAX
!-----------------------------------------------------------------------
      INDMAX=IND-1
      IF (INDMAX>IRMAX) THEN
        WRITE(*,*) &
     &  'internal ERROR: DEPLE:  IRDMAX must be increased to',INT(INDMAX*1.1)
        STOP
      ENDIF

      RETURN
    END SUBROUTINE

!************************ SUBROUTINE SETQIJ ****************************
!
! this subroutine setup QIJ(r) for a particular type of ion center
! Q_ij(r) = (\phi(r)\phi(r) - \phi_soft(r)\phi_soft(r))
! Note compare to GETQIJ, we do not divide charge by r^2
!
!***********************************************************************
    SUBROUTINE SETQIJ(GRIDC, L,LP,P,QIJ,NDIM)

      use mgrid
      USE prec
      USE pseudo

      IMPLICIT NONE

      TYPE (potcar) P

      TYPE (grid_3d)     GRIDC
      INTEGER L,LP,NDIM
      REAL(q) :: QIJ(NDIM)

      INTEGER I,N,LO,HI
      REAL(q) X(P%R%NMAX),Y(P%R%NMAX)
      REAL(q) DY

      DO I=1,P%R%NMAX
         QIJ(I)=(P%WAE(I,L)*P%WAE(I,LP)-P%WPS(I,L)*P%WPS(I,LP)) !/P%R%R(I)/P%R%R(I)
      ENDDO

      RETURN
    END SUBROUTINE

!******************** SUBROUTINE EXTERNAL_POT_ADD_ *********************
!
! this subroutine can be used to add an external potential
! the units of the potential are eV
! bifurcated from soubroutine EXTERNAL_POT in pot.F
!
!***********************************************************************
    SUBROUTINE EXTERNAL_POT_ADD_(GRIDC, LATT_CUR, CVTOT)
      USE prec
      USE base
      USE lattice
      USE mpimy
      USE mgrid
      USE poscar
      USE constant
      USE mextpot
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRIDC
      TYPE (latt)        LATT_CUR
      RGRID      CVTOT(DIMREAL(GRIDC%MPLWV))
      INTEGER :: NG, N1, N2, N3, N1MAX, N2MAX, N3MAX, NC

      NG=0

! by kuang
      IF (GRIDC%RL%NFAST==3) THEN
         ! mpi version: x-> N2, y-> N3, z-> N1
         N2MAX=GRIDC%NGX
         N3MAX=GRIDC%NGY
         N1MAX=GRIDC%NGZ

         DO NC=1,GRIDC%RL%NCOL
            N2= GRIDC%RL%I2(NC)
            N3= GRIDC%RL%I3(NC)
            DO N1=1,GRIDC%RL%NROW
               NG=NG+1
               CVTOT(NG)=CVTOT(NG) + EXTPT_GRID(N2,N3,N1)
            ENDDO
         ENDDO
      ELSE
         ! conventional version: x-> N1, y-> N2, z-> N3
         N1MAX=GRIDC%NGX
         N2MAX=GRIDC%NGY
         N3MAX=GRIDC%NGZ

         DO NC=1,GRIDC%RL%NCOL
            N2= GRIDC%RL%I2(NC)
            N3= GRIDC%RL%I3(NC)
            DO N1=1,GRIDC%RL%NROW
               NG=NG+1
               CVTOT(NG)=CVTOT(NG) + EXTPT_GRID(N1,N2,N3)
            ENDDO
         ENDDO
      ENDIF

    END SUBROUTINE EXTERNAL_POT_ADD_

    SUBROUTINE READ_ELEM_(GRIDC, CVTOT, I, RES)
      USE prec
      USE base
      USE lattice
      USE mpimy
      USE mgrid
      USE poscar
      USE constant
      implicit none

      TYPE (grid_3d)     GRIDC
      RGRID              CVTOT(DIMREAL(GRIDC%MPLWV))
      INTEGER :: I
      REAL(q) :: RES

      RES = CVTOT(I)
    END SUBROUTINE READ_ELEM_
