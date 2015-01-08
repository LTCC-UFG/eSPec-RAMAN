! interfaces to LAPACK procedures
! Pawel Salek, pawsa@ifm.liu.se
! NOTE:
! 980422: also some BLAS routines added at the end.
MODULE NumberTypes
  PUBLIC
  SAVE

#ifdef SYS_CRAY
  INTEGER, PARAMETER :: SP = KIND(1.0E0)
  INTEGER, PARAMETER :: mk = SP
#else
  INTEGER, PARAMETER :: SP = KIND(1.0E0)
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  INTEGER, PARAMETER :: mk = DP
#endif

  REAL(mk),PARAMETER :: PI = 3.14159265358979312
end MODULE NumberTypes

MODULE LapackInterface
  USE NumberTypes
  ! gesv - Solution to a Linear System in a General Matrix (Simple Driver)
  INTERFACE gesv
     SUBROUTINE sgesv(N, NRHS, SA, LDA, IPIVOT, SB, LDB, INFO)
       USE NumberTypes; IMPLICIT NONE
       INTEGER, INTENT(IN)     :: N, NRHS, LDA, LDB
       REAL(SP), INTENT(INOUT) :: SA, SB ! SA(LDA,N), SB(LDB,NRHS)
       INTEGER, INTENT(OUT)    :: IPIVOT !(N)
       INTEGER, INTENT(OUT)    :: INFO
     END SUBROUTINE sgesv
#ifndef SYS_CRAY
     SUBROUTINE dgesv(N, NRHS, DA, LDA, IPIVOT, DB, LDB, INFO)
       USE NumberTypes; IMPLICIT NONE
       INTEGER, INTENT(IN)     :: N, NRHS, LDA, LDB
       REAL(DP), INTENT(INOUT) :: DA, DB ! DA(LDA,N), DB(LDB,NRHS)
       INTEGER, INTENT(OUT)    :: IPIVOT !(N)
       INTEGER, INTENT(OUT)    :: INFO
     END SUBROUTINE dgesv
#endif
  END INTERFACE
  
  INTERFACE stebz
     SUBROUTINE sstebz(RANGE, ORDER, N, Vl, Vu, IL, IU, ABSTOL,&
          & D, E, M, NSPLIT, W, IBLOCK, ISPLIT, WORK, &
          & IWORK, INFO)
       USE NumberTypes; IMPLICIT NONE
       CHARACTER, INTENT(IN)      :: RANGE, ORDER
       INTEGER, INTENT(IN)        :: N, IL, IU
       REAL(SP), INTENT(IN)       :: VL, VU, ABSTOL
       REAL(SP)                   :: D, E, W, WORK
       INTEGER, INTENT(OUT)       :: M, NSPLIT, INFO
       INTEGER, INTENT(OUT)       :: IBLOCK, ISPLIT, IWORK
     end SUBROUTINE sstebz
#ifndef SYS_CRAY
     SUBROUTINE dstebz(RANGE, ORDER, N, Vl, Vu, IL, IU, ABSTOL, &
       D, E, M, NSPLIT, W, IBLOCK, ISPLIT, WORK, &
       IWORK, INFO)
       USE NumberTypes; IMPLICIT NONE
       CHARACTER, INTENT(IN)      :: RANGE, ORDER
       INTEGER,INTENT(IN)         :: N, IL, IU
       REAL(DP), INTENT(IN)       :: VL, VU, ABSTOL
       REAL(DP)                   :: D, E, W, WORK
       INTEGER, INTENT(OUT)       :: M, NSPLIT, INFO
       INTEGER, INTENT(OUT)       :: IBLOCK, ISPLIT, IWORK
     end SUBROUTINE dstebz
#endif
  end INTERFACE


  INTERFACE stein
       SUBROUTINE sstein(N, D, E, M, W, IBLOCK, ISPLIT, Z, LDZ, WORK, &
          IWORK, IFAIL, INFO)
       USE NumberTypes; IMPLICIT NONE
       INTEGER, INTENT(IN)                 :: N, M, LDZ
       REAL(SP), INTENT(IN)  :: D, E, W
       INTEGER, INTENT(IN)   :: iblock, isplit
       REAL(SP), INTENT(OUT) :: Z
       REAL(SP)              :: WORK
       INTEGER               :: IWORK
       INTEGER, INTENT(OUT)                :: ifail, info
     end SUBROUTINE sstein
#ifndef SYS_CRAY
     SUBROUTINE dstein(N, D, E, M, W, IBLOCK, ISPLIT, Z, LDZ, WORK, &
          IWORK, IFAIL, INFO)
       USE NumberTypes; IMPLICIT NONE
       INTEGER, INTENT(IN)                 :: N, M, LDZ
       REAL(DP),  INTENT(IN)  :: D, E, W
       INTEGER,  INTENT(IN)   :: iblock, isplit
       REAL(DP), INTENT(OUT)  :: Z
       REAL(DP)               :: WORK
       INTEGER                :: IWORK
       INTEGER, INTENT(OUT)   :: ifail, info
     end SUBROUTINE dstein
#endif
  end INTERFACE

  ! packed array eigenvalue problem, expert driver ---------------------------
  INTERFACE spevx
       SUBROUTINE sspevx(JOBZ, RANGE, UPLO, N, SA, SVL, SVU, IL, IU, SABTOL, &
            & NFOUND, SW, SZ, LDZ, SWORK, IWORK2, IFAIL, INFO)
       USE NumberTypes; IMPLICIT NONE
       CHARACTER, INTENT(IN)   :: JOBZ, RANGE, UPLO
       INTEGER, INTENT(IN)     :: N, IL, IU, LDZ
       REAL(SP), INTENT(INOUT) :: SA ! array
       REAL(SP), INTENT(IN)    :: SVL, SVU, SABTOL
       INTEGER, INTENT(OUT)    :: NFOUND, IWORK2, IFAIL, INFO
       REAL(SP), INTENT(OUT)   :: SW, SZ, SWORK
     END SUBROUTINE sspevx
     
#ifndef SYS_CRAY
       SUBROUTINE dspevx(JOBZ, RANGE, UPLO, N, DA, DVL, DVU, IL, IU, DABTOL, &
            & NFOUND, DW, DZ, LDZ, DWORK, IWORK2, IFAIL, INFO)
       USE NumberTypes; IMPLICIT NONE
       CHARACTER, INTENT(IN)   :: JOBZ, RANGE, UPLO
       INTEGER, INTENT(IN)     :: N, IL, IU, LDZ
       REAL(DP), INTENT(INOUT) :: DA ! array
       REAL(DP), INTENT(IN)    :: DVL, DVU, DABTOL
       INTEGER, INTENT(OUT)    :: NFOUND, IWORK2, IFAIL, INFO
       REAL(DP), INTENT(OUT)   :: DW, DZ, DWORK
     END SUBROUTINE dspevx
#endif
  end INTERFACE

  INTERFACE steqr
     SUBROUTINE SSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
       USE NumberTypes; IMPLICIT NONE
       CHARACTER, INTENT(IN)  :: compz
       INTEGER, INTENT(IN)    :: N, ldz
       REAL(SP)               :: D, E, work
       REAL(SP)               :: Z
       INTEGER, INTENT(OUT)   :: info
     end SUBROUTINE SSTEQR
#ifndef SYS_CRAY
      SUBROUTINE DSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
       USE NumberTypes; IMPLICIT NONE
       CHARACTER, INTENT(IN)  :: compz
       INTEGER, INTENT(IN)    :: N, ldz
       REAL(DP)               :: D, E, work
       REAL(DP)               :: Z
       INTEGER, INTENT(OUT)   :: info
     end SUBROUTINE DSTEQR
#endif
  end INTERFACE

  INTERFACE syev
       SUBROUTINE ssyev(JOBZ, UPLO, N, SA, LDA, SW, SWORK, LWORK, INFO)
       USE NumberTypes; IMPLICIT NONE
       CHARACTER, INTENT(IN)   :: JOBZ, UPLO
       INTEGER, INTENT(IN)     :: N, LDA, LWORK
       REAL(SP), INTENT(INOUT) :: SA ! array
       INTEGER, INTENT(OUT)    :: INFO
       REAL(SP), INTENT(OUT)   :: SW, SWORK
     END SUBROUTINE ssyev
     
#ifndef SYS_CRAY
       SUBROUTINE dsyev(JOBZ, UPLO, N, DA, LDA, DW, DWORK, LWORK, INFO)
       USE NumberTypes; IMPLICIT NONE
       CHARACTER, INTENT(IN)   :: JOBZ, UPLO
       INTEGER, INTENT(IN)     :: N, LDA, LWORK
       REAL(DP), INTENT(INOUT) :: DA ! array
       INTEGER, INTENT(OUT)    :: INFO
       REAL(DP), INTENT(OUT)   :: DW, DWORK
     END SUBROUTINE dsyev
#endif
  END INTERFACE

  INTERFACE syevx
       SUBROUTINE ssyevx(JOBZ, RANGE, UPLO, N, SA, LDA, SVL, SVU, IL, IU, &
            &SABTOL, NFOUND, SW, SZ, LDZ, SWORK, IWORK, IFAIL, INFO)
       USE NumberTypes; IMPLICIT NONE
       CHARACTER, INTENT(IN)   :: JOBZ, RANGE, UPLO
       INTEGER, INTENT(IN)     :: N, LDA, IL, IU, LDZ
       REAL(SP), INTENT(INOUT) :: SA ! array
       REAL(SP), INTENT(IN)    :: SVL, SVU, SABTOL
       INTEGER, INTENT(OUT)    :: NFOUND, IWORK, IFAIL, INFO
       REAL(SP), INTENT(OUT)   :: SW, SZ, SWORK
     END SUBROUTINE ssyevx
     
#ifndef SYS_CRAY
       SUBROUTINE dsyevx(JOBZ, RANGE, UPLO, N, DA, LDA, DVL, DVU, IL, IU, &
            &DABTOL, NFOUND, DW, DZ, LDZ, DWORK, LWORK, IWORK, IFAIL, INFO)
       USE NumberTypes; IMPLICIT NONE
       CHARACTER, INTENT(IN)   :: JOBZ, RANGE, UPLO
       INTEGER, INTENT(IN)     :: N, LDA, IL, IU, LDZ, LWORK
       REAL(DP), INTENT(INOUT) :: DA ! array
       REAL(DP), INTENT(IN)    :: DVL, DVU, DABTOL
       INTEGER, INTENT(OUT)    :: NFOUND, IWORK, IFAIL, INFO
       REAL(DP), INTENT(OUT)   :: DW, DZ, DWORK
     END SUBROUTINE dsyevx
#endif
  END INTERFACE
  
  INTERFACE heev
     SUBROUTINE cheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
     USE NumberTypes; IMPLICIT NONE
     CHARACTER, INTENT(IN)   :: jobz, uplo
     INTEGER, INTENT(IN)     :: lda, lwork, n
     REAL(SP), INTENT(INOUT) :: rwork, w
     COMPLEX(SP)             :: a, work
     INTEGER, INTENT(OUT)    :: info
   end SUBROUTINE cheev
#ifndef SYS_CRAY
   SUBROUTINE zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
     USE NumberTypes; IMPLICIT NONE
     CHARACTER, INTENT(IN)   :: jobz, uplo
     INTEGER, INTENT(IN)     :: lda, lwork, n
     REAL(DP), INTENT(INOUT) :: rwork, w
     COMPLEX(DP)             :: a, work
     INTEGER, INTENT(OUT)    :: info
   end SUBROUTINE zheev
#endif
  end INTERFACE
  
  INTERFACE ffti
     SUBROUTINE cffti(N, WSAVE)
       USE NumberTypes; IMPLICIT NONE
       INTEGER, INTENT(IN) :: N
       REAL(SP), INTENT(OUT) :: WSAVE
     end SUBROUTINE cffti
#ifndef SYS_CRAY
     SUBROUTINE zffti(N, WSAVE)
       USE NumberTypes; IMPLICIT NONE
       INTEGER, INTENT(IN) :: N
       REAL(DP), INTENT(OUT) :: WSAVE
     end SUBROUTINE zffti
#endif
  end INTERFACE

  INTERFACE fftf
     SUBROUTINE cfftf(n, a, wsave)
       USE NumberTypes; IMPLICIT NONE
       INTEGER, INTENT(IN)        :: n
       COMPLEX(SP), INTENT(INOUT) :: a
       REAL(SP)                   :: wsave
     end SUBROUTINE cfftf
#ifndef SYS_CRAY
     SUBROUTINE zfftf(n, a, wsave)
       USE NumberTypes; IMPLICIT NONE
       INTEGER, INTENT(IN)        :: n
       COMPLEX(DP), INTENT(INOUT) :: a
       REAL(DP)                   :: wsave
     end SUBROUTINE zfftf
#endif
  end INTERFACE
  
  INTERFACE fftb
     SUBROUTINE cfftb(n, a, wsave)
       USE NumberTypes; IMPLICIT NONE
       INTEGER, INTENT(IN)        :: n
       COMPLEX(SP), INTENT(INOUT) :: a
       REAL(SP)                   :: wsave
     end SUBROUTINE cfftb
#ifndef SYS_CRAY
     SUBROUTINE zfftb(n, a, wsave)
       USE NumberTypes; IMPLICIT NONE
       INTEGER, INTENT(IN)        :: n
       COMPLEX(DP), INTENT(INOUT) :: a
       REAL(DP)                   :: wsave
     end SUBROUTINE zfftb
#endif
  END INTERFACE

  ! ESSL ROUTINES -----------------------------------------------------
  INTERFACE cft
     SUBROUTINE scft(init, x, inc1x, inc2x, y, inc1y, inc2y, n,m, isign,&
          & scale, aux1, naux1, aux2, naux2)
       USE NumberTypes; IMPLICIT NONE
       INTEGER, INTENT(IN)        :: init, inc1x, inc2x, inc1y, inc2y, n,m
       INTEGER, INTENT(IN)        :: isign, naux1, naux2
       REAL(SP), INTENT(IN)       :: scale
       COMPLEX(SP), INTENT(INOUT) :: x(:), y(:)
       REAL(SP), INTENT(INOUT) :: aux1(:), aux2(:)
     end SUBROUTINE scft
#ifndef SYS_CRAY
     SUBROUTINE dcft(init, x, inc1x, inc2x, y, inc1y, inc2y, n,m, isign,&
          & scale, aux1, naux1, aux2, naux2)
       USE NumberTypes; IMPLICIT NONE
       INTEGER, INTENT(IN)        :: init, inc1x, inc2x, inc1y, inc2y, n, m
       INTEGER, INTENT(IN)        :: isign, naux1, naux2
       REAL(DP), INTENT(IN)       :: scale
       COMPLEX(DP), INTENT(INOUT) :: x(:), y(:)
       REAL(DP), INTENT(INOUT)    :: aux1(:), aux2(:)
     end SUBROUTINE dcft
#endif
  end INTERFACE
  
  ! BLAS ROUTINES -----------------------------------------------------
  INTERFACE GEMM
     SUBROUTINE cgemm(TRANSA, TRANSB, M, N, K, ALPHA,  A,  LDA,&
          & B, LDB, BETA, C, LDC )
       USE NumberTypes; IMPLICIT NONE
       CHARACTER, INTENT(IN)      :: TRANSA, TRANSB
       INTEGER, INTENT(IN)        :: M, N, K, LDA, LDB, LDC
       COMPLEX(SP), INTENT(IN)    :: ALPHA, BETA
       COMPLEX(SP), INTENT(INOUT) :: A,B,C        ! arrays
     END SUBROUTINE cgemm
#ifndef SYS_CRAY
     SUBROUTINE zgemm(TRANSA, TRANSB, M, N, K, ALPHA,  A,  LDA,&
          & B, LDB, BETA, C, LDC )
       USE NumberTypes; IMPLICIT NONE
       CHARACTER, INTENT(IN)      :: TRANSA, TRANSB
       INTEGER, INTENT(IN)        :: M, N, K, LDA, LDB, LDC
       COMPLEX(DP), INTENT(IN)    :: ALPHA, BETA
       COMPLEX(DP), INTENT(INOUT) :: A,B,C        ! arrays
     END SUBROUTINE zgemm
#endif
  END INTERFACE
  
  INTERFACE GEMV
     SUBROUTINE cgemv(TRANS, M, N,  ALPHA,  A,  LDA,  X,  INCX,&
          & BETA, Y, INCY )
       USE NumberTypes; IMPLICIT NONE
       CHARACTER, INTENT(IN)   :: TRANS
       INTEGER, INTENT(IN)     :: incX, incY, LDA, M, N
       COMPLEX(SP), INTENT(IN) :: ALPHA, BETA
       COMPLEX(SP), INTENT(IN) :: A, X       ! array, vector
       COMPLEX(SP), INTENT(INOUT) :: Y       ! vector

     END SUBROUTINE cgemv
#ifndef SYS_CRAY
     SUBROUTINE zgemv(TRANS, M, N,  ALPHA,  A,  LDA,  X,  INCX,&
          & BETA, Y, INCY )
       USE NumberTypes; IMPLICIT NONE
       CHARACTER, INTENT(IN)   :: TRANS
       INTEGER, INTENT(IN)     :: incX, incY, LDA, M, N
       COMPLEX(DP), INTENT(IN) :: ALPHA, BETA
       COMPLEX(DP), INTENT(IN) :: A, X       ! array, vector
       COMPLEX(DP), INTENT(INOUT) :: Y       ! vector
     END SUBROUTINE zgemv
#endif
  END INTERFACE
end MODULE LapackInterface




