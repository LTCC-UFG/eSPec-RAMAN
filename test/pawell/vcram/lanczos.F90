! lanczos.f90
! procedures for Lanczos evolution.
! Pawel Salek, pawsa@ifm.liu.se
! IMPLEMENTATION LEVEL:
! Lanczos was tested. The N must be factorizable into powers of 2, 3 and 5.
! 980617: ESSL Routines for IBM. Degraded performance for n2 problem
!         from 40s to 42s.      

MODULE Lanczos
  USE LapackInterface
  PUBLIC
  INTERFACE norm
     MODULE PROCEDURE normReal
     MODULE PROCEDURE normCplx
  end INTERFACE
CONTAINS
  ! -------------------------------------------------------------------
  ! FindGround finds the ground state in given potential.
  ! uses tridiagonal representation of hamiltonian
  ! this should be removed from this module to indep.f90
  LOGICAL FUNCTION FindGround(V, mass, dx, e0, Fi)
    IMPLICIT NONE
    REAL(mk), INTENT(IN)            :: V(:)
    REAL(mk), INTENT(IN)            :: mass, dx
    REAL(mk), INTENT(OUT)           :: e0
    REAL(mk), INTENT(OUT), OPTIONAL :: Fi(:)
    
    REAL(mk), dimension(size(V))   :: hd, ho, W
    INTEGER, dimension(size(V))    :: iblock, isplit
    REAL(mk), dimension(size(V)*5) :: work
    INTEGER,dimension(size(V)*3)   :: iwork
    INTEGER                        :: nsplit,ifail, info, N, M
    
    FindGround = .FALSE.
    N=SIZE(V)
    hd(:)=V(:)+1/(mass*(dx**2))
    ho(:)=-1/(2*mass*(dx**2))
    CALL stebz('I','B',N,1e0_mk,1e0_mk,1,1,1e-6_mk,hd(1), ho(1), M,nsplit,&
         & W(1), iblock(1), isplit(1), work(1), iwork(1), info)
    IF (info /=0 ) THEN 
       WRITE(0,*) 'FindGround: sstebz failed with info=', info
       RETURN
    ENDIF
    e0 = W(1)
    IF(PRESENT(Fi)) THEN
       CALL stein(N, hd(1), ho(1),1, W(1), iblock(1),isplit(1),Fi(1), N,&
            &work(1), iwork(1), ifail,info)
       IF (info /=0 .OR. ifail/=0) THEN 
          WRITE(0,*) 'FindGround: sstein falied with info=',info
          RETURN
       ENDIF
       Fi=Fi/norm(Fi)
       IF(ABS(MINVAL(Fi))>ABS(MAXVAL(Fi))) Fi = -Fi
    END IF
    FindGround = .TRUE.
  END FUNCTION FindGround
  
  ! =================================================================
  ! finds the ground state using the full matrix expansion of 
  ! the Hamiltonian.
  ! follows ~/HCl/Pack/findGround.m
  ! Uses packed storage.
  LOGICAL FUNCTION FindGroundFull(V, mass, dx, e0, Fi, eVec)
    IMPLICIT NONE
    REAL(mk), INTENT(IN)            :: V(:)
    REAL(mk), INTENT(IN)            :: mass, dx
    REAL(mk), INTENT(OUT), OPTIONAL :: e0
    REAL(mk), INTENT(OUT), OPTIONAL :: Fi(:)
    REAL(mk), INTENT(OUT), OPTIONAL :: eVec(:)
    
    INTEGER     :: N, Np, mi(1), loi(1), hii(1), i, j, idx, NFOUND, INFO
    REAL(mk)    :: rn, mn
    REAL(mk), ALLOCATABLE :: T(:), WORK(:)
    REAL(mk), PARAMETER   :: s(2) = (/1.0,-1.0/)
    INTEGER, ALLOCATABLE  :: IFAIL(:), IWORK2(:)

    FindGroundFull = .FALSE.
    N=SIZE(V)
    mn = minval(V); mi = minloc(V)
    rn = V(N) - mn

    ! THESE TWO ARE NOT DETERMINED OPTIMALLY....
    loi = MAX(1,MAXLOC(V(1:mi(1)), V(1:mi(1))-mn<rn*0.98))
    hii = MIN(N,MAXLOC(V(mi(1):N), V(mi(1):N)-mn<rn*0.95)+mi(1)-1)
    Np=hii(1)-loi(1)+1
    !WRITE (0,*)'FindGroundFull: from ',loi(1),' to ',hii(1),'(',Np,' points)'
    
    ALLOCATE(T((Np*Np+Np)/2), WORK(8*Np), IFAIL(Np),IWORK2(5*Np))
    ! set the matrix, first the diagonal element and then the 
    ! non-diagonal ones ...
    DO i=1,Np
       idx = (i-1)*(2*Np-(i-2))/2+1
       T(idx) = PI**2/(6*mass*(dx**2))+V(loi(1)+i-1)
       T(idx+1:idx+Np-i) = (/ (s(MOD(j,2)+1)/((j*dx)**2*mass),j=1,Np-i) /) 
    END DO
    
    IF(PRESENT(Fi)) THEN
       IF(.NOT.PRESENT(e0)) STOP "FindGroundFull: e0 not present"
       Fi = 0.0_mk      ! zeroize far elements
       CALL SPEVX('V','I','L',Np, T(1),0.0_mk,0.0_mk,1,1, -1.0_mk,&
            &NFOUND, e0, Fi(loi(1)), Np, WORK(1), IWORK2(1), IFAIL(1), INFO)
      IF( INFO /= 0 .OR. NFOUND /= 1 ) THEN
          WRITE(0,*) 'FindGroundFull - SPEVX: info=',INFO
          RETURN
       END IF
       ! make sure that the state us "up" - just for a nice visual effect
       IF(ABS(MINVAL(Fi))>ABS(MAXVAL(Fi))) Fi = -Fi
       Fi=Fi/norm(Fi)
    ELSE
       IF(PRESENT(eVec)) THEN
          CALL SPEVX('N','I','L',Np, T(1),0.0_mk,0.0_mk,1,SIZE(eVec),-1.0_mk,&
               &NFOUND, eVec(1), eVec(1), Np, WORK(1),IWORK2(1),IFAIL(1),INFO)
       ELSE
          IF(.NOT.PRESENT(e0)) STOP "FindGroundFull: no place to put energy to"
          CALL SPEVX('N','I','L',Np, T(1),0.0_mk,0.0_mk,1,1, -1.0_mk,&
               &NFOUND, e0, e0        , Np, WORK(1), IWORK2(1), IFAIL(1), INFO)
       END IF
       IF( INFO /= 0) THEN
          WRITE(0,*) 'FindGroundFull - SPEVX2: info=',INFO
          RETURN
       END IF
    END IF
    
    DEALLOCATE(T, WORK, IFAIL, IWORK2)
    FindGroundFull = .TRUE.
  END FUNCTION FindGroundFull
  

  ! -----------------------------------------------------------------
  REAL(mk) function normREAL(f)
    IMPLICIT NONE
    REAL(mk), dimension(:)    ::f
    normReal = sqrt(dot_product(f,f))
  END function normReal

  REAL(mk) function normCplx(f)
    IMPLICIT NONE
    COMPLEX(mk), dimension(:) ::f
    normCplx = sqrt(REAL(dot_product(f,f)))
  END function normCplx

  ! -----------------------------------------------------------------
  ! generates the Lanczos/Krylov subspace by finding the transf. matrix Q
  ! T(:,1) remembers diagonal elements of T, T(:,2) - upper diagonal
  ! V - the potential, F - the function
  SUBROUTINE GenKrylSp(V, F, dx, mass, Q, T, Z, D, WSAVE, Vspl)
! 			CALL GenKrylSp(Vde,Fde,dx,mass, Qd, Td, Zd, Dd, WSAVE,Vspl) !vk
    IMPLICIT NONE
    REAL(mk),INTENT(IN)      :: V(:,:)   ! the potential
    COMPLEX(mk), INTENT(IN)  :: F(:)   ! the wave function
    REAL(mk), INTENT(in)     :: dx, mass  ! space step, particle mass
    COMPLEX(mk), INTENT(out) :: Q(:,:) ! trnsf. from sp. to subsp.
    REAL(mk),INTENT(out)     :: T(:,:) ! tridiagonal H in subsp.
    REAL(mk),INTENT(OUT)     :: Z(:,:) ! the subspace basis
    REAL(mk),INTENT(out)     :: D(:) ! diagonal form of H in subsp.
    REAL(mk), INTENT(INOUT)  :: WSAVE(:)
    REAL(mk), INTENT(IN)     :: Vspl(:,:)

    INTEGER     :: N, i, info, k
    REAL(mk)    :: t1(SIZE(D)), work(2*SIZE(D)-2)
    COMPLEX(mk) :: tmp(SIZE(F))
    N=SIZE(V,1)
    Q(:,1)=F/norm(F)

    IF(SIZE(V,2) > 1) THEN
       tmp = HamOprCoupled(Q(:,1),V, mass,dx,WSAVE, Vspl); 
    ELSE
       tmp = HamOpr(Q(:,1),V(:,1), mass,dx,WSAVE); 
    ENDIF
    T(1,1)=DOT_PRODUCT(Q(:,1),tmp)
    Q(:,2)=tmp-T(1,1)*Q(:,1)
    Q(:,2)=Q(:,2)/norm(Q(:,2)) 

    IF(SIZE(V,2) > 1) THEN
       tmp = HamOprCoupled(Q(:,2),V, mass,dx,WSAVE, Vspl); 
    ELSE
       tmp = HamOpr(Q(:,2),V(:,1), mass,dx,WSAVE);
    ENDIF
    T(1,2)=DOT_PRODUCT(Q(:,1), tmp)
    T(2,1)=DOT_PRODUCT(Q(:,2), tmp)
    DO i=3,SIZE(D)
       Q(:,i)=tmp-T(i-1,1)*Q(:,i-1)-T(i-2,2)*Q(:,i-2)
       Q(:,i)=Q(:,i)/norm(Q(:,i)) 
       IF(SIZE(V,2) > 1) THEN
          tmp=HamOprCoupled(Q(:,i),V,mass,dx,WSAVE, Vspl);
       ELSE
          tmp=HamOpr(Q(:,i),V(:,1),mass,dx,WSAVE); 
       ENDIF
       T(i-1,2)=DOT_PRODUCT(Q(:,i-1), tmp)
       T(i,1) = DOT_PRODUCT(Q(:,i), tmp)
    ENDDO

    ! diagonalize H and Z
    D=T(:,1); t1 = T(:,2)
    CALL steqr('I', SIZE(D), D(1),t1(1),Z(1,1),SIZE(Z,1), work(1), info)
    if (info /= 0) then
       WRITE (0,*) 'GenKrylSp: steqr error, info=', info
       STOP
    endif
  END SUBROUTINE GenKrylSp

  ! estimStep --------------------------------------------------------
  ! estimates step of the Lanczos evolution in the subspace.
  ! Optimization introduced problems related to rounding errors. The algorithm
  ! had to be modified (the root is taken not of the product, but of 
  ! every term of the product).
  REAL(mk) FUNCTION estimStep(T, myEps)
    IMPLICIT NONE
    REAL(mk), DIMENSION(:,:)  :: T
    REAL(mk)                  :: myEps

    REAL(mk)                  :: p2
    INTEGER                   :: p, i
    p=SIZE(T,1)
    !write (0,*) 'estimStep: nondiagonal elements'
    !write (0,'(5F14.7)') ((i/T(i,2)),i=1,p-1)
    !p2 = PRODUCT( (/ ( i/T(i,2), i=1,p-1) /) )
    !estimStep = (p2*SQRT(myEps))**(1.0/(p-1))
    !p2 = PRODUCT( (/ ( sqrt(DBLE(i)/T(i,2)), i=1,p-1) /) )
    p2 = sqrt(1.0/ABS(T(1,2)))
    DO i =2, p-1
       p2 = p2*SQRT(DBLE(i)/ABS(T(i,2)))
    END DO
    estimStep = (p2*myEps)**(2.0_mk/REAL(p-1,mk))
    !write (0,*) 'p2=', p2
  end FUNCTION estimStep

  ! HamOpr -----------------------------------------------------------
  ! note that F cannot be modified!
  FUNCTION HamOpr(F, V, mass, dx, WSAVE)
    implicit none
    COMPLEX(mk), intent(IN) :: F(:)
    REAL(mk), intent(IN)    :: V(:)
    REAL(mk), INTENT(INOUT) :: WSAVE(:)
    COMPLEX(mk)             :: HamOpr(size(F))
    REAL(mk), INTENT(IN)    :: mass, dx

    INTEGER  :: i
!#if defined(SYS_CRAY)
!    REAL(mk), dimension(8*SIZE(F)) :: work
!    HamOpr=F
!    i=0;CALL ccfft(1,size(HamOpr),1.0,HamOpr(1),HamOpr(1),WSAVE(1),work(1),i)
!    HamOpr= HamOpr*(2*PI)**2/(2*mass*(dx*size(V))**2)*&
!           & (/(I*I,I=0,size(V)/2-1),(I*I,I=-size(V)/2,-1)/)
!    i=0;CALL ccfft(-1,size(HamOpr),1.0,HamOpr(1),HamOpr(1),WSAVE(1),work(1),i)
!    HamOpr=HamOpr/size(F)+ V*F
!#elif defined(SYS_ALPHA)
!    COMPLEX(mk) :: tmp(SIZE(F))
!    INTEGER ZFFT
!    i = ZFFT('C','C','B',F(1),HamOpr(1),SIZE(F),1)
!    IF(i /=0) STOP "lanczos.F90: ZFFT-B"
!    HamOpr = HamOpr*(2*PI)**2/(2*mass*(dx*SIZE(V))**2)*&
!           & (/(I*I,I=0,size(V)/2-1),(I*I,I=-size(V)/2,-1)/)
!    i = ZFFT('C','C','F',HamOpr(1),tmp(1),SIZE(F),1)
!    IF(i /=0) STOP "lanczos.F90: ZFFT-F"
!   HamOpr = tmp + V*F
!#elif defined(USE_ACML)
!    INTEGER :: N
!    N = SIZE(F)
!    HamOpr=F
!    CALL zfft1d(2,N,HamOpr,WSAVE,I)
!    IF(i /=0) STOP "lanczos.F90: ZFFT1D-B"
!    HamOpr= HamOpr*(2*PI)**2/(2*mass*(dx*size(V))**2)*&
!           & (/(I*I,I=0,size(V)/2-1),(I*I,I=-size(V)/2,-1)/)
!    CALL zfft1d(-2,N,HamOpr,WSAVE,I)
!    IF(i /=0) STOP "lanczos.F90: ZFFT1D-F"
!    HamOpr=HamOpr/size(F)+ V*F
!#elif defined(USE_ESSL)
!    REAL(mk) :: scale
!    ws=SIZE(WSAVE)/4; scale=1.0
!    HamOpr=F
!    CALL cft(0,HamOpr,1,1,HamOpr,1,1,size(F),1,1,scale,&
!         & WSAVE(1+2*ws:3*ws),ws,WSAVE(1+3*ws:4*ws), ws)
!    HamOpr= HamOpr*(2*PI)**2/(2*mass*(dx*size(V))**2)*&
!          & (/(I*I,I=0,size(V)/2-1),(I*I,I=-size(V)/2,-1)/)
!    CALL cft(0,HamOpr,1,1,HamOpr,1,1,size(F),1,-1,scale,&
!         & WSAVE(1:ws),ws,WSAVE(1+ws:2*ws), ws)
!    HamOpr=HamOpr/size(F)+ V*F
!#else
    HamOpr=F !/norm(F)!!!!!vka
!    write(0,*) size(V), norm(HamOpr)
    CALL fftf(size(HamOpr), HamOpr(1), WSAVE(1))
        HamOpr= HamOpr*(2*PI)**2/(2*mass*(dx*SIZE(F))**2)*&
           & (/(I*I,I=0,size(V)/2-1),(I*I,I=-size(V)/2,-1)/)
    CALL fftb(size(HamOpr), HamOpr(1), WSAVE(1))
    HamOpr= HamOpr/SIZE(F)+V*F
 !   write(0,*) norm(HamOpr)
!#endif
  end FUNCTION HamOpr

  ! HamOprCoupled ----------------------------------------------------
  ! note that F cannot be modified! Testing, testing...
  FUNCTION HamOprCoupled(F, V, mass, dx, WSAVE, Vspl)
    implicit none
    COMPLEX(mk), intent(IN) :: F(:)
    REAL(mk), INTENT(IN)    :: V(:,:), Vspl(:,:)
    REAL(mk), INTENT(INOUT) :: WSAVE(:)
    REAL(mk), INTENT(IN)    :: mass, dx

    INTEGER  :: i,k, r_low, r_high, c_low, c_high, N, dim
    COMPLEX(mk)             :: HamOprCoupled(SIZE(F))
    complex(mk) :: FF(size(V,1))

    N=SIZE(F)/SIZE(V,2)
	
    HamOprCoupled = (0.0, 0.0)
    dim = SIZE(V,2)
    DO k=1, dim !loop over the rows of the coupled Hamiltonian
       r_low = N*(k-1)+1   ! Indices for extraction of the relevant 
       r_high = N*k        ! region from the wave function.
       FF=F(r_low:r_high) !/norm(F(r_low:r_high))
       DO i=1, dim !loop over the  columns
 !     write(0,*) "k=",k, " i=", i, "NORM", norm(FF)
          IF (i==k) THEN
              HamOprCoupled(r_low:r_high) = HamOprCoupled(r_low:r_high) + &
                  & HamOpr(F(r_low:r_high), V(:,k), mass, dx, WSAVE)
          ELSE
  !        write(0,*) "ind", Ind(k,i,dim)
             c_low = N*(i-1)+1
             c_high = N*i
   !          write(0,*) "c:", c_low, c_high, norm(F(c_low:c_high))
    !         write(0,*) "r:", r_low, r_high, norm(HamOprCoupled(r_low:r_high))
             HamOprCoupled(r_low:r_high) = HamOprCoupled(r_low:r_high) + &
                  & Vspl(:,Ind(k,i,dim))*F(c_low:c_high)
     !        write(0,*) "c:", c_low, c_high, norm(HamOprCoupled(r_low:r_high))
          ENDIF
       ENDDO
    ENDDO
  END FUNCTION HamOprCoupled
  
  ! Function that converts the indices of the matrix 
  ! (i.e. of the coupled Hamiltonian) to a number
  ! (i.e to the index of the coupling vector from the input). 
  ! Check /home/yasve/pkg/ram/victor/test.m
  INTEGER FUNCTION Ind(i, j, dim)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i, j, dim
    
    INTEGER  :: ii, mini, shift
    
    mini = MIN(i,j)
    shift = 0
    DO ii=1, mini - 1
       shift = shift + MOD(dim - ii, dim)
    ENDDO

    Ind = ABS(i - j) + shift
   END FUNCTION Ind

end MODULE Lanczos

