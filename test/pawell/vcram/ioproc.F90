! ioproc.F90 is the core of input processing and program initialization.
! separates also typical sequential operations from the main driver 
! module.
!
! Pawel Salek, pawsa@theochem.kth.se
! 980226
! GLOBAL variables:
! outUnit - an unit the output is sent to.
! dumpPot    - dump potentials
! dumpGround - dump ground state initial wave function.
! dumpDOm    - dump |d(omega)>
! dumpFFT    - dump FFT stages of autocorrelation function.
! dumpQoper  - dump Q(R) operator.
! dumpAvgOps - dump average values of r and p, <r> and <p>.
! masterP - for MPP codes: if the instance should output anything.
! dumps are always off for slave PEs -except for Duration, which is
! gathered by the host and printed collectively.
! domOnly - only |d(omega)> is evaluated.
!
! CHANGE LOG:
! 990302: Raman-2
! 990308: parameters may directly follow keywords in the same line
! 990317: spherical bessel functions for Q
! 000822: CextFact added. This is to allow longer final state life time
!         but still reasonable size of data files. For no extension
!         it is 0, but sometimes it is useful to have 10 times longer
!         final life times than it would come out of energy scale 
!         calculations.
MODULE IOProcessor
  USE LapackInterface
  USE Lanczos
  SAVE
  INTEGER,  PARAMETER   :: outUnit = 0
  REAL(mk), PARAMETER   :: au2ev=27.2114, ev2au = 1.0/27.2114
  REAL(mk), PARAMETER   :: cm2au = 1.239842e-4*ev2au
  REAL(mk), PARAMETER   :: au2fs = 2.418884e-2
  
  ! ==================================================================
  ! TYPE DEFINITIONS
  ! ==================================================================
  ! spectrum%sigma(E,omega)
  TYPE spectrum
     REAL(mk) :: zeroInc, zeroOut, deltaOm
     INTEGER  :: inStep
     REAL(mk), POINTER :: sigma(:,:)
     LOGICAL  :: modified
  END TYPE spectrum
  
  INTEGER, PARAMETER :: ModeAbsorp = 1, ModeRaman = 2, ModeDecay = 3,&
       & ModeRamDir = 4, ModeRamSym = 5
  ! ==================================================================
  ! CALCULATION STATE
  ! ==================================================================
  INTEGER        :: Mode
  REAL(mk), POINTER  :: Vgr(:), Vde(:,:), Vfi(:,:), DOper(:,:), Vspl(:,:)
  COMPLEX(mk), POINTER :: QOper(:,:)                 ! Q operator
  REAL(mk)       :: MinX, MaxX, mass, Gamma, DecayThreshold, vcdip
  REAL(mk)       :: FinalGamma, FinalDecayThr
  CHARACTER      :: KinOrBin
  REAL(mk)       :: Emin, Emax, gPos, gWidth, StepPrec, convWidth, abscf
  INTEGER        :: EptNum, gStep, N, verbLvl, SubSize, CextFact, HomoParity
  INTEGER        :: absln
  COMPLEX(mk)    :: ACoeff
  CHARACTER(256) :: ofName='fort.19'
  LOGICAL        :: masterP=.TRUE.
  LOGICAL        :: dumpPot, dumpGround, dumpDOm, dumpFFT, dumpDur, dumpCoreE0
  LOGICAL        :: dumpFinalE0, dumpAvgOps = .FALSE., dumpSigmaT = .FALSE.
  LOGICAL        :: dumpQoper, DumpDoper, dumpSplit
  LOGICAL        :: optTestDissociation = .FALSE., optEvPotentials = .FALSE.
  LOGICAL        :: optE0Shifted = .FALSE., optNormalize = .FALSE.
  LOGICAL        :: optPrintCoreDumped = .FALSE.
  ! =================================================================
  REAL(mk), ALLOCATABLE    :: Fi(:)
  REAL(mk)                 :: E0

CONTAINS
  ! ===================================================================
  !               SPECTRUM HANDLING PROCEDURES
  ! ===================================================================
  ! the vector lengths are computed from requested accuracies.
  SUBROUTINE initSpectrum(spec, omCenter, omWidth, inMult)
    IMPLICIT NONE
    TYPE(spectrum), INTENT(OUT) :: spec
    REAL(mk), INTENT(IN)        :: omCenter, omWidth
    INTEGER, INTENT(IN)         :: inMult

    INTEGER :: clen, inCnt
    LOGICAL, SAVE :: doneOnce = .FALSE.

    ! find the data in the following sequence:
    clen = 2**nextpow2(4*EptNum)
    spec%deltaOm = (EMax-EMin)/REAL(EptNum)
    IF (omWidth==0.0) THEN
       inCnt = 1; spec%zeroInc = omCenter
       IF(.NOT.doneOnce) &
            &WRITE (0,*) 'WARNING: gWidth=0 -> narrow band calculation.'
    ELSE
       inCnt=NINT(2.0*gWidth/(inMult*spec%deltaOm))
       IF(inCnt <1) THEN
          WRITE(0,*)"spectrumInit: inCnt<1: either deltaOm or gStep too big.",&
               &" inCnt set to 1."
          inCnt = 1
       END IF
       spec%zeroInc = omCenter-inMult*spec%deltaOm*(inCnt-1)*0.5
    END IF
    spec%inStep = inMult
    
    ALLOCATE(spec%sigma(clen, inCnt))
    spec%zeroOut = INT((0.5*(EMax+EMin)-0.75*clen*spec%deltaOm)/spec%deltaOm)&
         & *spec%deltaOm
    IF(spec%zeroOut>EMin) THEN
       WRITE(0,"('spectrumInit: too low accuracy suspected. EMin=',F9.3,"//&
            &"' and computed zeroInc=',F9.3)") EMin, spec%zeroOut
       RETURN; !STOP "spectrumInit"
    END IF
    
    IF(.NOT.doneOnce .AND.masterP)&
         &WRITE(0,"('Spectra set alloc''d. Incoming energy range: ',F9.3,"//&
         &"' to ',F9.3,'eV, ',I3,' rows.',/,'Out energy range: (',"//&
         &"F9.3,', ',F9.3,') eV, ',I5,' points, dOmega=',F8.5)") &
         &spec%zeroInc/ev2au,&
         & (spec%zeroInc+inMult*spec%deltaOm*(inCnt-1))/ev2au, inCnt,&
         & spec%zeroOut/ev2au,(spec%zeroOut+spec%deltaOm*(clen-1))/ev2au,&
         &clen, spec%deltaOm/ev2au
    doneOnce = .TRUE.
    spec%modified = .FALSE.
  END SUBROUTINE initSpectrum

  ! placeSpectrum:
  ! places given spectrum C in pos column of spec, taking care of
  ! possible contractions from expCoef neighboring points in C to 1 in
  ! spec, and shifting also the array in the energy space.
  SUBROUTINE placeSpectrum(spec, pos, C, omZero)
    IMPLICIT NONE
    TYPE(spectrum), INTENT(INOUT) :: spec    ! C is placed in spec
    INTEGER, INTENT(IN)           :: pos     ! the column index in spec
    REAL(mk), INTENT(IN)          :: C(:)    ! the spectrum to be placed
    REAL(mk), INTENT(IN)          :: omZero  ! the freq. of C(1)

    INTEGER  :: shift, halfE, sz, l, u, i, sc

    ! contract...
    spec%sigma(:,pos) = 0
    halfE = CextFact/2
    sz = SIZE(spec%sigma,1)
    sc = SIZE(C)
    
    IF(sz*CextFact /= sc) STOP "Array size mismatch in placeSpectrum"
    DO i = 1, sz
       l = MAX(1, (i-1)*CextFact-halfE)
       u = MIN(sc,(i  )*CextFact-halfE)
       spec%sigma(i,pos) = SUM( C(l:u) )/(u-l+1)
    END DO

    ! calculate the shift and wrap it right...
    shift = NINT((spec%zeroOut-omZero)/spec%deltaOm)
    shift = MOD(shift-1, sz)+1

    ! ..and shift...
    spec%sigma(:,pos) = CSHIFT(spec%sigma(:,pos), shift)
    IF(optNormalize) spec%sigma(:,pos) = spec%sigma(:,pos)&
         & /SQRT( SUM(spec%sigma(:,pos)**2) )
    spec%modified = .TRUE.
  END SUBROUTINE placeSpectrum

  ! saveSpectrum -----------------------------------------------------
  ! saves a spectrum
  ! does the energy convolution, if the convWidth is set.
  SUBROUTINE saveSpectrum(spec, oFNm, mJ)
    IMPLICIT NONE
    TYPE(spectrum), INTENT(INOUT):: spec    ! modified modified
    CHARACTER(*), INTENT(IN)     :: oFNm
    INTEGER, INTENT(IN),OPTIONAL :: mJ

    INTEGER               :: Clen, i, j, k, maxJ
    REAL(mk), ALLOCATABLE :: convMat(:,:), vec(:)
    CHARACTER(40)         :: myFmt
    real(mk), allocatable :: vecR(:,:)
    
    maxJ = SIZE(spec%sigma,2)
    IF(PRESENT(mJ)) maxJ = mJ
    
    if(.NOT. spec%modified) RETURN
    IF(.NOT. ASSOCIATED(spec%sigma) ) STOP "saveSpec: spec%sigma not assigned."
    ALLOCATE(convMat(maxJ,maxJ), vec(maxJ))

    WRITE(0,*) "saveSpectrum: convWidth:  ", convWidth
    IF(convWidth>0) THEN
       convMat(1:maxJ,1) = EXP( -LOG(2.0_mk)*(/ ( &
            & (i-1)*spec%deltaOm*spec%inStep/convWidth, i=1, maxJ) /)**2 )
       convMat(1:maxJ,1) = convMat(1:maxJ,1)/norm( convMat(1:maxJ,1) )
       WRITE(0,"('Conv matrix:',/,(5F10.7))")  convMat(1:maxJ,1)
    ELSE
       convMat(1:maxJ,1) = (/ 1.0_mk, (0.0_mk, i=2,maxJ) /)
    END IF
    convMat(1,2:maxJ) = convMat(2:maxJ,1)
    DO j=2, maxJ
       convMat(j:maxJ,j) = convMat(1:maxJ-j+1,1)
       convMat(j,j+1:maxJ) = convMat(j+1:maxJ,j)
       convMat(1:maxJ,j)   = convMat(1:maxJ,j)/norm( convMat(1:maxJ,j) )
    END DO
    
    CLen = SIZE(spec%sigma,1)
    allocate(vecR(clen,maxj))
    ! open and write according to the new algorithm...
    ! ------------------------------------------------
    i=MAX(80,(maxJ+2)*12)
    OPEN(9,STATUS='REPLACE',ACTION='WRITE',RECL=i,FILE=oFNm)
    IF(maxJ==1) WRITE(9, FMT='(A1)', ADVANCE="NO") '#'

    WRITE(9,FMT="(I8,"//i2s(maxJ)//"(F12.4))")0,&
         & ((spec%zeroInc+i*spec%deltaOm*spec%inStep)/ev2au, i=0, maxJ-1)
    myFmt="(F8.3,"//i2s(maxJ)//"(E12.4))"
    
    DO i=1, CLen
       IF(KinOrBin == 'K') THEN
          vec = MATMUL( convMat, spec%sigma(i,1:maxJ) )
       ELSE
          ! ugly hack because the vectors corresponding to same excitation
          ! energy are shifted - see placeSpectrum.
          DO j=1, maxJ
             vec(j) = 0.0
             DO k = MAX( (j*spec%inStep -i+spec%inStep)/spec%inStep,1),&
                  & MIN( (j*spec%inStep+CLen-i)/spec%inStep, maxJ)
                vec(j) = vec(j) +convMat(k,j)*spec%sigma(i+(k-j)*spec%inStep,k)
             END DO
          END DO
       END IF
	     WRITE(9,FMT=myFmt) (spec%zeroOut+(i-1)*spec%deltaOm)/ev2au, vec
	  END DO
    CLOSE(UNIT=9)
    spec%modified = .FALSE.
  END SUBROUTINE saveSpectrum
  
  
  ! ===================================================================
  !               END of SPECTRUM HANDLING PROCEDURES
  ! ===================================================================
  SUBROUTINE confFreeMem
    DEALLOCATE(Fi, Vgr)
    IF(ASSOCIATED(Vde)   ) DEALLOCATE(Vde)
    IF(ASSOCIATED(Vfi)   ) DEALLOCATE(Vfi)
    IF(ASSOCIATED(QOper) ) DEALLOCATE(QOper)
    IF(ASSOCIATED(DOper) ) DEALLOCATE(DOper)
    IF(ASSOCIATED(Vspl)  ) DEALLOCATE(Vspl)	
  END SUBROUTINE confFreeMem
  
  ! ===================================================================
  ! reqAbsorp - check if all the data necessary for the absorbtion has been
  ! given.
  LOGICAL FUNCTION areReqMetP(augEne)
    IMPLICIT NONE
    REAL(mk), INTENT(IN) :: augEne
    REAL(mk)  :: recT, realT

    areReqMetP = .FALSE.
    IF(EMax<=Emin) THEN
       IF(masterP) WRITE (outUnit,*) &
            &"The energy range is not properly specified (KINET or BINDE)."
       RETURN
    END IF
    IF (.NOT.ASSOCIATED(Vgr)) THEN
       IF(masterP) WRITE (outUnit,*) &
            &"You can't even dream of starting without PES of initial state."
       RETURN
    END IF
    IF (.NOT.ASSOCIATED(Vde).AND.Mode/=ModeAbsorp) THEN
       IF(masterP) WRITE (outUnit,*) &
            &"Cannot do without an intermediate state."
       RETURN
    END IF
    IF(.NOT.ASSOCIATED(Vfi).AND.Mode/=ModeDecay) THEN
       IF(masterP) WRITE (outUnit,*) &
            &"Do not know a PES of the final state."
       RETURN
    END IF
    IF (augEne == 0 .AND.Mode == ModeRamSym) THEN
       IF(masterP) WRITE (outUnit,*) &
            &"Homonuclear molecule calculation (RAMSYM) requested but "//&
            &"no QOPER speficied."       
       RETURN
    END IF
    IF( ACoeff.EQ.0.0 .AND. Mode.EQ.ModeRamDir) THEN
       IF(masterP) WRITE (outUnit,*) &
            &"A=0 - it would be just an ordinary Raman (Use DIRCF)"
       RETURN
    END IF

    ! just WARNINGS and DEFAULTS ........
    IF(.NOT.ASSOCIATED(Vde).AND.dumpCoreE0) THEN
       IF(masterP) WRITE (outUnit,'(A)') &
       &"WARNING: Cannot dump the lowest vibr. energy of undefined core PES."
       dumpCoreE0 = .FALSE.
    END IF
    IF(.NOT.ASSOCIATED(Vfi).AND.dumpFinalE0) THEN
       IF(masterP) WRITE (outUnit,'(A)') &
       &"WARNING: Cannot dump the lowest vibr. energy of undefined final PES."
       dumpFinalE0 = .FALSE.
    END IF

    IF(Mode == ModeRamSym.AND. HomoParity ==0) THEN
       IF(masterP) WRITE (outUnit,'(A)') &
       &"DEFAULT: Parity for RamSym calculation not defined: defaulting to 1."
       HomoParity = 1
    END IF

    ! This is very crude estimation of the recommended time...
    IF(Mode.NE.ModeAbsorp) THEN
       realT = 2**(nextpow2(4*EptNum)-1)/(EMax-EMin)
       recT  = -LOG(DecayThreshold)/Gamma
       IF(recT>realT)THEN
          IF(masterP) WRITE (outUnit,"('The recommended final life time ',"//&
               &"F11.5,' is smaller than real: ',F11.5,/," //&
               &"'The Fourier transform is VERY LIKELY to be wrong.')") &
               &recT,realT
       END IF
    END IF
    areReqMetP = .TRUE.
  END FUNCTION areReqMetP

  ! initQOper ---------------------------------------------------------
    ! set up the Q operator, convert augEnergy to au and then to k vector...
    ! ====================
  SUBROUTINE initQOper(qOrder, qType, qAlpha, momentum,&
       & rho, rc, deltar, dx)
    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: qOrder
    CHARACTER, INTENT(IN)   :: qType
    REAL(mk), INTENT(IN)    :: qAlpha, momentum, rho, rc, deltar, dx
    REAL(mk)                :: x
    INTEGER                 :: i
    
    IF(ASSOCIATED(QOper)) STOP "ERROR: Make up your mind: QOPER or GAMVC"
    ALLOCATE(QOper(SIZE(Vgr),1))
    WRITE(0,*)'QOperator, k=',momentum, ' alpha=',qAlpha
    
    DO i=1, SIZE(Qoper,1)
       x = (MinX+ (i-1)*dx)
       SELECT CASE(qType)
       CASE('A')
          x = qAlpha*momentum*x
          QOper(i,1) = SIN(x-REAL(qOrder,mk)*PI*0.5)/x
       CASE('B')
          STOP "This option is temporarily disabled in ioproc.F90"
           !CALL DBESJ(x,0.5_mk+qOrder, 1,QOper(i), NZ) 
           !QOper(i) = QOper(i)*SQRT(0.5*PI/x)
        CASE('S')
           QOper(i,1) = EXP( CMPLX(0.0_mk, -qAlpha*momentum*x) )*&
                &(1.0+rho*(0.5+1.0/PI*ATAN((x-rc)/deltar)))
        CASE DEFAULT
           QOper(i,1) = EXP( CMPLX(0.0_mk, -qAlpha*momentum*x) )
        END SELECT
     END DO
   END SUBROUTINE initQOper

  ! -------------------------------------------------------------------
  !     initdata
  !
  ! the routine responsible for reading parameters from the stdio
  ! (which is most probably redirected from a file).
  ! Every parameter is preceeded by its name.
  ! assumes that it is called only once per instance and potentials are
  ! not yet allocated.
  ! KEYWORDS ARE ORDERED (almost) APLHABETICALLY.
  LOGICAL FUNCTION initdata(iUnit)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: iUnit
  
    CHARACTER(LEN=6) :: WORD
    REAL(mk)         :: aMod, aPhase, augEne, qAlpha, dx, xAvg, t(10), x
    REAL(mk)         :: SWITCH_rho, SWITCH_rc, SWITCH_deltar
    REAL(mk)         :: momentum, theta
    INTEGER          :: i, j, pos(1), qOrder, NZ, Num_spl, count
    CHARACTER        :: qType
    REAL(mk), POINTER :: tmp(:), qtmp(:,:)   ! qtmp for accomulating QOpers

    !     DEFAULTS ------------------------------------------------------
    initdata = .FALSE.
    NULLIFY(Vgr,Vde,Vfi,tmp,QOper,DOper,qtmp) 
    mass=-1.0
    N=-1; minX=1; maxX=-1
    SubSize = 16
    StepPrec = 1e-6
    FinalGamma = -1.0 ! disabled
    Gamma=1.0; Emin=0.0; Emax=0.0; EptNum=64
    verbLvl = 2; gStep=1
    absln = 5; abscf = 0.0007
    ACoeff = 0.0

    augEne = 0.0    ! disabled
    convWidth = 0.0 ! convolution disabled
    CextFact = 1    ! no extra extension of the final state lifetime
    dumpQoper = .FALSE.
    dumpDoper = .FALSE.
    dumpSplit = .FALSE.
    HomoParity = 0  ! use default value (see below)

    ! READ caclulation type: ABSORBtion, RAMAN scattering, RAMan with 
    ! a DIRect term.
    TypeLoop: DO
       READ(iUnit,'(A6)', END=222,ERR=222) WORD
       IF (WORD(1:1).NE.'#'.AND.WORD(1:1).NE.'!') EXIT TypeLoop
    END DO TypeLoop
    SELECT CASE(WORD)
    CASE("ABSORP")
       Mode = ModeAbsorp
       IF(masterP) WRITE (outUnit,*) "ABSORPTION to Vfinal computed"
    CASE("RAMANS")
       Mode = ModeRaman
       IF(masterP) WRITE (outUnit,*) "RAMANS scattering calculation"
    CASE("FDECAY")
       Mode = ModeDecay
       IF(masterP) WRITE (outUnit,*) "DECAY probability calculation"
    CASE("RAMDIR")
       Mode = ModeRamDir
       IF(masterP) WRITE (outUnit,*) "RAMAN scattering calculation",&
            &" with DIRect term included"
    CASE("RAMSYM")
       Mode = ModeRamSym
       IF(masterP) WRITE (outUnit,*) "RAMAN scattering calculation",&
            &" for a homonuclear molecule"
    CASE DEFAULT
       IF(masterP) WRITE (0,*) "Do not know calculation mode '",TRIM(WORD),&
            &"'. Known are ABSORP - absorption, RAMANS - Raman ",&
            &"scattering, FDECAY, RAMDIR - RAMan with DIRect term"
       RETURN
    END SELECT
    
    ! read other keywords...
    DO
!       WRITE(0,*) 'Reading keyword'
       READ(iUnit,'(A5)', ADVANCE="NO", END=222,ERR=222) WORD
       IF (WORD(1:1).EQ.'#' .OR.WORD(1:1).EQ.'!') THEN
          !WRITE(0,*) 'Skipped comment ',WORD
          READ(iUnit,"(A1)") WORD(1:1)
          CYCLE
       END IF
       SELECT CASE(WORD)
       CASE ("MASS ")
          READ(iUnit,*, ERR=101,END=101) mass
          IF(masterP) WRITE (outUnit,*) 'Read  mass ',  mass
       CASE("MAMU2")
          READ(iUnit,*, ERR=101,END=101) t(1), t(2)
          mass = t(1)*t(2)/(t(1)+t(2)) * 1836.15270
          IF(masterP) WRITE (outUnit,"('Atoms have masses ',F8.4,' and ',"//&
               &"F8.4,' a.m.u. Effective mass: ',F9.2,' a.u.')") t(1),t(2),mass
       CASE("DLINE")
          CALL loadLine(iUnit,"DLINE",tmp,'EV')
          CALL IncrDimension(Vde, tmp);DEALLOCATE(tmp)
       CASE("FLINE")
          CALL loadLine(iUnit,"FLINE",tmp,'EV')
          CALL IncrDimension(Vfi, tmp);DEALLOCATE(tmp)
       CASE("IHARM") 
          CALL loadHarm(iUnit,"IHARM",Vgr)
       CASE("DHARM")
          CALL loadHarm(iUnit,"DHARM",tmp)
          CALL IncrDimension(Vde, tmp);DEALLOCATE(tmp)
       CASE("FHARM")
          CALL loadHarm(iUnit,"FHARM",tmp)
          CALL IncrDimension(Vfi, tmp);DEALLOCATE(tmp)
       CASE ("IMORS") 
          CALL loadMors(iUnit,"IMORS",Vgr)
       CASE ("DMORS") 
          CALL loadMors(iUnit,"DMORS",tmp)
          CALL IncrDimension(Vde, tmp);DEALLOCATE(tmp)
       CASE ("FMORS") 
          CALL loadMors(iUnit,"FMORS",tmp)
          CALL IncrDimension(Vfi, tmp);DEALLOCATE(tmp)
       CASE("IBUCK")
          CALL loadBuck(iUnit,"IBUCK",Vgr)
       CASE("DBUCK")
          CALL loadBuck(iUnit,"DBUCK",tmp)
          CALL IncrDimension(Vde, tmp);DEALLOCATE(tmp)
       CASE("FBUCK")
          CALL loadBuck(iUnit,"FBUCK",tmp)
          CALL IncrDimension(Vfi, tmp);DEALLOCATE(tmp);
       CASE("IOMAR")
          CALL loadMorseOmar(iUnit,"IOMAR",Vgr)
       CASE("DOMAR")
          CALL loadMorseOmar(iUnit,"DOMAR",tmp)
          CALL IncrDimension(Vde, tmp);DEALLOCATE(tmp)
       CASE("FOMAR")
          CALL loadMorseOmar(iUnit,"FOMAR",tmp)
          CALL IncrDimension(Vfi, tmp);DEALLOCATE(tmp)
       CASE("ICMOR")
          CALL loadMorseCm(iUnit,"ICMOR",Vgr)
       CASE("DCMOR")
          CALL loadMorseCm(iUnit,"DCMOR",tmp)
          CALL IncrDimension(Vde, tmp);DEALLOCATE(tmp)
       CASE("FCMOR")
          CALL loadMorseCm(iUnit,"FCMOR",tmp)
          CALL IncrDimension(Vfi, tmp);DEALLOCATE(tmp)
       CASE("IHULH")
          CALL loadHulburtHirschfelder(iUnit,"IHULH",Vgr)
       CASE("DHULH")
          CALL loadHulburtHirschfelder(iUnit,"DHULH",tmp)
          CALL IncrDimension(Vde, tmp);DEALLOCATE(tmp)
       CASE("FHULH")
          CALL loadHulburtHirschfelder(iUnit,"FHULH",tmp)
          CALL IncrDimension(Vfi, tmp);DEALLOCATE(tmp)
       CASE("ISPLI")
          CALL loadSpline(iUnit,'ISPLI',Vgr)
       CASE("DSPLI")
          CALL loadSpline(iUnit,'DSPLI',tmp)
          CALL IncrDimension(Vde, tmp);DEALLOCATE(tmp)
       CASE("FSPLI")
          CALL loadSpline(iUnit,'FSPLI',tmp)
          CALL IncrDimension(Vfi, tmp);DEALLOCATE(tmp)
       CASE("CONVO") ! Convolution of final cross section
          READ(iUnit,*,ERR=101,END=101) convWidth
          IF(masterP) THEN
             WRITE(outUnit,"('Cross section convolution width:',F8.4,' eV')")&
                  & convWidth
             IF(Mode == ModeDecay .OR. Mode == ModeAbsorp) &
                  & WRITE (outUnit,"(A)") &
                  & "Cross section not available. Convolution will NOT be"//&
                  & " performed."
          END IF
          convWidth = convWidth * ev2au
       CASE("DIRCF") ! A - direct transition strength coefficient
          READ(iUnit,*,ERR=101,END=101) aMod, aPhase
          ACoeff = aMod*EXP(CMPLX(0.0_mk, aPhase*PI/180.0))
          IF(masterP) THEN
             IF(Mode /= ModeRamDir) WRITE (outUnit,"(A)",ADVANCE='NO') &
                  &"A coefficient IGNORED"
             WRITE(outUnit,*) ' Read DIRcoef=',ACoeff
          END IF
       CASE("DUMPS")
          CALL processDumps(iUnit)
       CASE("BINDE")
          READ(iUnit,*, ERR=101,END=101) Emin, Emax, EptNum
          IF(masterP) WRITE (outUnit,"('Binding energy, Emin:',F9.3,"//&
               &"' Emax:',F9.3,' (eV) scanned at ',I4,' points,')") &
               &Emin, Emax, ePtNum
          Emin=Emin*ev2au; Emax=Emax*ev2au
          KinOrBin = 'B'
       CASE("KINET")
          READ(iUnit,*, ERR=101,END=101) Emin, Emax, EptNum
          IF(masterP) WRITE (outUnit,"('Electron kinet.En., Emin:',F9.3,"//&
               &"' Emax:',F9.3,' (eV) scanned at ',I4,' points,')") &
               &Emin, Emax, ePtNum
          Emin=Emin*ev2au; Emax=Emax*ev2au
          KinOrBin = 'K'
       CASE("FIGAM")
          READ(iUnit,*, ERR=101,END=101) FinalGamma, FinalDecayThr
          x = -LOG(FinalDecayThr)/(FinalGamma*ev2au)
          IF(masterP) WRITE (outUnit,"('FORCING Final Gamma ',F10.5, "//&
               &"' eV DecayThreshold:',E10.3,/,8X,'Corresponding life "//&
               &"time:',F13.4,' au=',F10.3,' fs')") FinalGamma, FinalDecayThr,&
               & x, x*au2fs
          FinalGamma=FinalGamma*ev2au
       CASE("GAMMA")
          READ(iUnit,*, ERR=101,END=101) Gamma, DecayThreshold
          IF(masterP) WRITE (outUnit,"('Read Gamma ',F10.5, "//&
               &"'eV DecayThreshold ',E10.3)") Gamma, DecayThreshold
          Gamma=Gamma*ev2au
       CASE("GAMVC")
          NULLIFY(tmp); CALL loadSpline(iUnit,'GAMVC', tmp)
          CALL IncrDimension(qtmp, tmp); DEALLOCATE(tmp)
       CASE("GAMVL")
          NULLIFY(tmp); CALL loadLine(iUnit,'GAMVC', tmp, "au")
          CALL IncrDimension(qtmp, tmp); DEALLOCATE(tmp)
       CASE("SPLIT")
          NULLIFY(tmp); CALL loadSpline(iUnit,'SPLIT', tmp)
          CALL IncrDimension(Vspl, tmp); DEALLOCATE(tmp)
       CASE("SPLIL")
          NULLIFY(tmp); CALL loadLine(iUnit,'SPLIL', tmp, "au")
          CALL IncrDimension(Vspl, tmp); DEALLOCATE(tmp)
!vk differend Auger decays for vc-final states
       CASE("VCDIP")
          READ(iUnit,*, ERR=101,END=101) vcdip
!vk add Gaussian for splittings
       CASE("SPLIG")
          NULLIFY(tmp); CALL loadGauss(iUnit,'SPLIG', tmp, "au")
          CALL IncrDimension(Vspl, tmp); DEALLOCATE(tmp)      
       CASE("DIPOL")
          NULLIFY(tmp); CALL loadSpline(iUnit,'DIPOL', tmp)
          CALL IncrDimension(DOper, tmp); DEALLOCATE(tmp)
       CASE("DIPLI")
          NULLIFY(tmp); CALL loadLine(iUnit,'DIPOL', tmp, "au")
          CALL IncrDimension(DOper, tmp); DEALLOCATE(tmp)
       CASE("OFNAM")
          READ(iUnit,*, ERR=101,END=101) ofName
          IF(masterP)WRITE(outUnit,"('Output base file name: ',A)")TRIM(ofName)
       CASE("OMINC")
          READ(iUnit,*, ERR=101,END=101) gPos, gWidth, gStep
          IF(masterP) WRITE (outUnit,"('gPos:',F12.4,' eV gWidth:',"//&
               &"F12.4,' eV, gStep:',I4,' freq. units (dOmega)')")&
               & gPos,gWidth,gStep
          gPos=gPos*ev2au; gWidth=gWidth*ev2au
       CASE("OPTNS")
          CALL processOptions(iUnit)
       CASE("QOPER")
          READ(iUnit,*, ERR=101,END=101) AugEne, qAlpha, qOrder, qType
          IF(masterP) WRITE (outUnit,"('Q operator set up, assumed E_A=',"//&
               &"F9.3,' alpha mass coef=',F7.3,' ord=',I3,' -',A1,'-')") &
               &AugEne, qAlpha, qOrder, qType
          momentum = SQRT(2*AugEne*ev2au)
       CASE("QEXPT")
          qType  = 'E' ! exp(-ikr cos(theta))
          qOrder = 0   ! ignored
          READ(iUnit,*, ERR=101,END=101) AugEne, qAlpha, theta
          IF(masterP) WRITE (outUnit,"('Q =exp(-ikr*cos(t)) set up, E_A=',"//&
               &"F9.3,' alpha mass coef=',F7.3,' t=',F7.2,' degrees')") &
               &AugEne, qAlpha, theta
          momentum = SQRT(2*AugEne*ev2au)*COS(theta*PI/180.0)
       CASE("SWITC")
          qType  = 'S' ! exp(-ikr cos(theta))*(1+rho*switch(r;R_c,deltaR))
          qOrder = 0   ! where switch=0.5+1/pi*arctan((r-R_c)/deltaR)
          READ(iUnit,*, ERR=101,END=101) AugEne, qAlpha, theta,&
               &SWITCH_rho, SWITCH_Rc, SWITCH_deltar
          IF(masterP) &
               &WRITE (outUnit,&
               &"('Q =exp(-ikr*cos(t))*switch(R) set up,',/,'"//&
               &"  E_A=',F9.3,' eff. mass=',F7.3,' t=',F7.2,' degrees',"//&
               &"/,'  Switch parameters: rho=',F7.3,' Rc=',F6.3,' DeltaR='"//&
               &",F9.3)") &
               &AugEne, qAlpha, theta,SWITCH_rho, SWITCH_Rc, SWITCH_deltar
          momentum = SQRT(2*AugEne*ev2au)*COS(theta*PI/180.0)          
       CASE("PARIT")
          READ(iUnit,*, ERR=101,END=101) HomoParity
          IF(masterP) WRITE (outUnit,"('Parity factor for RAMSYM "//&
               &" calculation: ',I2)") HomoParity
       CASE("SPACE")
          READ(iUnit,*, ERR=101,END=101) minX, maxX, N
          dx=(maxX-minX)/(N-1)
          IF(masterP) WRITE (outUnit,"('Read minX:',F7.3,' au maxX:',F7.3,"//&
               &"' au defined on ',I4,' points, dx=',F6.3)")  minX, maxX, N, dx
          IF( minX>=maxX .OR. N<2 ) THEN
             WRITE (0,*) 'minX >= maxX. or N < 2. STOP'
             RETURN
          END IF
!vk test Absorbin gboundary conditions
       case("ABSVC")
       		read(iUnit,*,ERR=101,END=101) absln, abscf
    			IF(masterP) THEN
             WRITE(outUnit,"('Read ABSVC: length=',I5,'/',I2,', coeff=',E9.3)") &
             & N, absln, abscf
          END IF   
!vk           
       CASE("SUBSZ")
          READ(iUnit,*, ERR=101,END=101) SubSize
          IF(masterP) WRITE (outUnit,*) 'Read SubSize ',  SubSize
       CASE("SPREC")
          READ(iUnit,*, ERR=101,END=101) StepPrec
          IF(masterP) WRITE (outUnit,"(' Read StepPrec ',E8.2)") StepPrec
       CASE("VELEV")
          READ(iUnit,*, ERR=101,END=101) verbLvl
       CASE DEFAULT
          IF(masterP) WRITE(outUnit,"((A))")&
               & 'Unknown keyword (',WORD,') encountered. Known are:',&
               & "MASS MAMU2 IHARM DHARM FHARM ISPLI DSPLI FSPLI GAMMA KINET",&
               & " BINDE OMINC SPACE SUBSZ SPREC VELEV OFNAM DUMPS SPLIS SPLIL"
          RETURN
     END SELECT
  ENDDO
222 CONTINUE
  
  IF (ASSOCIATED(qtmp)) THEN !Needed because QOper is complex pointer
     ALLOCATE(QOper(SIZE(qtmp, 1), SIZE(qtmp, 2)))
     QOper = qtmp;  DEALLOCATE(qtmp)
  ENDIF

  !Check if all defined decay potentials have corresponding DOpers or QOpers
  !Check if all finals states have corresponding DOpers.
  !This checks must be organized differently.
  !FIXME:
  IF (ASSOCIATED(Vde)) THEN
     WRITE(outUnit, *) "Number of decay states: ", SIZE(Vde,2)
     IF(ASSOCIATED(QOper) .AND. (SIZE(Vde,2) .NE. SIZE(QOper,2))) THEN
        WRITE(outUnit,*) "The number of decay operators should", &
             " equal to the number of decay states. STOP"
        WRITE(outUnit,*) SIZE(Vde,2) , SIZE(QOper,2)
        STOP
     ENDIF
     IF(ASSOCIATED(DOper) .AND. (SIZE(Vde,2) .NE. SIZE(DOper,2))) THEN
        WRITE (outUnit,*) "The number of dipole moment operators should", &
             " equal to the number of decay states. STOP"
        STOP
     ENDIF
  ELSEIF (ASSOCIATED(QOper) .AND. (SIZE(Vfi,2) .NE. SIZE(QOper,2))) THEN
     WRITE(outUnit,*) "The number of dipole moment operators should", &
             " equal to the number of final states. STOP"
     WRITE(outUnit,*) SIZE(Vfi,2) , SIZE(QOper,2)
     STOP
     ! Check if splitting between the potentials is defined.
     ! If not set it to 0.0.
  ELSE
     IF(ASSOCIATED(Vspl))  THEN
        Num_spl = SIZE(Vspl,2)
     ELSE
        Num_spl = 0
     ENDIF     ! find the number of the splittings defined so far
     count = 0 ! counter for the splittings between the potentials
     DO i=1, SIZE(Vde,2)
        DO j=i+1, SIZE(Vde,2)
           count = count + 1
           IF (count <= Num_spl) CYCLE
           ALLOCATE(tmp(N)); tmp = 0.0
           CALL IncrDimension(Vspl, tmp); DEALLOCATE(tmp)
        ENDDO
     ENDDO
  ENDIF
  
  
  IF(.NOT. areReqMetP(augEne) ) RETURN

  ! POST-INITIALIZATION PART
  ! ========================
  t(1:1) = MINVAL(Vgr)
  Vgr = Vgr - t(1)
  IF(ASSOCIATED(Vde)) Vde = Vde - t(1)
  IF(ASSOCIATED(Vfi)) Vfi = Vfi - t(1)

  ALLOCATE(Fi(N))

  IF(FinalGamma>0) THEN
     CextFact = MAX( 1, &
          &INT( (8.5/REAL(2**nextpow2(8*EptNum))) / FinalGamma ) )
     WRITE(0,*) "C extension factor set to ", CextFact
  END IF
  

  IF(augEne /=0.0) CALL initQOper(qOrder, qType, qAlpha, momentum, &
       & SWITCH_rho, SWITCH_rc, SWITCH_deltar, dx)
  ! Find the initial state
  ! =====================================
  IF(.NOT. FindGroundFull(Vgr, mass, dx, E0, Fi) ) RETURN

  If(dumpGround.AND.masterP)&
       & CALL dumpFunction(Fi, minX, dx, TRIM(ofName)//".fGround")

  IF(dumpPot.AND.masterP) THEN
     t(1) = 1.0
     t(2) = 0.0
     IF(optEvPotentials) t(1) = au2ev
     IF(optE0Shifted) t(2) = E0
     CALL dumpFunction((Vgr-t(2))*t(1), minX, dx, TRIM(ofName)//".init")
     IF(ASSOCIATED(Vde)) THEN
        DO i=1,SIZE(Vde,2)
           WRITE(WORD,'(I6)') i
           CALL dumpFunction((Vde(:,i)-t(2))*t(1), minX, dx, &
                &TRIM(ofName)//".decay"//TRIM(ADJUSTL(WORD)))
        ENDDO
     ENDIF

     IF(ASSOCIATED(Vfi)) THEN
        DO i=1,SIZE(Vfi,2)
           WRITE(WORD,'(I6)') i
           CALL dumpFunction((Vfi(:,i)-t(2))*t(1), minX, dx, &
                &TRIM(ofName)//".final"//TRIM(ADJUSTL(WORD)))
        ENDDO
     ENDIF
  END IF
  
  IF (verbLvl>0.AND.masterP) THEN
     xAvg = SUM( (/ (Fi(i+1)**2*(MinX+i*dx), i=0,N-1) /) )
     pos=MAXLOC(ABS(Fi))
     WRITE(0,"('Initial E=',F11.5,' (eV) Packet at ',F7.3,"//&
          &"' au')")E0*au2ev, xAvg
!     CALL PrintFVBra(Fi, Vgr, minX, dx) !vkr-->
!     IF(ASSOCIATED(Vde)) THEN
!        DO i=1,SIZE(Vde,2)
!           CALL PrintFVBra(Fi, Vde(:,i), minX, dx)
!        ENDDO
!     ENDIF !<--vkr
  END IF
  IF(masterP.AND.dumpCoreE0) THEN
     DO i=1,SIZE(Vde,2)
        IF(.NOT. FindGroundFull(Vde(:,i), mass, dx, eVec=t)) RETURN
        WRITE(0,"(A18,'lowest vibr. states in core PES, E-E0 (eV)',/A19,6F10.4)")&
             & ' ', ' ',(t-e0)*au2ev
     ENDDO
  END IF
  IF(masterP.AND.dumpFinalE0) THEN
     DO i=1,SIZE(Vfi,2)
        IF(.NOT. FindGroundFull(Vfi(:,i), mass, dx, eVec=t)) RETURN
        WRITE(0,"(A18,'lowest states in final PES, E-E0 (eV):',/,A19,10F10.4)")&
             & ' ', ' ',(t-e0)*au2ev
     ENDDO
  END IF
  IF(verbLvl>1 .AND. ASSOCIATED(QOper)) THEN
     DO i=1,SIZE(QOper,2)
        CALL PrintFVBra(REAL(QOper(:,i)), AIMAG(QOper(:,i)), minX, dx)
     ENDDO
  ENDIF
  IF(ASSOCIATED(QOper).AND.dumpQoper) THEN
     DO i=1,SIZE(QOper,2)
        WRITE(WORD,'(I6)') i
        CALL dumpFunction(REAL(QOper(:,i)), minX, dx, TRIM(ofName)//".rQoper"&
             &//TRIM(ADJUSTL(WORD)))
     ENDDO
  ENDIF
  IF(ASSOCIATED(DOper).AND.dumpDoper) THEN
     DO i=1,SIZE(DOper,2)
        WRITE(WORD,'(I6)') i
        CALL dumpFunction(DOper(:,i), minX, dx, TRIM(ofName)//".rDoper"&
             &//TRIM(ADJUSTL(WORD)))
     ENDDO
  ENDIF
  ! Dump the splittings.
  IF(ASSOCIATED(Vspl).AND.dumpSplit) THEN
     DO i=1,SIZE(Vspl,2)
        WRITE(WORD,'(I6)') i
        CALL dumpFunction(Vspl(:,i), minX, dx, TRIM(ofName)//".rSplit"&
             &//TRIM(ADJUSTL(WORD)))
     ENDDO
  ENDIF
  ! WATCH OUT - you'll get wrong energy of the ground state if you substract!
  !ep=MINVAL(Vgr); Vgr=Vgr-ep; Vde=Vde-ep; Vfi=Vfi-ep
  ! the code below is much better.
  IF(masterP) WRITE(outUnit,*)&
       &'---------------------------------------------------------------------'
  initdata = .TRUE.
  RETURN
101 IF(masterP) WRITE(outUnit,*) &
         &'Parameter reading error. Last keyword: ', WORD
  RETURN
END FUNCTION initdata


  ! ===================================================================
  ! checks how wide is the potential well at given energy
  SUBROUTINE checkPotWidth(V, om)
    REAL(mk), INTENT(IN) :: V(:)
    REAL(mk), INTENT(IN) ::  om

    REAL(mk) :: dx
    INTEGER  :: mn(1), loi(1), hii(1)

    dx = (MaxX-MinX)/(N-1)
    mn = minloc(V)
    IF(V(mn(1))<om) THEN
       loi = MAXLOC(V(1:mn(1)), V(1:mn(1))<om)
       hii = MAXLOC(V(mn(1):SIZE(V)), V(mn(1):SIZE(V))<om)+mn(1)-1
       WRITE(0,"('The well width at the energy ',F8.3,' eV is ',F6.3,"//&
            &"' au (',F9.2,',',F9.2,')')") om*au2ev, (hii(1)-loi(1))*dx,&
            &(loi(1)-1)*dx+MinX, (hii(1)-1)*dx+MinX
    ELSE
       WRITE(0,"('Tuning to level ',F8.3,' eV below the well bottom at ',"//&
            &"F8.3,' eV')") om*au2ev, V(mn(1))*au2ev
    END IF
  END SUBROUTINE checkPotWidth

  ! dumpPotential ---------------------------------------------------
  SUBROUTINE dumpPotential(V, fName)
    IMPLICIT NONE
    REAL(mk), INTENT(IN)     :: V(:)
    CHARACTER(*), INTENT(IN) :: fName
    !REAL(mk) :: dx; INTEGER  :: i

    CALL dumpFunction(V,minX, (maxX-minX)/(N-1), fName)
    ! OPEN(20, FILE=fName)
    ! WRITE(20, "(2(E14.6))") ( (/ minX+i*dx, V(i+1) /), i=0,SIZE(V)-1)
    ! CLOSE(20)
  END SUBROUTINE dumpPotential

  SUBROUTINE dumpFunction(V, x0, dx, fName)
    IMPLICIT NONE
    REAL(mk), INTENT(IN)     :: V(:), x0, dx
    CHARACTER(*), INTENT(IN) :: fName

    INTEGER  :: i

    OPEN(20, FILE=fName)
    WRITE(20, "(2(E14.6))") ( (/ x0+i*dx, V(i+1) /), i=0,SIZE(V)-1)
    CLOSE(20)
  END SUBROUTINE dumpFunction
  
  SUBROUTINE dumpFunctionC(f, x0, dx, fName)
    IMPLICIT NONE
    COMPLEX(mk), INTENT(IN)  :: f(:)
    REAL(mk), INTENT(IN)     :: x0, dx
    CHARACTER(*), INTENT(IN) :: fName

    INTEGER  :: i

    OPEN(20, FILE=fName)
    WRITE(20, "(3(E14.6))") ( (/ x0+i*dx, REAL(f(i+1)), AIMAG(f(i+1)) /), &
         &i=0,SIZE(f)-1)
    CLOSE(20)
  END SUBROUTINE dumpFunctionC

  ! loadLine -------------------------------------------------------
  ! loads data of linear potential and sets given vector vVec
  SUBROUTINE loadLine(iUnit, name, vVec, units)
    IMPLICIT NONE
    INTEGER, INTENT(IN)                :: iUnit
    CHARACTER(5), INTENT(IN)           :: name
    character(2), intent(IN)           :: units
    REAL(mk), DIMENSION(:),POINTER     :: vVec
    REAL(mk)                :: x0, e0, slope, dx
    INTEGER                 :: i

    IF (ASSOCIATED(vVec)) CALL twicePot(name)
    IF (N<=0 .OR. minX>=maxX) CALL potDef(name,'N    ')
    IF (MASS<0) CALL potDef(name,'mass ')
    ALLOCATE(vVec(N))
    READ(iUnit,*, ERR=101,END=101) x0, e0, slope
    IF(masterP) WRITE (outUnit,123) name, x0, e0, slope
    if(units.eq."EV") then
        e0=e0*ev2au; slope=slope*ev2au
    endif
    dx=(maxX-minX)/(N-1)
    vVec=(/ (slope*(I*dx+minX-x0)+e0, I=0, N-1) /)
    RETURN
101 IF(masterP) WRITE(outUnit,*) name,': parameter reading error.'
  STOP
123 FORMAT(A,': Value at r=',F7.3,' au: ',F9.4,' au. Grows ',F8.4,'au/au')
  END SUBROUTINE loadLine

!vk add gaussian
  ! loadLine -------------------------------------------------------
  ! loads data of linear potential and sets given vector vVec
  SUBROUTINE loadGauss(iUnit, name, vVec, units)
    IMPLICIT NONE
    INTEGER, INTENT(IN)                :: iUnit
    CHARACTER(5), INTENT(IN)           :: name
    character(2), intent(IN)           :: units
    REAL(mk), DIMENSION(:),POINTER     :: vVec
    REAL(mk)                :: x0, e0, hwhm, dx
    INTEGER                 :: i

    IF (ASSOCIATED(vVec)) CALL twicePot(name)
    IF (N<=0 .OR. minX>=maxX) CALL potDef(name,'N    ')
    IF (MASS<0) CALL potDef(name,'mass ')
    ALLOCATE(vVec(N))
    READ(iUnit,*, ERR=101,END=101) x0, e0, hwhm
    IF(masterP) WRITE (outUnit,123) name, x0, e0, hwhm
    dx=(maxX-minX)/(N-1)
!    vVec=(/ (slope*(I*dx+minX-x0)+e0, I=0, N-1) /)
		vVec=(/ (epsilon(e0)+e0*exp(-log(2.0)*((I*dx+minX-x0)/hwhm)**2), I=0, N-1) /)
    RETURN
101 IF(masterP) WRITE(outUnit,*) name,': parameter reading error.'
  STOP
123 FORMAT(A,': Value at r=',F7.3,' au: ',F9.4,' au. Width (HWHM) ',F8.4,'1/au')
  END SUBROUTINE loadGauss

  
  ! loadHarm -------------------------------------------------------
  ! loads data of a harmonic potential and sets given vector vVec
  SUBROUTINE loadHarm(iUnit, name, vVec)
    IMPLICIT NONE
    CHARACTER(5), INTENT(IN)           :: name
    REAL(mk), DIMENSION(:),POINTER     :: vVec
    INTEGER, INTENT(IN)                :: iUnit
    REAL(mk)                :: om, rp, ep, dx
    INTEGER                 :: i

    IF (ASSOCIATED(vVec)) CALL twicePot(name)
    IF (N<=0 .OR. minX>=maxX) CALL potDef(name,'N    ')
    IF (MASS<0) CALL potDef(name,'mass ')
    ALLOCATE(vVec(N))
    READ(iUnit,*, ERR=101,END=101) om, rp, ep
    IF(masterP) WRITE (outUnit,123) name,  om, rp, ep
    om=om*ev2au; ep=ep*ev2au; dx=(maxX-minX)/(N-1)
    !WRITE(0,*) rp, minX,(N-1)*dx
    vVec=0.5*mass*(om**2)* (/((I*dx-rp+minX)**2,I=0,N-1)/)+ep
    !WRITE(0,*) vVec(1), vVec(N)
    RETURN
101 IF(masterP) WRITE(outUnit,*) name,': parameter reading error.'
  STOP
123 FORMAT(A,': om=',E10.3,' eV, r0=',E10.3,' au, e0=',E10.3,' eV')
  END SUBROUTINE loadHarm

  ! loadMors -------------------------------------------------------
  SUBROUTINE loadMors(iUnit, name, vVec)
    IMPLICIT NONE
    CHARACTER(5), INTENT(IN)       :: name
    REAL(mk), DIMENSION(:),POINTER :: vVec
    INTEGER, INTENT(IN)            :: iUnit
    REAL(mk)                :: D0,alp, rp, e0, dx
    INTEGER                 :: i

    IF (ASSOCIATED(vVec)) CALL twicePot(name)
    IF (N<=0 .OR. minX>=maxX) CALL potDef(name,'N    ')
    IF (MASS<0) CALL potDef(name,'mass ')
    ALLOCATE(vVec(N))
    READ(iUnit,*, ERR=101,END=101) D0, alp, rp, e0
    IF(masterP) WRITE (outUnit,123) name, D0,alp,rp,e0
    D0=D0*ev2au; e0=e0*ev2au; dx=(maxX-minX)/(N-1)
    vVec=D0*(1-EXP(-alp*( (/ (i*dx+minX-rp, i=0,N-1) /))))**2+e0     !-D0
    RETURN
101 IF(masterP) WRITE(outUnit,*) name,': parameter reading error.'
    STOP
123 FORMAT(A,' D0:',E10.3,' eV, alpha:',E10.3,', r0:',E10.3,' au, e0',E10.3)
  END SUBROUTINE loadMors
  
  ! loadMorseOmar -------------------------------------------------------
  ! loads data of a Morse potential with format: omega (eV), anharmonicity
  ! constant (eV), r0 (AA), energy of the bottom of the potential (eV)
  SUBROUTINE loadMorseOmar(iUnit, name, vVec)
    IMPLICIT NONE
    CHARACTER(5), INTENT(IN)       :: name
    REAL(mk), DIMENSION(:),POINTER :: vVec
    INTEGER, INTENT(IN)            :: iUnit
    REAL(mk)                :: D0,alp, rp, e0, dx, anh, om
    INTEGER                 :: i

    IF (ASSOCIATED(vVec)) CALL twicePot(name)
    IF (N<=0 .OR. minX>=maxX) CALL potDef(name,'N    ')
    IF (MASS<0) CALL potDef(name,'mass ')
    ALLOCATE(vVec(N))
    READ(iUnit,*, ERR=101,END=101) om, anh, rp, e0
    D0  = om**2/(4*anh)*ev2au
    alp = SQRT(2*mass*anh*ev2au)
    rp  = rp/0.52917725
    e0  = e0*ev2au
    dx=(maxX-minX)/(N-1)
    IF(masterP) WRITE (outUnit,123) name,D0*au2ev,alp,rp,e0*au2ev
    vVec=D0*(1-EXP(-alp*( (/ (i*dx+minX-rp, i=0,N-1) /))))**2+e0     !-D0
    RETURN
101 IF(masterP) WRITE(outUnit,*) name,': parameter reading error.'
    STOP
123 FORMAT(A,' D0:',E10.3,' eV, alpha:',E10.3,', r0:',E10.3,' au, e0',E10.3)
  END SUBROUTINE loadMorseOmar
  
  ! loadMorseCm -------------------------------------------------------
  ! loads data of a Morse potential with format: omega (cm-1), anharmonicity
  ! constant (cm-1), r0 (AA), energy of the bottom of the potential (eV)
  SUBROUTINE loadMorseCm(iUnit, name, vVec)
    IMPLICIT NONE
    CHARACTER(5), INTENT(IN)       :: name
    REAL(mk), DIMENSION(:),POINTER :: vVec
    INTEGER, INTENT(IN)            :: iUnit
    REAL(mk)                :: D0,alp, rp, e0, dx, anh, om
    INTEGER                 :: i

    IF (ASSOCIATED(vVec)) CALL twicePot(name)
    IF (N<=0 .OR. minX>=maxX) CALL potDef(name,'N    ')
    IF (MASS<0) CALL potDef(name,'mass ')
    ALLOCATE(vVec(N))
    READ(iUnit,*, ERR=101,END=101) om, anh, rp, e0
    D0  = om**2/(4*anh)*cm2au
    alp = SQRT(2*mass*anh*cm2au)
    rp  = rp/0.52917725
    e0  = e0*ev2au
    dx=(maxX-minX)/(N-1)
    IF(masterP) WRITE (outUnit,123) name,D0*au2ev,alp,rp,e0*au2ev
    vVec=D0*(1-EXP(-alp*( (/ (i*dx+minX-rp, i=0,N-1) /))))**2+e0     !-D0
    RETURN
101 IF(masterP) WRITE(outUnit,*) name,': parameter reading error.'
    STOP
123 FORMAT(A,' D0:',E10.3,' eV, alpha:',E10.3,', r0:',E10.3,' au, e0',E10.3)
  END SUBROUTINE loadMorseCm

  ! loadHulburtHirschfelder ----------------------------------------------
  ! loads data of a Hulburt-Hirschfelder potential, as given in Herzberg 
  ! book.
  SUBROUTINE loadHulburtHirschfelder(iUnit, name, vVec)
    IMPLICIT NONE
    CHARACTER(5), INTENT(IN)       :: name
    REAL(mk), DIMENSION(:),POINTER :: vVec
    INTEGER, INTENT(IN)            :: iUnit
    REAL(mk)                :: De, ome, omexe, r0, be, alphe, e0
    REAL(mk)                :: dx, cc, cb, cbeta
    INTEGER                 :: i

    IF (ASSOCIATED(vVec)) CALL twicePot(name)
    IF (N<=0 .OR. minX>=maxX) CALL potDef(name,'N    ')
    IF (MASS<0) CALL potDef(name,'mass ')
    ALLOCATE(vVec(N))
    READ(iUnit,*, ERR=101,END=101) De, ome, omexe, r0, be, alphe, e0

    ! I don't like the mass conversion constant in the equation below...
    cbeta = 0.06442*ome * SQRT(mass/(De*1836.15270))
    cc = 1.0 - (1.0_mk + (alphe * ome)/(6.0* be**2))/(cbeta*r0)
    cb = 1.0 + (7.0/12.0- &
         & (5.0/4.0 +5.0*alphe*ome/(12.0*be**2) &
         &+ 5*(alphe*ome)**2/(144*be**4)&
         &-2.0*omexe/(3.0*be))/(cbeta*r0)**2)/cc

    dx=(maxX-minX)/(N-1)
    IF(masterP) WRITE (outUnit,123) name, De*au2ev, r0, e0
    !WRITE(outUnit, *) "3.cBeta:", cbeta, " cc: ", cc," cb:", cb
    !WRITE(outUnit, *) "3.cBeta2", ome*cm2au * SQRT(mass/(De*cm2au))


    vVec = (/ (i*dx+minX-r0, i=0,N-1) /) ! the R-R0 vector...
    vVec=De*cm2au*(&
         & (1-EXP(-cbeta*vVec))**2 &
         &+ cc*(cbeta*vVec)**3 * EXP(-2.0*cbeta*vVec)*(1+cb*cbeta*vVec)) &
         &+e0*ev2au
    RETURN
101 IF(masterP) WRITE(outUnit,*) name,': parameter reading error.'
    STOP
123 FORMAT(A,' D0:',E10.3,' eV, r0:',E10.3,' au, e0',E10.3)
  END SUBROUTINE loadHulburtHirschfelder

  ! loadBuck -------------------------------------------------------
! vk -- introduce beta
  SUBROUTINE loadBuck(iUnit, name, Vvec)
    IMPLICIT NONE
    CHARACTER(5), INTENT(IN)        :: name
    REAL(mk), DIMENSION(:), POINTER :: vVec
    INTEGER, INTENT(IN)             :: iUnit
    REAL(mk)                :: A,alp,B,e0,dx,beta
    INTEGER                 :: i

    IF (ASSOCIATED(vVec)) CALL twicePot(name)
    IF (N<=0 .OR. minX>maxX) CALL potDef(name,'N    ')
    ALLOCATE(vVec(N))
    READ(iUnit,*, ERR=101,END=101) A, alp, B, beta, e0
    IF(masterP) WRITE (outUnit,123) name,": A[eV], alpha[1/a.u], B , beta [a.u.], e0[eV]",&
         &A, alp, B, beta, e0
    A=A*ev2au; B=B*ev2au; e0=e0*ev2au
    dx=(maxX-minX)/(N-1)
    vVec=A*EXP(-alp*( (/ (i*dx+minX, i=0,N-1) /)))+e0&
         &-B/(/ ((i*dx+minX)**beta, i=0,N-1) /)
    RETURN
101 IF(masterP) WRITE(outUnit,*) name,': parameter reading error.'
  STOP
123 FORMAT(A,A,5E10.3)
!123 FORMAT(A,A,E10.3,E10.3,E10.3,E10.3)
  END SUBROUTINE loadBuck
  
  ! error handling procedures ---------------------------------------
  SUBROUTINE potDef(pott, var)
    CHARACTER(5), INTENT(IN) :: pott, var
    IF(masterP) WRITE (outUnit,*) &
         &pott,' encountered but ',var,' not defined. STOP'
    STOP
  END SUBROUTINE potDef
  
  SUBROUTINE twicePot(pott)
    CHARACTER(LEN=5), INTENT(IN) :: pott
    IF(masterP) WRITE (outUnit,*) pott,' defined twice. STOP'
    STOP
  END SUBROUTINE twicePot

  ! processOptions --------------------------------------------------
  ! processes different options to RAM.
  SUBROUTINE processOptions(iUnit)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: iUnit
    CHARACTER(256) :: optionsLine
    INTEGER        :: idx

    READ(iUnit,"(A)",ERR=101,END=101) optionsLine
    DO
       optionsLine = ADJUSTL(optionsLine)
       idx = INDEX(optionsLine,' ')
       if(idx == 1) EXIT
       SELECT CASE(optionsLine(1:idx-1))
       CASE("CORE_DUMPED")
          optPrintCoreDumped = .TRUE.
          IF(masterP) WRITE(outUnit,"(A)", ADVANCE="NO") "core evol. dumped, "
       CASE("E0SHIFTED")
          optE0Shifted = .TRUE.
          IF(masterP) WRITE(outUnit,"(A)", ADVANCE="NO") "E0_shifted "
       CASE("EVPOTS")
          optEvPotentials = .TRUE.
          IF(masterP) WRITE(outUnit,"(A)", ADVANCE="NO") "eV_Potentials "
       CASE("NORMALIZE")
          optNormalize = .TRUE.
          IF(masterP) WRITE(outUnit,"(A)", ADVANCE="NO") "Normalize "
       CASE("TEST_DISS")
          optTestDissociation = .TRUE.
          IF(masterP) WRITE(outUnit,"(A)", ADVANCE="NO") "TEST_DISS "
       CASE DEFAULT
          IF(masterP) WRITE (outUnit,*)"processOptions: unknown option ", &
               &optionsLine(1:idx-1),&
               &" known are: E0SHIFTED EVPOTS NORMALIZE TEST_DISS"
       END SELECT
       optionsLine = optionsLine(idx:LEN(optionsLine))
    ENDDO
    IF(masterP) WRITE(outUnit,*)": option(s) selected."
    RETURN
101 IF(masterP) WRITE(outUnit,*)'processOptions: line with options expected.'
    STOP
  END SUBROUTINE processOptions
  
  ! processDumps ----------------------------------------------------
  ! loads subsequent line and determines what intermediate data 
  ! shall be dumped to a file.
  ! SETS: dump* variables, used in other parts of the program
  SUBROUTINE processDumps(iUnit)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: iUnit
    CHARACTER(256)  :: dumpLine
    INTEGER         :: idx
    READ (iUnit,"(A)",ERR=101,END=101) dumpLine
    DO
       dumpLine=ADJUSTL(dumpLine)
       IF (dumpLine(1:1)==" ") EXIT
       idx=INDEX(dumpLine," ")
       SELECT CASE(dumpLine(1:idx-1))
       CASE("AVERAGES")
          IF(masterP) THEN 
             dumpAvgOps=.TRUE.; WRITE (outUnit,9,ADVANCE="NO") "r&p Averages "
          END IF
       CASE("DOMEGA")
          IF(masterP) THEN
             dumpDOm=.TRUE.; WRITE (outUnit,9,ADVANCE="NO") "dumpDOm " 
          END IF
       CASE("GROUND")
          IF(masterP) THEN
             dumpGround=.TRUE.; WRITE (outUnit,9,ADVANCE="NO") "dumpGround " 
          END IF
       CASE("FFT")
          IF(masterP) THEN
             dumpFFT=.TRUE.; WRITE (outUnit,9,ADVANCE="NO") "dumpFFT "
          END IF
       CASE("DURATION")
          dumpDur=.TRUE.; IF(masterP) WRITE(outUnit,9,ADVANCE="NO")"Duration "
       CASE("CORE_E0")
          dumpCoreE0=.TRUE.;IF(masterP) WRITE(outUnit,9,ADVANCE="NO")"coreE0 "
       CASE("FINAL_E0")
          dumpFinalE0=.TRUE.;IF(masterP) WRITE(outUnit,9,ADVANCE="NO")"finalE0 "
       CASE("POTENTIALS")
          IF(masterP) THEN 
             dumpPot=.TRUE.; WRITE (outUnit,9,ADVANCE="NO") "dumpPot "
          END IF
       CASE("QOPER")
          dumpQoper=.TRUE.;IF(masterP) WRITE(outUnit,9,ADVANCE="NO")"dump Qr "
       CASE("SPLIT")
          dumpSplit=.TRUE.;IF(masterP) WRITE(outUnit,9,ADVANCE="NO")"dump Vr "
       CASE("DIPOL")
          dumpDoper=.TRUE.;IF(masterP) WRITE(outUnit,9,ADVANCE="NO")"dump Dr "
       CASE("SIGMA_T")
          IF(masterP) THEN 
             dumpSigmaT=.TRUE.; WRITE (outUnit,9,ADVANCE="NO") "sigma(t) "
          END IF

       CASE DEFAULT
          IF(masterP) WRITE (outUnit,"('processDumps: unknown dump ',A,/"//&
               &"'Known are: AVERAGES CROSSCOMP POTENTIALS GROUND FFT "//&
               &"DURATION CORE_E0 FINAL_E0')") dumpLine(1:idx-1)
       END SELECT
       dumpLine=dumpLine(idx:LEN(dumpLine))
    END DO
    IF(masterP) WRITE(outUnit,*)": selected."
    RETURN
101 IF(masterP) WRITE(outUnit,*)'processDumps: line with parameters expected.'
    STOP
9   FORMAT (A)
  END SUBROUTINE processDumps
      
  ! loadSpline -----------------------------------------------------
  ! loads set of points from stdio and then fits a spline to it.
  ! uses netlib/fmm routines
  SUBROUTINE loadSpline(iUnit,nm,vVec)
    IMPLICIT NONE
    INTEGER, INTENT(IN)      :: iUnit
    CHARACTER(5), INTENT(IN) :: nm
    REAL(mk), POINTER        :: vVec(:)

    INTEGER                             :: ptCnt, i
    REAL(mk), ALLOCATABLE               :: pts(:,:)
    REAL(mk), DIMENSION(:), ALLOCATABLE :: ac, bc, cc
    REAL(mk)                            :: dx
#ifndef SYS_CRAY
    DOUBLE PRECISION seval
#else
    REAL seval
#endif
    IF (ASSOCIATED(vVec)) CALL twicePot(nm)
    IF (N<=0 .OR. MinX>MaxX) CALL potDef(nm,'N    ')
    ALLOCATE(vVec(N))
    READ(iUnit,*, ERR=101,END=101) ptCnt
    IF (masterP) WRITE (outUnit,*) nm," defined on ", ptCnt," points."
    ALLOCATE(pts(ptCnt,2), ac(ptCnt), bc(ptCnt),cc(ptCnt))
    READ(iUnit,*, ERR=101,END=101) (pts(I,:),I=1,ptCnt)
    !WRITE(0,'(2(F17.7))') (pts(I,:),I=1,ptCnt)
    CALL spline(ptCnt,pts(:,1),pts(:,2),ac,bc,cc)
    dx=(maxX-minX)/(N-1)

    DO i=1, N
       vVec(i)=seval(ptCnt,minX+(i-1)*dx,pts(:,1),pts(:,2),ac,bc,cc)
    END DO
    ! DEALLOCATE(pts,ac,bc,cc) - not necessary for allocatable
    RETURN
101 IF(masterP) WRITE (outUnit,*) 'loadSpline: failed for ',nm
    ! DEALLOCATE(pts,ac,bc,cc) - not necessary for allocatable
  END SUBROUTINE loadSpline

  ! IncrDimension -----------------------------------------------------
  ! For multidimensional input arrays like Vde, QOper, DOper
  ! Adds a vector at the and of a matrix. Deallocate "add" outside.
  SUBROUTINE IncrDimension(array,add)
    IMPLICIT NONE
    REAL(mk), POINTER :: array(:,:), add(:)
    
    REAL(mk), POINTER :: tmp(:,:)    

    IF (.NOT. ASSOCIATED(array)) THEN 
       ALLOCATE(array(SIZE(add),1))
       array(:,1) = add;
       RETURN
    END IF
    NULLIFY(tmp)
    ALLOCATE(tmp(SIZE(array,1),SIZE(array,2) + 1))
    tmp(:,1:SIZE(array,2)) = array
    tmp(:,SIZE(array,2) + 1) = add
    DEALLOCATE(array)
    array => tmp
    RETURN
  END SUBROUTINE IncrDimension


  ! OUTPUT ROUTINES ==================================================
  ! ==================================================================
  ! PrintFVBra -------------------------------------------------------
  ! must use Max(maxval) to handle zero potential properly.
  SUBROUTINE PrintFVBra(F,V,minX,dx)
    implicit none
    REAL(mk), dimension(:) :: F,V
    REAL(mk)  :: minX, dx

    REAL(mk) a, b, c
    integer :: i, stp

    a=MINval(V(:)); b=MAX(0.01_mk,MAXval(V)-a)
    c=MAXval(abs(F))**2
    stp=MAX(1,SIZE(F)/1024)
  
 !!!VKA 
    do i=1, size(F),stp
      write (*,'(3E14.5)') (((i-1)*dx)+minX), (V(i)-a)*c/b,F(i)**2 
    enddo
 !!!VKA
    
 !   print *, '(plot 2 '
 !   do i=1, size(F),stp
 !      write (*,'(3E12.4)') (((i-1)*dx)+minX), (V(i)-a)*c/b,F(i)**2 
 !   enddo
 !   print *,') (Wait 1000)'
  end SUBROUTINE PrintFVBra

!   ! OUTPUT ROUTINES ==================================================
!   ! ==================================================================
!   ! PrintFVBra -------------------------------------------------------
!   ! must use Max(maxval) to handle zero potential properly.
!   SUBROUTINE PrintFVBra(F,V,minX,dx)
!     implicit none
!     REAL(mk), DIMENSION(:) :: F
!     REAL(mk), DIMENSION(:,:) :: V
!     REAL(mk)  :: minX, dx

!     REAL(mk) a(SIZE(F,2)), b(SIZE(F,2)), c(SIZE(F,2))
!     INTEGER :: i(SIZE(F,2)), stp(SIZE(F,2)), k
! !    CALL PrintFVBra(ABS(Fact(low:high))**2,V(:,i), vSt, dx)
!     PRINT *, '(plot 2 '
!     DO k=1,SIZE(F,2) !loop over states
!        a(k)=MINVAL(V(:,k)); b(k)=MAX(0.01_mk,MAXVAL(V(:,k))-a(k))
!        c(k)=MAXVAL(ABS(F(:,k)))**2
!        stp(k)=MAX(1,SIZE(F:,k)/512)
       
!        DO i(k)=1, SIZE(F,k),stp(k)
! !          WRITE (*,'(E12.4)') &
!           WRITE (*,321) &
!                &(((i(k)-1)*dx)+minX), (V(i,:)-a)*c/b,F(i,:)**2 
!        ENDDO
!     ENDDO
!     PRINT *,') (Wait 1000)'
! 321 FORMAT()
!   end SUBROUTINE PrintFVBra

  ! saveDuration -----------------------------------------------------
  SUBROUTINE saveDuration(maxJ, deltaOm, durAr)
    IMPLICIT NONE
    INTEGER,INTENT(IN)   :: maxJ
    REAL(mk), INTENT(IN) :: deltaOm, durAr(:)
    INTEGER  i
    
    OPEN(UNIT=9,STATUS='REPLACE',ACTION='WRITE',FILE=TRIM(ofName)//".duration")
    WRITE(9,*) "omega (eV) duration (a.u)"
    WRITE(9,"((F9.3,E13.4))") ( &
         & (/(gPos-gWidth+gStep*deltaOm*i)/ev2au,durAr(i+1)/),&
         & i=0,maxJ-1)
    CLOSE(UNIT=9)
  END SUBROUTINE saveDuration

  ! saveSigma --------------------------------------------------------
  ! saves sigma Array - if it is allocated.
  ! chooses right shift that (Emin+Emax) is roughly in the middle. Uses 
  ! the periodicity of C(\omega)
  SUBROUTINE saveSigma(maxJ, deltaOm, Sigma, eneShift)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: maxJ
    REAL(mk), INTENT(IN) :: deltaOm, Sigma(:,:)
    REAL(mk), OPTIONAL   :: eneShift

    REAL           :: om1
    INTEGER        :: i, Clen, meIdx, shift
    CHARACTER(40)  :: myFmt

    CLen = SIZE(Sigma,1)
    ! find suitable shifts
    !---------------------
    meIdx = INT((EMin-gPos+gWidth)/(2*deltaOm))
    shift = Clen/2-INT((0.5*(EMax+EMin)-gPos+gWidth)/deltaOm)
    om1 = gPos-gWidth-shift*deltaOm
    IF(PRESENT(eneShift)) om1 = om1+eneShift
    shift = MOD( MOD(-shift, CLen)+CLen, CLen)

    ! open and write according to the new algorithm...
    ! ------------------------------------------------
    i=MAX(80,(maxJ+2)*12)
    OPEN(9,STATUS='REPLACE',ACTION='WRITE',RECL=i,FILE=TRIM(ofName)//".spec")
    WRITE(9,FMT="(I8,"//i2s(maxJ)//"(F12.4))") maxJ,&
         & ((gPos-gWidth+gStep*deltaOm*i)/ev2au,I=0,maxJ-1)
    
    myFmt="(F8.3,"//i2s(maxJ)//"(E12.4))"
    DO i=1, SIZE(Sigma,1)
       WRITE(9,FMT=myFmt) (om1+(i-1)*deltaOm)/ev2au,&
       & Sigma(MOD(i+shift-1,Clen)+1,1:maxJ)
    END DO
    CLOSE(UNIT=9)

  END SUBROUTINE saveSigma


  ! printCArr --------------------------------------------------------
  ! The frequency axis IS correct.
  ! however the code is to BE IMPROVED! The handling of wrapped plots is 
  ! uncertain.
  ! Y - the cross-section array, where the first element corresponds to
  ! freq. om0 and the spacing between elements is deltaOm
  ! only part between frLo and frHi is to be printed.
  SUBROUTINE printCArr(Y, om0, frLo, frHi, deltaOm, omInc)
    IMPLICIT NONE
    REAL(mk), INTENT(IN) :: Y(:)
    REAL(mk), INTENT(IN) :: om0, frLo, frHi, deltaOm, omInc

    INTEGER  :: loIdx,hiIdx, cLen, i
  
    cLen = SIZE(Y)
    loIdx = 1+MAX(0.0_mk,(frLo-om0)/deltaOm)
    hiIdx = MIN(REAL(cLen,mk),(frHi-om0)/deltaOm+1.0)
    
    ! the printing part, without previous splining ----------------
    PRINT *,'(Plot 1'
    DO i=loIdx,hiIdx
       PRINT '(2F15.8)', ((i-1)*deltaOm+om0)/ev2au, Y(i)
    ENDDO
    PRINT *,')(Title The spectra for omInc=',omInc/ev2au,'eV)(Wait 2000)'

  END SUBROUTINE printCArr

  ! CONVERSION ROUTINES
  ! ---------------------------------------------------------------
  ! i2s converts given integer number to a string.
  FUNCTION i2s(i)
    INTEGER, INTENT(IN) :: i
    CHARACTER(7)        :: i2s
    INTEGER :: j,k

    j=i
    k = 7
    IF (j<0) THEN
       WRITE (0,*)'Cannot convert negative integers. STOP'
       STOP
    END IF
    IF(j==0) THEN
       i2s='0'
    ELSE
       i2s=''
       DO WHILE (j>0.AND.k>0)
          i2s=ACHAR(MOD(j,10)+48)//i2s
          k=k-1; j=j/10
       END DO
    END IF
  END FUNCTION i2s

  ! ---------------------------------------------------------------
  INTEGER FUNCTION nextpow2(l)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: l
 
  nextpow2=0
  DO WHILE(L > 2**nextpow2)
     nextpow2=nextpow2+1
  END DO
  END FUNCTION nextpow2

END MODULE IOProcessor
