! ram.f90 solves the Raman scattering with algorithm extended to
! non-monochromatic incident radiation. At the moment creates file with 
! the data and the convolution is performed by a postprocessor.
! Notation according to TheoBackup.tex
! Pawel Salek, pawsa@ifm.liu.se
! 980318
! NOTES:
! 990225: do not dump more often than every DumpFreq minutes.
! 990302: Version two started.
! 991004: duration was not saved.
! 200810xx: VC core-excited state implemented by Yasen Velcov
! 201002xx: Raman scattering including core-excited VC
!
PROGRAM RamanTwo
  USE LapackInterface
  USE PFSolver
  USE IOProcessor
  IMPLICIT NONE
  REAL(mk), PARAMETER                 :: DumpFreq = 60.0

  COMPLEX(mk), ALLOCATABLE :: C(:), cA(:), cB(:), cC(:)
  REAL(mk)                 :: dtSamp, omInc, lastDump, omZero
  REAL(mk), ALLOCATABLE    :: durAr(:)
  INTEGER                  :: stp, cnt, crt, Clen, i
  TYPE(spectrum)           :: Sigma, SigDir, SigCross, Sig1, Sig2,Sig12
	
  lastDump=0.0
  IF(.NOT. initdata(iUnit=5)) STOP "Initialization error."

  ! ---------------------------------------------------------------------
  ! do the calculations. The scan in the energy space is done by a ******
  ! Fourier transform.                                             ******
  ! the scan for the incoming frequency range is done as follows:
  ! the points are chosen in such a way that spacing is a integer 
  ! multiplicity of deltaOm. The number of points is limited to make 
  ! sure that the distance between points in omega space is not less than
  ! deltaOm. The calculated date is stored in Sigma, which row j correspond
  ! to incoming frequency gPos-gWidth-(Clength-j)*deltaOm. 
  CALL initSpectrum(Sigma, gPos, gWidth, gStep)
  Clen = SIZE(Sigma%sigma,1)*CextFact
  ALLOCATE(C(Clen),durAr(SIZE(Sigma%sigma,2)))
  SELECT CASE(Mode)
  CASE(ModeRamDir) 
     CALL initSpectrum(SigDir,  gPos, gWidth, gStep)
     CALL initSpectrum(SigCross,gPos, gWidth, gStep)
      ALLOCATE(cA(Clen), cB(Clen))
  CASE(ModeRamSym)
     CALL initSpectrum(Sig1,  gPos, gWidth, gStep)
     CALL initSpectrum(Sig2,  gPos, gWidth, gStep)
     CALL initSpectrum(Sig12, gPos, gWidth, gStep)
     ALLOCATE( cA(Clen), cB(Clen), cC(Clen) )
  END SELECT
  
  dtSamp=2*PI/(Sigma%deltaOm*(SIZE(Sigma%sigma,1)-1))
                    ! the time sampling rate; a.u.


  ! THE MAIN LOOP OVER INCOMING FREQUENCIES ---------------------------
  DO stp=1, SIZE(Sigma%sigma,2)
     omInc = Sigma%zeroInc + (stp-1) * Sigma%deltaOm * Sigma%inStep
     WRITE (0,"(/,'S. no.',I3,'    : omega=',F8.4,' eV')")stp, omInc*au2ev

     IF(dumpCoreE0) THEN
        DO i=1, SIZE(Vde,2)
           CALL checkPotWidth(Vde(:,i), E0+omInc)
        ENDDO
     ENDIF

     SELECT CASE(Mode)
     CASE(ModeAbsorp)
        CALL calcAbsorption(Fi, E0, omInc, dtSamp, C, omZero)
     CASE(ModeDecay)
        CALL calcDecay (Fi, E0, omInc, durAr(stp))
     CASE(ModeRaman)
        CALL calcRaman (Fi, E0, omInc, dtSamp, C, omZero, durAr(stp))
     CASE(ModeRamDir)
        CALL calcRamDir(Fi, E0, omInc, dtSamp, C, cA, cB, omZero,&
             & durAr(stp))
        CALL placeSpectrum(SigDir,   stp, REAL(cA), omZero)
        CALL placeSpectrum(SigCross, stp, REAL(cB), omZero)
     CASE(ModeRamSym)
        CALL calcRamSym(Fi, E0, omInc, dtSamp, C, cA,cB,cC, omZero,&
             & durAr(stp))
        CALL placeSpectrum(Sig1, stp, REAL(cA), omZero)
        CALL placeSpectrum(Sig2, stp, REAL(cB), omZero)
        CALL placeSpectrum(Sig12,stp, REAL(cC), omZero)
     CASE DEFAULT
        STOP "RAM: SOMETHING UNPREDICTED"
     END SELECT

     WRITE(0,*) 'omZero=',omZero,' C(1)=',REAL(C(1))
     IF(Mode /= ModeDecay) CALL placeSpectrum(Sigma, stp, REAL(C), omZero)
 !    CALL printCArr(Sigma%sigma(:,stp), Sigma%zeroOut, Emin, Emax, &
 !         &Sigma%deltaOm, omInc)
     CALL SYSTEM_CLOCK(COUNT=cnt,COUNT_RATE=crt)
     IF(lastDump+DumpFreq<REAL(cnt)/REAL(60*crt) .AND. Mode /= ModeDecay) THEN
        CALL saveSpectrum(Sigma,TRIM(ofName)//".spec", stp)
        lastDump = REAL(cnt)/REAL(60*crt)
     END IF
  END DO
  
  WRITE(0,*)'Saving to: ', TRIM(ofName)
  
  IF(Mode /= ModeDecay) CALL saveSpectrum(Sigma,TRIM(ofName)//".spec")
		
  IF(dumpDur) CALL dumpFunction(durAr, Sigma%zeroInc*au2ev,&
       &Sigma%deltaOm * Sigma%inStep*au2ev, TRIM(ofName)//".duration")

  WRITE(0,*) 'Done.'
  DEALLOCATE(durAr)
  CALL confFreeMem

  SELECT CASE(Mode)
  CASE(ModeRamDir)
     CALL saveSpectrum(SigDir,  TRIM(ofName)//".specDir")
     CALL saveSpectrum(SigCross,TRIM(ofName)//".specCross")
     SigDir%sigma = Sigma%sigma+SigDir%sigma+SigCross%sigma
     SigDir%modified = .TRUE.
     CALL saveSpectrum(SigDir, TRIM(ofName)//".specTotal")
     DEALLOCATE(SigDir%sigma, SigCross%Sigma, cA, cB)
  CASE(ModeRamSym)
     CALL saveSpectrum(Sig1, TRIM(ofName)//".spec1")
     CALL saveSpectrum(Sig2, TRIM(ofName)//".spec2")
     CALL saveSpectrum(Sig12,TRIM(ofName)//".spec12")
     DEALLOCATE(Sig1%sigma, Sig2%Sigma, Sig12%sigma, cA, cB, cC)
  END SELECT
  

DEALLOCATE(Sigma%sigma, C)

END PROGRAM RamanTwo
! end of main ------------------------------------------------------
