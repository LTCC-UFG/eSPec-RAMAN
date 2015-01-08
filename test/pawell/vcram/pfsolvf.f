! pfsolv.f90
! Solver procedures to TD improved Faris'es method for Raman scattering
! within narrow-band approximation. 
! Pawel Salek, pawsa@theochem.kth.se
! 980219
! 980422 : BLAS routines. Huge boost on Suns.
! 980617 : ESSL fft routines for IBM. In my case it degraded performance 
!          from 40s to 42s.
! NOTE
! 1. the array declarations are gradually upgraded from DIMENSION(:) :: p
!    to the form recommended by HIRLAM team  :: p(:)
!    The BLAS routines make code a little messy, always try to debug first 
!    the pure F90 code (undefine USE_BLAS).
! 2. dumping of dOmega shall be generalized and moved to IOProcessor. 
!    Probably unified also with dumpPotential.
! NOTE MODIFIED definition of duration time (only n1 returned, search for ###)
! 3. separate subroutines for calculations of: Absorption, Raman scattering,
!    Raman with direct term.
! REMARKS: 990311 - fixed formula mistake in the RamDir module.
!    The fourier transform has to be changed: the creation of c.c. data 
!    should happen in one place only.
! 990322 - handling of homonuclear molecules.
! 990518 - two energy scales introduced: kinetic energy and the binding 
!          energy.
MODULE PFSolver
  USE LapackInterface
  USE Lanczos
CONTAINS
  SUBROUTINE convertEnergyScale(omInc, vec, eneZero)
    USE IOProcessor, ONLY: KinOrBin
    IMPLICIT NONE
    REAL(mk), INTENT(IN)       :: omInc
    COMPLEX(mk), INTENT(INOUT) :: vec(:)
    REAL(mk), INTENT(OUT)      :: eneZero

    IF(KinOrBin == 'K') THEN
       ! workaround COMPAQ F90 compiler bug, DIGITAL Fortran 90, 5.2, ECO 01
       vec = vec(SIZE(vec):1:-1); vec = CSHIFT(vec,-1)
       !vec = CSHIFT(vec(SIZE(vec):1:-1),-1)
       eneZero = omInc
    ELSE
       eneZero = 0.0_mk
    END IF
  END SUBROUTINE convertEnergyScale
  
  ! calcAbsorption ===================================================
  ! calculates absorption from initial to a final state. Starts 
  ! calcCTLanczos. 
  ! the enrgy conversion is kind of messy, I haven't figured out any better 
  ! way to do it. Be _CAUTIONOUS_ when changing this stuff.
  ! energy scale must be in BINDE mode.
  ! VK: Coupling in the final states
  SUBROUTINE calcAbsorption(Fi, E0, omInc, dtSamp, C, omZero)
    USE IOProcessor, ONLY: MinX, MaxX, Vfi, QOper, Mass,VerbLvl,DumpSigmaT,&
         &ofName, KinOrBin, DumpFFT, DOper,&
         &dumpFunctionC, dumpFunction
    IMPLICIT NONE
    REAL(mk), INTENT(IN)     :: Fi(:)
    REAL(mk), INTENT(IN)     :: dtSamp, E0, omInc
    COMPLEX(mk), INTENT(OUT) :: C(:)
    REAL(mk), INTENT(out)    :: omZero

    INTEGER     :: k, low, high, N
    COMPLEX(mk) :: ft(SIZE(Fi)*SIZE(Vfi,2))
    REAL(mk)    :: dx

    N=SIZE(Vfi,1)
    dx = (MaxX-MinX)/(SIZE(Fi)-1)
    ft = (0.0, 0.0)
    omZero = MINVAL(Vfi(:,1))+E0 ! arbitrary shift, this might cause troubles!
    DO k=1,SIZE(Vfi,2) ! loop over states
       low = N*(k-1)+1
       high = k*N
       ft(low:high) = Fi ! define the N-th state wave packet
       IF(ASSOCIATED(DOper)) ft(low:high) = ft(low:high) * DOper(:,k)!
    ENDDO ! end of loop over states

    CALL calcCTLanczos(Vfi - omZero - E0, ft, dtSamp, minX, dx,mass, C, verbLvl)

    IF(dumpSigmaT) CALL dumpFunctionC(C, 0.0_mk, dtSamp,TRIM(ofName)//".Ct")
    CALL doFFT(C, dtSamp);     CALL convertEnergyScale(omInc,C, dx)

    IF(KinOrBin == 'K') WRITE(0,*) "KINET selected. Results uninterpretable."
    !omZero = omZero + dx
    IF(dumpFFT) CALL dumpFunction(ABS(C), 0.0_mk, dtSamp,&
         & TRIM(ofName)//".Cfft")
  END SUBROUTINE calcAbsorption
  
   ! calcDecay ========================================================
   ! dumpDur - says if the duration time should be calculated and dumped.
   ! the result is placed in duration
   SUBROUTINE calcDecay(Fi, E0, omInc, duration, dOmOut)
     USE IOProcessor, ONLY: MaxX, MinX, MasterP, au2ev, Vde, Mass, Gamma,&
          &DecayThreshold, VerbLvl, DumpDur, QOper, ev2au, DumpDOm, ofName,&
          &i2s, PrintFVBra, dumpFunction, DOper
     IMPLICIT NONE
     REAL(mk), INTENT(IN)  :: Fi(:)
     REAL(mk), INTENT(IN)  :: omInc, E0
     REAL(mk), INTENT(OUT) :: duration
     COMPLEX(mk),OPTIONAL  :: dOmOut(Size(Fi))

     COMPLEX(mk) :: dOm(SIZE(Fi)), acc(SIZE(Fi)*SIZE(Vde,2))
     REAL(mk)    :: dx, durationK, nrm
     INTEGER     :: i, pos(1), k
	   INTEGER     :: low, high, N !vks
     COMPLEX(mk) :: ft(SIZE(Fi)*SIZE(Vde,2)) !vk
     
     !; dOm=1.0
     dx = (MaxX-MinX)/(SIZE(Fi)-1)
     pos = MAXLOC(Fi)
     dOm = (0.0,0.0)
     duration = 0.0
     N=SIZE(Vde,1) !vk 

    DO k = 1, SIZE(Vde,2) ! VK
        IF(masterP) THEN
           WRITE(0,"('Resonance frequency V(R_0)-E0=',F9.4,' eV (pos=',I4,')')")&
                &(Vde(pos(1),k)-E0)*au2ev, pos
        END IF
        low = N*(k-1)+1 !vk
        high = k*N !vk
        ft(low:high) = Fi ! define the N-th state wave packet !vk
        IF(ASSOCIATED(DOper)) then
        	ft(low:high) = Fi * DOper(:,k) ! R-dependent transition matrix element
        endif
   		ENDDO ! end of loop over states

      nrm = norm(ft)
      
	    CALL integr(ft, Vde-E0-omInc, mass, minX, dx, Gamma,& !vk
          & DecayThreshold, acc, verbLvl, dumpDur, durationK) !vk

!============================= test block
!	do i=1, SIZE(ft)
!		write(11,*) i, abs(ft(i)), abs(acc(i))
!	enddo
!	write(0,*) ' >>> NORM FT, NORM(acc):', norm(ft), norm(acc)
!  write(0,*) '-----------------------------------------------------------'
!=============================

      duration = MAX(duration, durationK) !? where it is used?

			DO k = 1, SIZE(Vde,2) 
	      ! apply Q operator if this was requested in the input file      
				IF(ASSOCIATED(QOper)) THEN !R-dependent Qoper
		      low = N*(k-1)+1 !vk
 		      high = k*N !vk
         	WRITE(0,*) '<d|Q|d>/<d|d>=',&
               &DOT_PRODUCT(acc(low:high), QOper(:,k)*acc(low:high))/DOT_PRODUCT(acc(low:high), acc(low:high))
					acc(low:high) = acc(low:high)*QOper(:,k)
				END IF
			
				! print out wavepacket |\Psi(0)>	
!  	    PRINT *, '(Title "|d(omega)>, omega=',omInc*ev2au,'")'
!    	  print *, '(Title "State #',k,'")'
!      	IF (verbLvl/=0) CALL PrintFVBra(ABS(dOm(low:high))**2, Vde(:,k), minX, dx)
				dOm = dOm + acc(low:high)
			ENDDO

!============================= test block
!	do i=1, SIZE(Fi)*SIZE(Vde,2)
!		write(12,*) i, abs(acc(i))
!	enddo
!	write(0,*) ' >>> NORM dOm:', norm(acc), norm(dOm)
!=============================

			 ! dump result finally as a complex function
  	   IF(dumpDOm) THEN
    	    OPEN(21, ACTION="WRITE", &
      	       &FILE=TRIM(ofName)//".dOmega"//i2s(NINT(omInc/ev2au*100)))
        	WRITE(21,"('# |d(omega)>, omega=',F8.3,'eV')") omInc/ev2au
					WRITE(21,"((3(E13.4)))") ( &
!          	   &(/ minX+(i-1)*dx, REAL(dOm(i)), AIMAG(dOm(i)) /), i=1, SIZE(dOm))
                  &(/ minX+(i-1)*dx, abs(dOm(i)), AIMAG(dOm(i)) /), i=1, SIZE(dOm))
        	CLOSE(21)
	     END IF
!		 ENDDO  

     IF( PRESENT(dOmOut)) dOmOut = dOm
     !IF(verbLvl>1) THEN
     !   WRITE(0,"('The average position:',F10.3)") &
     !        & SUM( (/ ( (MinX+i*dx)*ABS(dOm(i+1))**2, i=0,N-1 ) /) )/&
     !              & SUM( ABS(dOm)**2 )
     !END IF
   
   END SUBROUTINE calcDecay
  
!   ! calcRaman ========================================================
   SUBROUTINE calcRaman(Fi, E0, omInc, dtSamp, C, omZero, duration)
     USE IOProcessor, ONLY: MaxX, MinX, Vfi, Mass, VerbLvl, &
          &DumpSigmaT, DumpFFT, EMax, EMin, EptNum, ofName, &
          &dumpFunction, vcdip !vk differ in Auger rate added
     IMPLICIT NONE
     REAL(mk), INTENT(IN)  :: Fi(:)
     REAL(mk), INTENT(IN)     :: omInc, E0, dtSamp
     COMPLEX(mk), INTENT(OUT) :: C(:)
     REAL(mk), INTENT(OUT)    :: omZero, duration

     COMPLEX(mk) :: dOm(SIZE(Fi)), dOm1(size(fi)*size(vfi,2))
     !INTEGER    ::  
	   INTEGER     :: low, high, N, k !vks
     REAL(mk)    :: dx
     dx = (MaxX-MinX)/(SIZE(Fi)-1)

     CALL calcDecay(Fi, E0, omInc, duration, dOm)
     ! THE SECOND PHASE....

!vibronic coupling in the final states: !vk
    N=SIZE(Vfi,1)
    dOm1 = (0.0, 0.0)
    DO k=1,SIZE(Vfi,2) ! loop over states
       low = N*(k-1)+1
       high = k*N
	! vk : different Auger decay to vc-final states
	IF(k==2) dOm=dOm*VCdip
       dOm1(low:high) = dOm ! define the N-th state wave packet
!       IF(ASSOCIATED(DOper)) ft(low:high) = ft(low:high) * DOper(:,k)!
    ENDDO ! end of loop over states

!     CALL calcCTLanczos(Vfi-E0, dOm, dtSamp, minX, dx,mass, C, verbLvl)
     CALL calcCTLanczos(Vfi-E0, dOm1, dtSamp, minX, dx,mass, C, verbLvl)
     IF(dumpSigmaT) CALL dumpFunction(ABS(C), 0.0_mk, dtSamp, &
          &TRIM(ofName)//".Ct")

     CALL doFFT(C, dtSamp); CALL convertEnergyScale(omInc,C,omZero)

     IF(dumpFFT) CALL dumpFunction(ABS(C), omInc,(EMax-EMin)/(eptNum-1),&
          &TRIM(ofName)//".Cfft")

   END SUBROUTINE calcRaman
  
!   ! calcRamSym =======================================================
!   ! computes Raman scattering for a symmetric, homonuclear molecule.
   SUBROUTINE calcRamSym(Fi, E0, omInc, dtSamp, C, C1,C2,C12, omZero, duration)
     USE IOProcessor, ONLY: MaxX, MinX, QOper, Vfi, Mass, VerbLvl, HomoParity
     IMPLICIT NONE
     REAL(mk), INTENT(IN)  :: Fi(:)
     REAL(mk), INTENT(IN)     :: omInc, E0, dtSamp
     COMPLEX(mk), INTENT(OUT) :: C(:), C1(:), C2(:), C12(:)
     REAL(mk), INTENT(OUT)    :: omZero, duration

     COMPLEX(mk)          :: dOm(SIZE(Fi)), c12Tmp(SIZE(C))
     REAL(mk)             :: dx
     COMPLEX(mk), POINTER :: myQ(:)
     
     dx = (MaxX-MinX)/(SIZE(Fi)-1)

     ! disable Q first...
     ! ==================
     IF(.NOT.ASSOCIATED(QOper)) THEN
        WRITE(0,"(A)") "ERROR: Q is not associated in RamSym"
        RETURN
     END IF
    
     myQ => QOper(:,1)! In this case assume that RamDyr will use one decay state;
     NULLIFY(QOper)
    
     ! proceed with Q switched off...
     ! ==================
     CALL calcDecay(Fi, E0, omInc, duration, dOm)

     ! compute C1 and the second half of the contribution to C12
     ! ================================
     CALL calcCTLanczos(Vfi-E0, myQ*dOm, dtSamp, minX, dx,mass, C1, verbLvl,&
          & C12, CONJG(myQ)*dOm)

     ! compute C2 and the first half of the contribution to C12
     ! ================================
     CALL calcCTLanczos(Vfi-E0, CONJG(myQ)*dOm, dtSamp, minX, dx,mass, C2, &
          &verbLvl, C12tmp, myQ*dOm)


     C12 = C12tmp+C12

     C = C1 +C2 + C12*HomoParity
     ! fourier transform and dumping....
     ! ====================
     CALL doFFT(C, dtSamp);       CALL convertEnergyScale(omInc, C,   omZero)
     CALL doFFT(C1, dtSamp);      CALL convertEnergyScale(omInc, C1,  omZero)
     CALL doFFT(C2, dtSamp);      CALL convertEnergyScale(omInc, C2,  omZero)
     CALL doFFT(C12, dtSamp);     CALL convertEnergyScale(omInc, C12, omZero)

     ! restore QOper - just to make it clean...
     ! ===============
     ! QOper(:,1) => myQ ! This ofcourse doesn't work. Copy again
     ALLOCATE(QOper(SIZE(myQ),1))
     QOper(:,1) = myQ
     DEALLOCATE(myQ)
   END SUBROUTINE calcRamSym

   ! calcRamDir ========================================================
   ! Must put cCross(-t) = cCross(t)^*
   ! note that omZero is set TWICE to the same value.
   !
   SUBROUTINE calcRamDir(Fi, E0, omInc, dtSamp, C, cDir, cCross, omZero, duration)
     USE IOProcessor, ONLY: MaxX, MinX, Emax, Emin, EptNum, Vfi, Mass,&
          & VerbLvl, ACoeff, DumpSigmaT, ofName, DumpFFT, ev2au, &
          & dumpFunction
     IMPLICIT NONE
     REAL(mk), INTENT(IN)     :: Fi(:)
     REAL(mk), INTENT(IN)     :: omInc, E0, dtSamp
     COMPLEX(mk), INTENT(OUT) :: C(:), cDir(:), cCross(:)
     REAL(mk), INTENT(OUT)    :: omZero, duration

     COMPLEX(mk) :: dOm(SIZE(Fi)), fun(SIZE(Fi)), tmpC(SIZE(C))
     REAL(mk)    :: potFin(SIZE(Fi,1)), dx, dltOm 

     INTEGER     :: lc, i

     lc = SIZE(C)/2+1
     dx = (MaxX-MinX)/(SIZE(Fi)-1)
     dltOm = (EMax-EMin)/EptNum
    
     CALL calcDecay(Fi, E0, omInc, duration, dOm)

     ! THE SECOND PHASE... 
     !====================
     fun    = FI
     !potFin = Vfi - E0
     
     !enSh = MINVAL(Vfi)-E0 !+deltaE
write(0,*)'NORM: RES, DIR ===============>', norm(dOm), norm(fi)

     ! compute wave-packet evolution of Psi(0) on final PES
     CALL calcCTLanczos(Vfi - E0 , dOm, dtSamp, minX, dx, mass, C,    &
     			&	verbLvl, cCross, fun)

     ! compute evolution of |0> on final PES
     CALL calcCTLanczos(Vfi - E0 , fun, dtSamp, minX, dx, mass, cDir, &
          & verbLvl, tmpC,   dOm)

     ! find out the cross term from cCross and tmpC
     cCross(1:lc) = CMPLX(0.0_mk,1.0_mk)*&
          &(ACoeff*tmpC(1:lc) -CONJG(ACoeff)*cCross(1:lc))
     ! FINAL FFT - all the crucial parts must be in the SECOND HALF
     !=============

     cDir = ABS(ACoeff)**2*cDir

     IF(dumpSigmaT) THEN
        CALL dumpFunction(ABS(C),       0.0_mk, dtSamp,&
          &TRIM(ofName)//".Ct")
        CALL dumpFunction(REAL(cDir),   0.0_mk, dtSamp,&
             &TRIM(ofName)//".CtDir")
        CALL dumpFunction(REAL(cCross), 0.0_mk, dtSamp,&
             &TRIM(ofName)//".CtCrossRe")
        CALL dumpFunction(AIMAG(cCross),0.0_mk, dtSamp,&
             &TRIM(ofName)//".CtCrossIm")
     END IF
    
     CALL doFFT(C, dtSamp); CALL convertEnergyScale(omInc,C,omZero)

     CALL doFFT(cCross, dtSamp);  CALL convertEnergyScale(omInc,cCross, omZero)
     CALL doFFT(cDir, dtSamp);    CALL convertEnergyScale(omInc,cDir,   omZero)

     IF(dumpFFT) THEN
        CALL dumpFunction(REAL(C),      omInc/ev2au,dltOm/ev2au,&
             &TRIM(ofName)//".Cfft")
        CALL dumpFunction(REAL(cDir),   omInc/ev2au,dltOm/ev2au,&
             &TRIM(ofName)//".CfftDir")
        CALL dumpFunction(REAL(cCross), omInc/ev2au,dltOm/ev2au,&
             &TRIM(ofName)//".CfftCross")
        CALL dumpFunction(AIMAG(cCross),omInc/ev2au,dltOm/ev2au,&
             &TRIM(ofName)//".CfftCrossIm")
     END IF
   END SUBROUTINE calcRamDir
  
   ! ===================================================================
   ! the result of integration is placed in res==dOm
   ! dumpDur - says if the duration is to be evaluated and dumped.
   ! the result is placed in duration
   ! decays thresholds: 0.005 is enough for the calculation of |d(omega)>
   ! is divided by 10 for calculation of duration.
   SUBROUTINE integr(fi, Vde, mass, minX, dx, Gamma, DecayThreshold, res, &
        & verb, dumpDur, duration)
     USE IOProcessor, ONLY: SubSize, au2ev, StepPrec, &
          & optPrintCoreDumped, dumpAvgOps,           &
          & PrintFVBra, Vspl !VK: put couplins here
     IMPLICIT NONE
     COMPLEX(mk), INTENT(IN)  :: Fi(:)
     REAL(mk), INTENT(IN)     :: Vde(:,:)
     REAL(mk), INTENT(IN)     :: mass, minX, dx, Gamma, DecayThreshold
     COMPLEX(mk), INTENT(OUT) :: res(SIZE(FI))
     INTEGER, INTENT(IN)      :: verb
     LOGICAL, INTENT(IN)      :: dumpDur
     REAL(mk), INTENT(OUT)    :: duration
  
     COMPLEX(mk) :: Qd(SIZE(Fi),SubSize)
     REAL(mk)    :: Td(SubSize,2)
     REAL(mk)    :: Zd(SubSize, SubSize)
     REAL(mk)    :: Dd(SubSize)
     COMPLEX(mk) :: Fde(SIZE(Fi)), tVec(SIZE(Fi))
     COMPLEX(mk) :: vl(Subsize), vc(Subsize)
     REAL(mk)    :: actt, dt, maxT, lastt, n1, n2, xavg, resAvg, nrm
     INTEGER     :: i, stepCnt, k, low, high, N
     complex(mk) :: FdeR(size(Vde,1)),resR(size(Vde,1))      
#ifdef USE_BLAS
     COMPLEX(mk), PARAMETER :: ONE=(1,0), ZERO=(0,0)
     COMPLEX(mk)            :: qTmp(SIZE(Qd,1),SIZE(Qd,2)),tmp(SubSize)
     COMPLEX(mk)            :: zTmp(SIZE(Zd,1),SIZE(Zd,2))
#endif
#ifdef SYS_CRAY
     REAL(mk),DIMENSION(12*SIZE(Fi))    :: wsave
     REAL(mk),DIMENSION(4*SIZE(Fi))     :: work
     i=0; CALL ccfft(0,size(Fde),1.0, Fde,Fde, WSAVE, work,i)
#elif defined(USE_ACML)
     REAL(mk), DIMENSION(2*(3*SIZE(Fi)+100)) :: WSAVE! size uncertain
     CALL zfft1d(0,SIZE(Fi),Fde,WSAVE(1),i)
#else
!     REAL(mk), DIMENSION(4*SIZE(Fi)+15) :: WSAVE !vkr
!     CALL ffti(SIZE(Fi),WSAVE(1)) !vkr
     REAL(mk), DIMENSION(4*SIZE(Fi)/SIZE(Vde,2)+15) :: WSAVE !vka
     CALL ffti(SIZE(Fi)/SIZE(Vde,2),WSAVE(1)) !vka
#endif
   	 N=size(Vde,1)
     res = 0.0
     tVec= 0.0
     stepCnt = 0
     maxT=-LOG(DecayThreshold)/Gamma
     actt=0.0; lastt=0.0
     IF (verb>=0) WRITE (0,&
          &"('integr       : decay threshold:',E9.2,' maxT=',F8.2)") &
          &DecayThreshold, maxT
     IF(verb>1) WRITE(0,"('V_min=',F13.3,' eV, V_max=',F13.3,' eV')")&
          &MINVAL(Vde)*au2ev, MAXVAL(Vde)*au2ev
     nrm=norm(Fi)
     write(0,*) '=====================================', nrm !vka
     Fde=Fi/nrm
     ! the run -----------------------------
     IF(verb>1) &
          & WRITE (0,"(A9,A9,A6,A8,A7,A7,A6,A8)")&
          &'mt',' t ','dt','||d>|',' <x>:','<res>','  stps','E'
     DO WHILE (actt<=maxT)     

       CALL GenKrylSp(Vde,Fde,dx,mass, Qd, Td, Zd, Dd, WSAVE,Vspl) !Vk 
        dt = estimStep(Td, StepPrec)
        vl = (/ (CMPLX(Gamma,Dd(i)),i=1,SubSize) /)
        vc = Zd(1,:)*EXP(-Gamma*actt)*((1-EXP(-vl*dt))/vl)
!#if defined(USE_BLAS)
!        tmp = MATMUL(Zd, vc )
!        CALL GEMV('N',SIZE(Qd,1), SIZE(Qd,2),ONE, Qd(1,1), SIZE(Qd,1),&
!             &tmp(1),1,ONE, res(1), 1)
!        tmp=MATMUL(Zd, Zd(1,:)*EXP(Dd*CMPLX(0.0_mk,-dt)) )
!        CALL GEMV('N',SIZE(Qd,1), SIZE(Qd,2),ONE, Qd(1,1), SIZE(Qd,1),&
!             &tmp(1),1,ZERO, Fde(1), 1)
!        IF (dumpDur) THEN
!           !vc= vc*actT+&
!           !    &EXP(-Gamma*actT)*(1.0-(1.0+vl*dt)*EXP(-vl*dt))/(vl**2)*Zd(1,:)
!           vc = (vc + EXP(-Gamma*actT)*(actT-(actT+dt)*EXP(-vl*dt))*Zd(1,:))/vl
!           vc = MATMUL(Zd, vc)
!           CALL GEMV('N',SIZE(Qd,1), SIZE(Qd,2),ONE, Qd(1,1), SIZE(Qd,1),&
!                &vc(1),1,ONE, tVec(1), 1)
!        END IF
!#else      
        res = res+MATMUL(Qd, MATMUL(Zd, vc ) )
          
        Fde = MATMUL(Qd, MATMUL(Zd,Zd(1,:)*EXP(Dd*CMPLX(0.0,-dt)) ) )
        IF (dumpDur) THEN
           vc = (vc + EXP(-Gamma*actT)*(actT-(actT+dt)*EXP(-vl*dt))*Zd(1,:))/vl
           tVec=tVec+MATMUL(Qd,vc)
           tVec=tVec+MATMUL(Qd,MATMUL(Zd,vc))
        END IF
!#endif
        actt=actt+dt
        stepCnt=stepCnt+1
				

        ! do some reporting...
				IF (actt>lastt) THEN
					lastt=lastt+10.0
					FdeR=(0.0,0.0); resR=(0.0,0.0)
          do k=1,size(Vde,2)
	        	IF(verb>=0.AND.EXP(-Gamma*actT)*ABS(Fde(k*SIZE(Fde)/size(Vde,2)))>0.5/(SIZE(Fde)/size(Vde,2))) &
  	        	&WRITE (0,"('PACKET ',I2,' ARRIVED TO the right EDGE! t=',F7.1,"//&
    	        &"', value=',E8.2)") k, actT, ABS(Fde(SIZE(Fde)))
						enddo
          IF(verb>1) THEN
 						DO k = 1, SIZE(Vde,2) !vk
			      	low = N*(k-1)+1; high = k*N
			        FdeR = FdeR + Fde(low:high)
			        resR = resR + res(low:high)
		        enddo  		
!              xavg = DOT_PRODUCT(Fde, &
!                   &(/ ((minX+(i-1)*dx)*Fde(i),i=1,SIZE(Fde)) /) )             
!              resAvg = SUM( (/ ( (minX+i*dx)*ABS(res(i+1))**2, i=0,size(Fi)-1 ) /) )/& !vkr
!                   & DOT_PRODUCT(res,res) !vkr
            xavg = DOT_PRODUCT(FdeR,(/ ((minX+(i-1)*dx)*FdeR(i),i=1,N) /) )             
            resAvg = SUM( (/ ( (minX+i*dx)*ABS(resR(i+1))**2, i=0,N-1 ) /) )/& !vkr
                 & DOT_PRODUCT(resR,resR) !vkr
            WRITE (0,"(F9.1,F9.2,F6.3,F8.3,F7.3,F7.3,I6,F10.4)") &
                 &maxT,actt,dt, norm(resR), xavg, resAvg, &
                 & stepCnt, Td(1,1)*au2ev
            PRINT "('(Title integr t=',F7.2,' dt=',F7.4,' E=',F8.4,'eV)')",&
                 & actt, dt, Td(1,1)*au2ev
            IF(optPrintCoreDumped) THEN
               CALL PrintFVBra(ABS(FdeR)**2 * EXP(-Gamma*actT),Vde(:,1), minX, dx) !vk
            ELSE
               CALL PrintFVBra(ABS(FdeR)**2,Vde(:,1), minX, dx) !vk
            END IF
          END IF
        ENDIF

        IF(dumpAvgOps) CALL dumpAverages(Fde)
     ENDDO
     IF (dumpDur) THEN
        n1=norm(tVec); n2=norm(res)
        duration=n1/n2 ! ### was n1/n2
        IF(verb>0) WRITE(0,"(A,F6.1,A,F10.2,A,E10.3)") &
             &"||d>|=",n2," |tVec|=",n1," tDur=",n1
     END IF
   END SUBROUTINE integr

  ! -----------------------------------------------------------------
  ! calcCTLanczos
  ! never shortens the time step in order to avoid unnecessary 
  ! generation of Lanczos space. In this way, the calculation of the 
  ! autocorrelation FUNCTION point takes ONLY constant overhead: 
  ! one matrix-matrix multiplication.
  ! DO NOT assume that Fi is normalized.
  SUBROUTINE calcCTLanczos(V, Fi, dtSamp, vSt, dx,mass, C, verb,&
       & auxCorr, auxV)
    USE IOProcessor, ONLY: SubSize, StepPrec, optTestDissociation,&
         & PrintFVBra, Vspl, minx, absln, abscf ! Get the interactions here.
    IMPLICIT NONE
    REAL(mk),    INTENT(IN) :: V(:,:)
    COMPLEX(mk), INTENT(IN) :: Fi(:)
    REAL(mk),    INTENT(IN) :: dtSamp, vSt, dx, mass
    COMPLEX(mk),INTENT(OUT) :: C(:)
    INTEGER,     INTENT(IN) :: verb
    COMPLEX(mk), OPTIONAL   :: auxCorr(SIZE(C)), auxV(size(Fi))
    
    
    COMPLEX(mk)            :: Fact(SIZE(Fi))
    REAL(mk)               :: actt, lastT, nrm
    REAL(mk)               :: dts, dt
    INTEGER                :: j, i, stepCnt,i1, N, Clen, low, high
    COMPLEX(mk)            :: Qf(size(Fi),SubSize)
    REAL(mk)               :: Tf(SubSize,2)   
    REAL(mk)               :: Zf(SubSize, SubSize)
    REAL(mk)               :: Df(SubSize)
    COMPLEX(mk)            :: tmp(SubSize)
    REAL(mk),ALLOCATABLE   :: absVec(:)
    LOGICAL                :: auxCorrFunc, doExit
#ifdef USE_BLAS
    COMPLEX(mk), PARAMETER :: ONE=(1,0), ZERO=(0,0)
    COMPLEX(mk)            :: vl(SubSize), tmpFi(SIZE(Qf,1))
    COMPLEX(mk)            :: zTmp(SIZE(Zf,1),SIZE(Zf,2))
#endif    
#ifdef SYS_CRAY
    REAL(mk),DIMENSION(12*SIZE(Fi)/SIZE(V,2))    :: wsave
    REAL(mk),DIMENSION(4*SIZE(Fi)/SIZE(V,2))     :: work
   j=0; CALL ccfft(0,size(Fact)/SIZE(V,2),1.0, Fact,Fact, WSAVE, work,j)
#elif defined(USE_ACML)
    REAL(mk), DIMENSION(2*(3*SIZE(Fi)/SIZE(V,2)+100)) :: WSAVE! size uncertain
    CALL zfft1d(0,SIZE(Fi)/SIZE(V,2),Fi,WSAVE(1),i)
#else   
!    REAL(mk), DIMENSION(4*SIZE(Fi)+15) :: WSAVE
!    CALL ffti(SIZE(Fi),WSAVE(1))
    REAL(mk), DIMENSION(4*SIZE(Fi)/SIZE(V,2)+15) :: WSAVE
    CALL ffti(SIZE(Fi)/SIZE(V,2),WSAVE(1))
#endif

    nrm=SQRT(DOT_PRODUCT(Fi,Fi))
    IF(nrm<1e-8)&
         & STOP "Norm of |d> less than 1e-10 -> numerical noise on output."
    auxCorrFunc = PRESENT(auxCorr) .AND. PRESENT(auxV)
    IF(auxCorrFunc) auxCorr = 0.0

    C=0.0; actt=0; lastT=0; stepCnt=0; Clen=SIZE(C)
    N=SIZE(Fi); i1=N/absln
!vk test    ALLOCATE(absVec(i1)); absVec=EXP(-0.008*(/ (REAL(j)/i1,j=0,i1-1) /) )
 ALLOCATE(absVec(i1)); absVec=EXP(-abscf*(/ (REAL(j)/(i1-REAL(j)),j=0,i1-1) /)) 
 !/(i1-(/ (REAL(j),j=0,i1-1) /)) )
    Fact=Fi/nrm
!vk test
		write(0,*) "=======>>>>>>", size(absVec), i1
!		write(30,*) (/ (REAL(j),j=0,i1-1) /)
		do j=1,i1
		write(20,*) j, absVec(j)
		enddo
!stop

    C(1)=nrm**2
    IF(auxCorrFunc) auxCorr(1) = DOT_PRODUCT(auxV,Fi)
    IF(verb>=0) WRITE(0,"('calcCTLanczos: clen=',I6,' dtSamp=',F7.3,"//&
         & "' maxT=',F11.1,' |d(0)|=',F11.3)") &
         &Clen, dtSamp,Clen*dtSamp/2.0, nrm
    ! calcCTLanczos: MAIN LOOP -------------------------------------------
    ! print header if required...

    IF(verb>1) WRITE(0,"(A8,A7,A6,A7,A11,A9, 2A10)")&
         & 'Time ', 'dt ','step ',' <x> ','|C(t)| ','Arg C(t)','D(1)=E', 'NRM'

    TheLoop: DO j=2,Clen/2+1
       DO WHILE(actt< (j-1)*dtSamp)
          nrm=nrm*norm(Fact) 
          CALL GenKrylSp(V,Fact, dx, mass, Qf,Tf,Zf,Df, WSAVE, Vspl)
					Qf(N-i1+1:N,:)= Qf(N-i1+1:N,:)*SPREAD(absVec,2,SubSize)
          dt=estimStep(Tf, StepPrec)
!#if defined(USE_BLAS)
!          zTmp=Zf; vl=Zf(1,:)*EXP(CMPLX(0,-dt)*Df)
!          CALL GEMV('N',SIZE(Zf,1), SIZE(Zf,2),ONE, zTmp(1,1), SIZE(Zf,1),&
!               &  vl(1), 1, ZERO, tmp(1),1)
!          CALL GEMV('N',SIZE(Qf,1), SIZE(Qf,2),ONE,Qf(1,1),SIZE(Qf,1),&
!               &tmp(1), 1, ZERO, Fact(1),1)
!#else
          Fact = MATMUL(Qf,MATMUL(Zf, Zf(1,:)*EXP(CMPLX(0,-dt)*Df)))
!#endif
          stepCnt=stepCnt+1
          actt=actt+dt
       END DO
       dts=(j-1)*dtSamp-(actt-dt)
!#if defined(USE_BLAS)
!       zTmp=Zf; vl=Zf(1,:)*EXP(CMPLX(0,-dts)*Df)
!       CALL GEMV('N',SIZE(Zf,1), SIZE(Zf,2), ONE, zTmp(1,1), SIZE(Zf,1),&
!            & vl(1), 1, ZERO, tmp(1),1)
!       CALL GEMV('N',SIZE(Qf,1), SIZE(Qf,2), ONE, Qf(1,1), SIZE(Qf,1),&
!            & tmp(1),1, ZERO, tmpFi(1),1)
!       C(j)=nrm*DOT_PRODUCT(Fi,tmpFi)
!#else
       tmp=MATMUL(Zf, Zf(1,:)*EXP(CMPLX(0,-dts)*Df))
       C(j)=nrm*DOT_PRODUCT(Fi, MATMUL(Qf,tmp))
!#endif
       IF(auxCorrFunc) auxCorr(j) = nrm*DOT_PRODUCT(auxV, MATMUL(Qf,tmp))
       IF (verb>1 .AND. lastT < (j-1)*dtSamp) THEN
! vk : change output step: 50.0 => 16.0
          lastT=lastT+50.0
          WRITE (0,"(F8.2,F7.3,I6,F7.3,E11.3,F9.2,F10.2, E11.3)") actt, dt, stepCnt, &
               & SUM( (/ ( (vSt+i*dx)*ABS(Fact(i+1))**2, i=0,N-1 ) /) ), &
               & ABS(C(j)), ATAN2(AIMAG(C(j)),REAL(C(j)))*180.0/PI,&
               & Df(1), nrm
!!!!!!!!!!!!
          PRINT ("('(Title ""CTLanczos t=',F8.1,' dt=',F7.4,' E[eV]=',F8.3,"//&
               &"'|C(t)|=',E11.3,' nrm=',F10.6,'"")')"), actt, dt, Tf(1,1),&
               & ABS(C(j)), nrm
          DO i=1, SIZE(V,2) !loop over states
            low = (i-1)*SIZE(Fi)/SIZE(V,2)+1
            high = i*SIZE(Fi)/SIZE(V,2)
    !        ! Print N times. Must be changed. gbrat can take more 
    !        ! than three columns.
	         CALL PrintFVBra(ABS(Fact(low:high)),V(:,i), vSt, dx)
   !         CALL sleep(1)
         ENDDO
!!!!!!!!!!!!         
       ENDIF
       IF(optTestDissociation) THEN
          doExit = ABS(C(j)) < 5e-6*nrm
          ! more strict rules are needed when the cross term is computed...
          IF(doExit .AND. auxCorrFunc) &
               & doExit = ABS(auxCorr(j)) < MAXVAL(ABS(auxCorr(1:j-1)))*5e-6
          IF(doExit) THEN
             IF (verb>=0) THEN
                WRITE (0,*) 'Dissociation observed. Loop finished.'
                WRITE (0,"('Time=',F8.1,' C(t)=',E9.3)") actt,ABS(C(j))
             END IF
             EXIT TheLoop
          END IF
       END IF
    END DO TheLoop

    DEALLOCATE(absVec)

    IF(verb>=0) WRITE(0,"('Last dt:',E10.3,' Avg. dt:',E10.3,"//&
         & "' Total',I5,' steps. steps/clen: ',F7.4)")&
         & dt, actt/stepCnt, stepCnt, (2.0*stepCnt)/Clen

    !these last two lines _theoretically are redundant, since the doFFT
    ! takes care about the second halv of the arrays but 
    ! they remain here just for nice sigma(t) dumps...
    C(Clen/2+2:Clen)=CONJG(C(Clen/2:2:-1))
    IF(auxCorrFunc)  auxCorr(Clen/2+2:Clen)=CONJG(auxCorr(Clen/2:2:-1))

  END SUBROUTINE calcCTLanczos


  ! doFFT ============================================================
  ! the execution time is not crucial and we can can simplify it (in hope
  ! to improve portability) by calling only the standard fftPack routine.
  ! SETS the second half of vector.
  SUBROUTINE doFFT(C, dtSamp)
    USE IOprocessor, ONLY: FinalGamma, MasterP, au2ev
    IMPLICIT NONE
    COMPLEX(mk), INTENT(INOUT) :: C(:)
    REAL(mk), INTENT(IN)       :: dtSamp

    INTEGER  :: Clen, i
    REAL(mk) :: wnd(SIZE(C)), gm
    LOGICAL, SAVE :: msgP = .TRUE.
    
!#ifdef SYS_CRAY
!    REAL(mk),DIMENSION(12*SIZE(C))    :: wsave
!    REAL(mk),DIMENSION(4*SIZE(C))     :: work
!    i=0; CALL ccfft(0,size(C),1.0, C,C, WSAVE, work,i)
!#else
    REAL(mk), DIMENSION(4*SIZE(C)+15) :: WSAVE
    CALL ffti(SIZE(C),WSAVE(1))
!#endif

    Clen = SIZE(C)
    C(Clen/2+2:clen) = CONJG(C(clen/2:2:-1))

    IF(FinalGamma <= 0.0) THEN
       gm = 17.0/REAL(Clen,mk)
       IF(masterP.AND.msgP) &
            &WRITE(0,"('The FinalGamma has been set to:',F12.7,' eV')")&
            &gm/dtSamp*au2ev
       msgP = .FALSE.
    ELSE
       gm = FinalGamma*dtSamp
    ENDIF
    wnd(1:Clen/2+1)=EXP( (/ (-REAL(i)*gm, i=0, Clen/2) /) )
    !wnd(1:Clen/2+1)=EXP( (/ (-31.0*(REAL(i)/Clen)**2, i=0, Clen/2) /) )
    wnd(Clen/2+2:Clen) = wnd(Clen/2:2:-1)

    C=C*wnd
    !WRITE (0,*) 'The 0th element is ',ABS(SUM(C))
!#ifdef SYS_CRAY
!    i=0;CALL ccfft(1,Clen,1.0,C(1),C(1),WSAVE(1),work(1),i)
!    C=C*dtSamp
!#else
    CALL fftb(Clen, C(1), WSAVE(1))
    C=C*dtSamp !C*dtSamp verified on angus, 2000.04.26
!#endif
  END SUBROUTINE doFFT

  ! dumpAverages =====================================================
  ! dumps average values of some operators
  SUBROUTINE dumpAverages(Fi)
    USE IOprocessor, ONLY: MinX, MaxX
    IMPLICIT NONE
    COMPLEX(mk), INTENT(IN) :: Fi(:)

    INTEGER i
    REAL(mk) :: dx
    dx = (MaxX-MinX)/(SIZE(Fi)-1)
    WRITE(0,"('<r>: ',F10.6,' au')") DOT_PRODUCT(Fi,&
         & (/ ( (MinX+dx*i)*Fi(i+1), i=0, SIZE(FI)-1) /) )
  END SUBROUTINE dumpAverages
  
end MODULE PFSolver
