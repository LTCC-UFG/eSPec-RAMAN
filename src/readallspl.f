      subroutine readallspl(Fil,IFLTH,X,Y,T,BCOEFRe,BCOEFIm,
     &    XKNOT,YKNOT,TKNOT,NP,nf,KX,NORM,stfil,dim)
      IMPLICIT NONE
c------------------------------------------
c
c Purpose:
c This routine reads a wavepacket propagation into matrices,
c   then generates a 3D spline representation for it
c
c
c     Goiania, 25th of september of 2014
c     Vinicius Vaz da Cruz
c-------------------------------------------
      INTEGER I,J,NF,K,IFLTH,LDF,MDF,stfil
      INTEGER II,JJ,KK,FPR(4),PRTRMSWHOLE,NOREG
      INTEGER NP(*),KX,WHOLEGRID
      REAL*8 X(*), Y(*), T(*)
      REAL*8 WORK(5),NORM
      REAL*8 checkRe, checkIm,DBS3VL,RMSD,RMSDRE,RMSDIM
      LOGICAL file_exists,file_open
      CHARACTER*14 wwork,wp,wo,pol
      CHARACTER*6 DIM
      CHARACTER*30 Fil, Filen
c--- 3D SPLINE PARAMETERS
      REAL*8 TKNOT(*), YKNOT(*), XKNOT(*)
      REAL*8 Re(NP(1),NP(2),NF),Im(NP(1),NP(2),NF)
      REAL*8 BCOEFRe(NP(1),NP(2),NF),BCOEFIm(NP(1),NP(2),NF)
c      REAL*8 BCOEFRe(LDF,MDF,*),BCOEFIm(LDF,MDF,*)
c----
      INTEGER TBEGIN,TEND
      REAL*8 RATE,TDIFF

      CALL SYSTEM_CLOCK(TBEGIN,RATE)
      write(*,*) "Reading wavepacket files: "
      DO K = 1, NF, 1
c--------reads Psi(r,t) for a fixed R value.
         CALL FILENAME(Filen,Fil,IFLTH,K+stfil-1)
         INQUIRE(FILE=Filen, EXIST=file_exists)
         IF (file_exists) THEN
            OPEN(unit=10,file=Filen)
         ELSE
            write(*,*) "could not open file ", Filen
            STOP
         ENDIF   

d     write(*,*)'>>',fpr(1),fpr(2),fpr(3),fpr(4)

         READ(10,*) wwork, wp, T(K), pol
         T(K) = T(K) * 41.3413745758D+0 ! fs -> au
         write(*,'(A20,A3,ES17.10,A5)')FILEN,',',T(K),'a.u.'
         DO I=1,NP(2),1
            DO J=1, NP(1),1
               IF(DIM(1:3).EQ.'.2D')THEN
                  READ(10,*) WORK(1),WORK(2),WORK(3),WORK(4)
                  X(I)=WORK(1)
               ELSE
                  READ(10,*)WORK(2),WORK(3),WORK(4)
               ENDIF
               Y(J)=WORK(2)
               Re(J,I,K)=WORK(3)*NORM
               Im(J,I,K)=WORK(4)*NORM
            ENDDO
         ENDDO
         CLOSE(10)
c---------------------end of K loop
      ENDDO

      CALL SYSTEM_CLOCK(TEND,RATE)
      TDIFF = REAL(TEND - TBEGIN)/REAL(RATE)
     
c      WRITE(*,*) 'the read data will be printed to files, (debug)'
c      DO K=1,NF,1
c         CALL FILENAME(Filen,'check_read',9,K)
c         CALL PRTWP(Re(1:LDF,1:MDF,K),Im(1:LDF,1:MDF,K),X,Y,T(K),NXR,NYR
c     &        ,FILEN)
c      ENDDO

      write(*,*)
      write(*,*)'All data has been read!'
      write(*,'(A9,F12.5,A10)') 'it took ',TDIFF, 'seconds.'
      write(*,*)'------------------------------------------------------'
      write(*,*)


c--------spline procedure
      write(*,*)
      write(*,*)'Generating spline coefficients matrix'
      write(*,*)
      
      CALL SYSTEM_CLOCK(TBEGIN,RATE)
d      WRITE(*,*)'sizes', NXR,NYR,NF
d      WRITE(*,*)'sizes', SIZE(X),SIZE(Y),SIZE(T)
   
      CALL DBSNAK(NP(1), Y, KX, YKNOT)
      IF(DIM(1:3).EQ.'.2D') CALL DBSNAK(NP(2), X, KX, XKNOT)
      CALL DBSNAK(NF, T, KX, TKNOT)

d      DO I =1, NYR,1
d         write(1,*) Y(I)
d      ENDDO
d      write(1,*)
d      DO I =1, NXR,1
d         write(1,*) X(I)
d      ENDDO
d      write(1,*)
d      DO I =1, NF,1
d         write(1,*) T(I)
d      ENDDO

      
         
      IF(DIM(1:3).EQ.'.2D')THEN
c     write(*,*) 'real part'
         CALL DBS3IN(NP(1),Y,NP(2),X,NF,T,Re,NP(1),NP(2),
     &        KX,KX,KX,YKNOT,XKNOT,TKNOT,BCOEFRe)
c     write(*,*) 'imaginary part'
         CALL DBS3IN(NP(1),Y,NP(2),X,NF,T,Im,NP(1),NP(2),
     &        KX,KX,KX,YKNOT,XKNOT,TKNOT,BCOEFIM)
      ELSE
c     write(*,*) 'real part'
         CALL DBS2IN(NP(1),Y,NF,T,Re(1:NP(1),1,1:NF),NP(1),
     &        KX,KX,YKNOT,TKNOT,BCOEFRe(1:NP(1),1,1:NF))
c     write(*,*) 'imaginary part'
         CALL DBS2IN(NP(1),Y,NF,T,Im(1:NP(1),1,1:NF),NP(1),
     &        KX,KX,YKNOT,TKNOT,BCOEFIM(1:NP(1),1,1:NF))  
      ENDIF

      CALL SYSTEM_CLOCK(TEND,RATE)
      TDIFF = REAL(TEND - TBEGIN)/REAL(RATE)

      write(*,*)
      write(*,*)'3D spline matrices calculated!'
      write(*,'(A9,F12.5,A10)') 'it took ',TDIFF, 'seconds.'
      write(*,*)
c      write(*,*)'estimating error in interpolation:'
c      write(*,*)
c      CALL ERRORSPL(FIL,IFLTH,Re,Im,BCOEFRE,BCOEFIM,Y,X,T
c     &    ,NYR,NXR,NF,YKNOT,XKNOT,TKNOT,KX,PRTRMSWHOLE)

      write(*,*)
      write(*,*)'Finished 3D spline section!'
      write(*,*)'------------------------------------------------------'
      
      END

c----------- routine to perform the error analysis in the spline interpolation
      SUBROUTINE ERRORSPL(FIL,IFLTH,Re,Im,BCOEFRE,BCOEFIM,Y,X,T
     &     ,NYR,NXR,NF,YKNOT,XKNOT,TKNOT,KX,PRTRMSWHOLE)
      IMPLICIT NONE
      INTEGER I,J,K,NF,NXR,NYR,IFLTH,M,PRTRMSWHOLE,KX
      CHARACTER*30 Filen,Fil
      REAL*8 RMSD,RMSDRE,RMSDIM,ERRMAXRE,ERRMAXIM,RMSDALL
      REAL*8 Re(NYR,NXR,*),Im(NYR,NXR,*),BCOEFRE(NYR,NXR,*)
      REAL*8 YKNOT(*),XKNOT(*),TKNOT(*)
c      REAL*8 BCOEFIM(LDF,MDF,*),X(*),Y(*),T(*)
      REAL*8 BCOEFIM(NYR,NXR,*),X(*),Y(*),T(*)

      IF(NF.GT.5)M = INT(NF/5)

      write(*,*)'computing RMSE in relation to individual files'
      write(*,'(A38, A15,A3,A15)')'','Real part','','Imag. part'
      DO K=1,NF,M
d      k=1
d        write(1,*)K,'real'
         RMSDRE=RMSD(Re(1:NYR,1:NXR,K),BCOEFRE,Y,X,T(K),YKNOT,XKNOT
     &        ,TKNOT,KX,NYR,NXR,NF,ERRMAXRE)
d         write(1,*)K,'imaginary' 
         RMSDIM=RMSD(Im(1:NYR,1:NXR,K),BCOEFIM,Y,X,T(K),YKNOT,XKNOT
     &        ,TKNOT,KX,NYR,NXR,NF,ERRMAXIM)
         write(*,'(A35,I4,X2 2ES17.10)')'Root mean square errors 
     &for file',K,RMSDRE,RMSDIM
         write(*,'(A35,I4,X2 2ES17.10)')'maximum error for file',
     &        K,ERRMAXRE,ERRMAXIM
c         CALL FILENAME(Filen,'check_spl',9,K)
c         CALL PRTWPSPL(BCOEFRE,BCOEFIM,X,Y,T(K),NXR,NYR,NF,FILEN,
c     &        YKNOT,XKNOT,TKNOT,KX)
      ENDDO

      IF(PRTRMSWHOLE.EQ.0)THEN
         write(*,*)
         write(*,*)'computing RMSE in relation to the whole data'
         RMSDRE=RMSDALL(Re(1:NYR,1:NXR,1:NF),BCOEFRE,Y,X,T,YKNOT,XKNOT
     &        ,TKNOT,KX,NYR,NXR,NF,ERRMAXRE)
         RMSDIM=RMSDALL(Im(1:NYR,1:NXR,1:NF),BCOEFIM,Y,X,T,YKNOT,XKNOT
     &        ,TKNOT,KX,NYR,NXR,NF,ERRMAXIM)
         write(*,'(A35,X2 2ES17.10)')'Root mean square errors' 
     &        ,RMSDRE,RMSDIM
         write(*,'(A35,X2 2ES17.10)')'maximum error',
     &        ERRMAXRE,ERRMAXIM
         
         IF(RMSDRE.GT.1.0D-5 .OR.RMSDIM.GT.1.0D-5)THEN
            write(*,*) 'Error in the interpolation is too large
     &(> 1.0D-5), please check your calculation!'
            STOP
         ENDIF
      ENDIF

      END


c-------------- routine to print two matrices M1(J,I) and M2(I,J) into a file
      SUBROUTINE PRTWP(M1,M2,X,Y,T,NX,NY,FILEN)
      IMPLICIT NONE
      INTEGER I,J,NX,NY,LDF,MDF
      CHARACTER*30 FILEN
      REAL*8 M1(NY,*),M2(NY,*),X(*),Y(*),T
      OPEN(unit=12,FILE=FILEN)

      WRITE(12,*)'# time =',T
      DO I=1,NX,1
         DO J=1,NY,1
            WRITE(12,'(4ES20.10)')X(I),Y(J),M1(J,I),M2(J,I)
         ENDDO
      ENDDO
      CLOSE(12)
      END

c-------------- routine to print wavepacket into a file from its spline representation
      SUBROUTINE PRTWPSPL(BCOEFRE,BCOEFIM,X,Y,T,NX,NY,NF,FILEN,
     &  YKNOT,XKNOT,TKNOT,KX)
      IMPLICIT NONE
      INTEGER I,J,NX,NY,NF,KX
      CHARACTER*30 FILEN
      REAL*8 BCOEFRE(NY,NX,*),BCOEFIM(NY,NX,*),DBS3VL
      REAL*8 YKNOT(*),XKNOT(*),TKNOT(*),WORK(2)
      REAL*8 X(*),Y(*),T
      OPEN(unit=12,FILE=FILEN)

      WRITE(12,*)'# time =',T
      DO I=1,NX,1
         DO J=1,NY,1
            WORK(1)=dbs3vl(Y(J),X(I),T,KX,KX,KX,
     &           yknot,xknot,tknot,ny,nx,nf,bcoefre)
            WORK(2)=dbs3vl(Y(J),X(I),T,KX,KX,KX,
     &           yknot,xknot,tknot,ny,nx,nf,bcoefim)
            WRITE(12,'(4ES20.10)')X(I),Y(J),WORK(1),WORK(2)
         ENDDO
         WRITE(12,*)
      ENDDO
      CLOSE(12)
      END

c------------------ routine to compute the root mean square between
c------------------ a section of a matrix M1, and the values generated by the
c------------------ spline coefficients matrix
      FUNCTION RMSD(M,MBCOEF,Y,X,T,YKNOT,XKNOT,TKNOT,
     &     KX,NYR,NXR,NF,ERRMAX)
      IMPLICIT NONE
      INTEGER I,J,NT,KX,KY,KZ
      INTEGER NXR,NYR,NF
      REAL*8 VAR,RMSD,MSPL,M(NYR,*),ERRMAX
      REAL*8 X(*),Y(*),T,XKNOT(*),YKNOT(*),TKNOT(*)
      REAL*8 MBCOEF(NYR,NXR,*),DBS3VL,VAL
    
      NT=NYR*NXR
      VAR = 0.0D+0
      ERRMAX = 0.0D+0
    
      DO I=1,NXR,1
         DO J=1,NYR,1
            MSPL=dbs3vl(Y(J),X(I),T,kx,kx,kx,
     &           yknot,xknot,tknot,nyr,nxr,nf,mbcoef)
            VAL = (M(J,I) - MSPL)**2 
d            write(1,'(5ES27.15)')X(I),Y(J),M(J,I),MSPL,VAL
            IF(DSQRT(VAL).GT.ERRMAX) ERRMAX=DSQRT(VAL/(NT*1.0D+0))
            VAR = VAR + VAL
         ENDDO
d         write(1,*)
      ENDDO

      RMSD = DSQRT(VAR/(NT*1.0D+0))

      RETURN
      END

 
c------------------ routine to compute the root mean square between
c------------------ a matrix M1, and the values generated by the
c------------------ spline coefficients matrix
      FUNCTION RMSDALL(M,MBCOEF,Y,X,T,YKNOT,XKNOT,TKNOT,
     &     KX,NYR,NXR,NF,ERRMAX)
      IMPLICIT NONE
      INTEGER I,J,K,NT,KX,KY,KZ
      INTEGER NXR,NYR,NF
      REAL*8 VAR,RMSDALL,MSPL,M(NYR,NXR,*),ERRMAX
      REAL*8 X(*),Y(*),T(*),XKNOT(*),YKNOT(*),TKNOT(*)
      REAL*8 MBCOEF(NYR,NXR,*),DBS3VL,VAL
      
      NT=NYR*NXR*NF
      VAR = 0.0D+0
      ERRMAX = 0.0D+0

      DO K=1,NF
         DO I=1,NXR,1
            DO J=1,NYR,1
               MSPL=dbs3vl(Y(J),X(I),T(K),kx,kx,kx,
     &              yknot,xknot,tknot,nyr,nxr,nf,mbcoef)
               VAL = (M(J,I,K) - MSPL)**2 
d     write(1,'(5ES27.15)')X(I),Y(J),M(J,I),MSPL,VAL
               IF(DSQRT(VAL).GT.ERRMAX) ERRMAX=DSQRT(VAL/(NT*1.0D+0))
               VAR = VAR + VAL
            ENDDO
d     write(1,*)
         ENDDO
      ENDDO

      RMSDALL = DSQRT(VAR/(NT*1.0D+0))

      RETURN
      END



c------------------- routine to determine filename

      SUBROUTINE FILENAME(Filen,Fil,IFLTH,K)
      IMPLICIT NONE
      INTEGER K,IFLTH,LEN_TRIM
      CHARACTER*30 Filen,Fil
      CHARACTER*4 CHNUM

      WRITE(CHNUM,'(I4)')K
      IF(K.LT.10)THEN
         Filen = Fil(1:IFLTH)//'000'//CHNUM(4:4)//'.dat'
      ELSEIF(K.LT.100)THEN
         Filen = Fil(1:IFLTH)//'00'//CHNUM(3:4)//'.dat'
      ELSEIF(K.LT.1000)THEN
         Filen = Fil(1:IFLTH)//'0'//CHNUM(2:4)//'.dat'
      ELSEIF(K.LT.10000)THEN
         Filen = Fil(1:IFLTH)//CHNUM(1:4)//'.dat'
      ENDIF

      END

c------------------- routine to read whole wavepacket file

      SUBROUTINE READWHOLE(Fil,IFLTH,K,Re,Im,X,Y,NX,NY,T)
      IMPLICIT NONE
      INTEGER I,J,K,IFLTH,NX,NY
      REAL*8 Re(NY,NX), Im(NY,NX), X(NX),Y(NY),T
      CHARACTER*30 Fil, FILEN,wwork,wp,pol
      CHARACTER*4 CHNUM

      CALL FILENAME(FILEN,FIL,IFLTH,K)

      OPEN(UNIT=25,FILE=FILEN)
      READ(25,*) wwork, wp, T, pol
      DO I=1,NX,1
         DO J=1,NY,1
            READ(25,*) X(I), Y(J), Re(J,I),Im(J,I)
         ENDDO
      ENDDO
      CLOSE(25)
      END

      


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   this was just to check whether the c pointer was being properly arranged into a fortran 3D matrix
      SUBROUTINE CHECKMATRIX(WPERE,WPEIM,NXR,NYR,NE,X,Y,E)
      implicit none
      INTEGER I,J,K,L,NXR,NYR,NE
      REAL*8 WPERE(NYR,NXR,*),WPEIM(NYR,NXR,*),X(NXR),Y(NYR),E(NE)
      
      DO I=1,NXR,1
         DO J=1,NYR,1
            WRITE(1,*) 'X,Y ', X(I),Y(J) 
            DO K=1,NE,1
               WRITE(1,'(3ES15.7)')E(K),WPERE(J,I,K),WPEIM(J,I,K)
            ENDDO
            WRITE(1,*)
         ENDDO
      ENDDO

      END
