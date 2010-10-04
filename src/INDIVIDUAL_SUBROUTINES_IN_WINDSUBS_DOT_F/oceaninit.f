      SUBROUTINE OCEANINIT
C
         INCLUDE 'comblk.h'    
         INCLUDE 'comblk1.h'    
c
      REARTH=6371.E3
      LATMIN=10.
      LATMAX=47.5
      LONGMIN=-98.5
      LONGMAX=-50.
c
      PI=3.1415927E0
      RAMP=1.E0
      SMALL=1.E-10
      TIME0=0.E0
      BETA=1.98E-11
      GRAV=9.806E0
      UMOL=2.E-5
C
      IINT=0
      TIME=0.
c
C--- Model contral parameters----------------------------
C--- see file PARAMETERS.inp for meanings of parameters
191   format(a15)
c---------- falk 07-06-01 use description of file in script
C     OPEN(10,FILE='PARAMETERS.inp',STATUS='OLD')
      READ(10,*) MODE
      READ(10,*) NBC
      READ(10,*) TIME0
      READ(10,*) IMAY
      READ(10,*) sstsource
      READ(10,*) startdate
      READ(10,191) pathname
      READ(10,191) pathname1
      READ(10,*) DTI
      READ(10,*) ISPLIT
      READ(10,*) NREAD
      READ(10,*) IPRTH1
      READ(10,*) INOWINDH,VC
      READ(10,*) IDAMP,ISMOTH,SMH
      READ(10,*) wndg
      READ(10,*) IHOURS
      IF(NREAD.GT.0) READ(10,'(A)') FRSTI
      rewind(10)
C
      ISPADV=5
      TPRNU=1.   
      SMOTH=0.1
      HORCON=0.1
      TBIAS=10.
      SBIAS=35.
      RHO_0=1024.e0
      ISWTCH=100000
      IPRTD2=1
c
C--------------------------------------------------------------------------
C     READ NAMELIST (DISCONNECTED HERE) FOR PROBLEM PARAMETERS: 
C       MODE = 2; 2-D CALCULATION (BOTTOM STRESS CALCULATED IN ADVAVE)
C              3; 3-D CALCULATION (BOTTOM STRESS CALCULATED IN PROFU,V)
C              4; 3-D CALCULATION WITH T AND S HELD FIXED
C       NBC = TYPE OF THERMO. B.C.; SEE SUBROUTINE PROFT
C       DTI = INTERNAL TIME STEP
C       DTE = EXTERNAL TIME STEP
C       ISPLIT = DTI/DTE
C       NREAD=0, NO RESTART INPUT FILE; NREAD=1, RESTART
C       IPRTD1 = PRINT INTERVAL IN DAYS; AFTER ISWTCH DAYS, PRINT
C                  INTERVAL = IPRTD2
C       IDAYS = LENGTH OF RUN IN DAYS
C       ISPADV = STEP INTERVAL WHERE EXTERNAL MODE ADVECTIVE TERMS ARE
C                  NOT UPDATED
C       HORCON = CONSTANT IN SMAGORINSKY HORIZONTAL DIFFUSIVITY
C       TPNU = HORIZONTAL DIFFUSIVITY PRANDTL NUMBER
C       SMOTH = CONSTANT IN TIME SMOOTHER TO PREVENT SOLUTION SPLITTING
C
C     READ(5,PRMTR)
C--------------------------------------------------------------------------
c
      DTE=DTI/FLOAT(ISPLIT)
      DTE2=DTE*2
      DTI2=DTI*2
      IPRINT=IPRTH1*3600/INT(DTI)
      ISWTCH=ISWTCH*24*3600/INT(DTI)   
      IEND=IHOURS*3600/INT(DTI)
C
C----------------------------------------------------------------------
C             ESTABLISH PROBLEM CHARACTERISTICS
C          ****** ALL UNITS IN M.K.S. SYSTEM ******
C      F,BLANK AND B REFERS TO FORWARD,CENTRAL AND BACKWARD TIME LEVELS.
C----------------------------------------------------------------------
C
      DAYI=1.E0/86400.E0
C
C--- ASSIGNING SPACING AND CORIOLIS
C--- STEP IN LONGITUDE(PHI) IS (LONGMAX-LONGMIN)/(IM-1)
C--- STEP IN LATTITUDE(LAMBDA) IS (LATMAX-LATMIN)/(JM-1)
C
      DO 14 J=1,JM
      DO 14 I=1,IM
      PHI=LONGMIN+(I-1)*(LONGMAX-LONGMIN)/(IM-1)
      PHI=2.*3.1415927*PHI/360.
      LAMBDA=LATMIN+(J-1)*(LATMAX-LATMIN)/(JM-1)
      LAMBDA=2.*3.1415927*LAMBDA/360
      COR(I,J)=4.*3.1415927/(24.*60.*60.)*SIN(LAMBDA)
      DX(I,J)=REARTH*COS(LAMBDA)*(LONGMAX-LONGMIN)/(IM-1)
      DX(I,J)=DX(I,J)*2.*3.1415927/360.
      DY(I,J)=REARTH*(LATMAX-LATMIN)/(JM-1)
      DY(I,J)= DY(I,J)*2.*3.1415927/360.
      ART(I,J)=DX(I,J)*DY(I,J)
  14  CONTINUE
C
      DO 13 J=2,JM
      DO 13 I=2,IM
      ARU(I,J)=.25*(DX(I,J)+DX(I-1,J))*(DY(I,J)+DY(I-1,J))
  13  ARV(I,J)=.25*(DX(I,J)+DX(I,J-1))*(DY(I,J)+DY(I,J-1))
      DO 121 J=1,JM
      ARU(1,J)=ARU(2,J)
      ARV(1,J)=ARV(2,J)
 121  continue
      DO 122 I=1,IM
      ARU(I,1)=ARU(I,2)
      ARV(I,1)=ARV(I,2)
 122  continue
C
      call date2day(year,julday,startdate)
      print*,'startdate=',startdate,';  julday=',julday
      julday=julday+ihours/24.
      call day2date(year,julday,enddate)
c
      TIME0=0
      TIME=TIME0
      CALL OUTPUT
      print *,' OCEANINIT end'
c
      return
      end

