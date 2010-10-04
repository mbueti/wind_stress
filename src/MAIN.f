      PROGRAM MAIN
C$$$  MAIN PROGRAM DOCUMENTATION BLOCK
C
C MAIN PROGRAM: WINDONLY
C   PRGMMR: Richard Yablonsky          ORG: NP22         DATE: 2009-12-30
C
C ABSTRACT: This program runs only the wind stress generation component
C   of the Princeton Ocean Model configured for the United domain.
C   The entire code is split into two main subroutines: OCEANINIT and
C   OCEANSTEP. OCEANINIT reads all the input data and initializes model
C   variables. OCEANSTEP run the model for one time step.
C
C PROGRAM HISTORY LOG:
C   --------  BLUMBERG/MELLOR - Original implementation for POM
C   99-06-02  frolov - restructured, reconfigured and modified
C   09-12-30  yablonsky - modified for wind stress generation only
C
C INPUT FILES:
C   UNIT   10    PARAMETERS.inp file, contains configuration parameters
C   UNIT   15    storm history file, file name specified in PARAMETERS.inp
C
C OUTPUT FILES:
C   UNIT   39    binary wind stress data TXY.YYMMDDHH
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 77
C
C$$$
      include 'comblk.h'
      include 'comblk1.h'

c-------------- falk 01-05-04 specify param. for avrtau
      tauavr=0.
      taumax=0.
      awucon=1.
      bwucon=0.
      migtau=0


      CALL OCEANINIT
      DO IINT=1,IEND
        CALL OCEANSTEP(0)
      END DO
      
      STOP
      END
