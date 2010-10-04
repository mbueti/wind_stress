      PARAMETER (IMA=66,JMA=66)
c      COMMON /TOCN/  TOSTEP,TOCEAN,TOHOUR,TOCAT,ISM,PRT,PRD,IOSTEP 
c      REAL MINUTS1,MINUTS2
c      REAL DTE2,DTI2,DAYI, julday
c      real DT1,DT2,DT4,dayi,days,GRAV,SMOTH,RDRAG,time0
c      REAL ISPI,ISP2I,TMP,RHO_0
c      integer iend,year 
c      integer*4 startdate,enddate
c      CHARACTER*40 FRSTI, FRSTO
c      CHARACTER*15 pathname,pathname1

      COMMON/connect1/SWRAD(im,jm),SWRADs(im,jm),
     1     UG(im,jm,KB),VG(im,jm,KB),RDRAG
      common/misc/NBC,nbc2d(im,jm),nbc2ds(im,jm),mode,SMOTH
     1     ,IHOURS,INOWINDH,IPRTH1,TIME0,dayi,DT1,DT2,DT4,days
      common/dating/julday,startdate,enddate,year,pathname,pathname1

      real DT1,DT2,DT4,dayi,days,GRAV,SMOTH,RDRAG,time0
c      REAL ISPI,ISP2I,TMP,RHO_0
c      integer iend,year
      integer*4 startdate,enddate
      CHARACTER*40 FRSTI, FRSTO
      CHARACTER*15 pathname,pathname1
c
c      COMMON/connect1/UTB(IM,JM),VTB(IM,JM),UTF(IM,JM),VTF(IM,JM),
c     1     ADVX(IM,JM,KB),ADVY(IM,JM,KB),ADVUA(IM,JM),ADVVA(IM,JM),
c     2     TSURF(IM,JM),SSURF(IM,JM),WINDX2(IM,JM),WINDY2(IM,JM),
c     2     DRHOX(IM,JM,KB),DRHOY(IM,JM,KB),TRNU(IM,JM),TRNV(IM,JM),
c     3     ADVUU(IM,JM),ADVVV(IM,JM),WINDX1(IM,JM),WINDY1(IM,JM),
c     4     SWRAD(IM,JM),SSTIN(IM,JM),IDAMP,ISMOTH,IHOURS,
c     5     MINUTS1,MINUTS2,DTE2,DTI2,IEND,NREAD,IPRTH1,INOWINDH,
c     6     DAYI,ISPI,ISP2I,FRSTI,FRSTO,SMH,MODE,TIME0,ISPLIT,IMAY
c	    
      COMMON/ATMOCN1/XSU(IM),XSV(IM),XST(IM),YSU(JM),YSV(JM),YST(JM), 
     1                               TMA(IM,JM)   

      COMMON/ATMOCN2/SST(IMA,JMA),XSTR(IMA,JMA),YSTR(IMA,JMA),                   
     1           HFLX(IMA,JMA),XLON(IMA),YLAT(JMA),
     2           RFSW(IMA,JMA)               
      real*4 SST,XSTR,YSTR,HFLX,XLON,YLAT,RFSW
C
c      REAL LAMBDA,PHI,REARTH,LATMIN,LATMAX,LONGMIN,LONGMAX,CMP
c
c      COMMON/sphere/REARTH,LATMIN,LATMAX,LONGMIN,LONGMAX
c
c      common/misc/NBC,nbc2d(im,jm),nbc2ds(im,jm),RHO_0,pi,small,
c     1            beta,ispadv,smoth,horcon,iswtch,iprtd2
c
c      common/dating/julday,startdate,enddate,year,pathname,pathname1
