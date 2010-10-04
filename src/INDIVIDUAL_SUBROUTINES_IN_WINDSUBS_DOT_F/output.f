      SUBROUTINE OUTPUT
C=========================================================
C save the output for graphics/analysis
C Written by Erxuan Fu, GSO,URI,11/3/94
C=========================================================
      INCLUDE 'comblk.h'
      include 'comblk1.h'

      REAL TMP1(IM,JM),TMP2(IM,JM),TMP3(IM,JM,KB),jjulday
      REAL DB1(KB-1)
      DIMENSION TB1(IM,JM,KB),U1(IM,JM,KB),V1(IM,JM,KB)
      DIMENSION SB1(IM,JM,KB),OHC(IM,JM),MLDTH(IM,JM)
      DIMENSION TB1S(IM,JM)
      INTEGER yyear
      CHARACTER FN*15, DOUT*8
      integer*4 sstartdate,date
C
      CALL DATE2DAY(year,julday,startdate)
CC
       julday=julday+time+1.e-5
CC
      call day2date(year,julday,date)
      WRITE(DOUT,'(I8.8)') date
C
      FN = 'TXY.'//DOUT
      OPEN(39,FILE=FN,STATUS='UNKNOWN',form='unformatted')
      WRITE(39) TAUX
      WRITE(39) TAUY
      CLOSE(39)
      WRITE(6,*) ' --> ',FN
C
      RETURN
      END
