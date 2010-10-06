      subroutine avrtau(xln,ylt,mig)
      parameter(timesm4=12.,timesm=12.,RAVR=100.e3)
c--------- This subr. is called from WIND in phase4 (mig=0)
c---------         and from atmos2ocean.f (mig=1) in coupled run
      include 'comblk.h'
      include 'TVARY.h'
      REAL LATMIN,LATMAX,LONGMIN,LONGMAX
      COMMON/sphere/REARTH,LATMIN,LATMAX,LONGMIN,LONGMAX
      real xln,ylt,tavr,counter,x1,y1,x0,y0,r,deltax,deltay
      real xlnc,yltc,tauavrp,taumaxp
      real tauavr,taumax,awucon,bwucon
      real RAVR,RRCT,xrct,yrct
      integer irct,jrct
      integer migtau
c
c------------- in comblk.h included common/tau/
c------------- also the next 5 parameters incuded in RST file
c     common/tau/ tauavr,taumax,awucon,bwucon,migtau
c-------------
      if(mig.eq.1) then
       print *,'begin avrtau: couped model'
      else
       print *,'begin avrtau: phase4'
      end if
c-------------
c
c     write(6,201) LATMIN,LATMAX,LONGMIN,LONGMAX
 201  format('avrtau: LATMIN,LATMAX,LONGMIN,LONGMAX=',4f7.2)
c
      pi=3.1415927
c------- save previous tauavr, taumax
      tauavrp=tauavr
      taumaxp=taumax
c-------  for coupled run use TC position from TVARY.h
c-------  for phase4 run use TC position: (xln,ylt)
      if(mig.eq.1) then
       xlnc=poslon
       yltc=poslat
      else
       xlnc=xln
       yltc=ylt
      end if
c
      tavr=0.0
      counter=0.0
      taumax=0.0
      RRCT=1.e8
      irct=1000
      jrct=1000
      do j=1,jm
       do i=1,im
        x1=(LONGMIN+float(I-1)*(LONGMAX-LONGMIN)/float(IM-1))*pi/180.
        y1=(LATMIN+float(J-1)*(LATMAX-LATMIN)/float(JM-1))*pi/180.
        x0=xlnc*pi/180.
        y0=yltc*pi/180.
        DELTAX=REARTH*COS(y0)*(x1-x0)
        DELTAY=REARTH*(y1-y0)
        r=SQRT(DELTAX**2+DELTAY**2)
        if(r.lt.RAVR) then
          tauabs=sqrt(wusurf(i,j)**2+wvsurf(i,j)**2)
          if(tauabs*fsm(i,j).gt.taumax) taumax=tauabs
          tavr=tavr+tauabs*fsm(i,j)
          counter=counter+fsm(i,j)
        end if
        if(r.lt.RRCT) then
         RRCT=r
         irct=i
         jrct=j
         xrct=x1*180./pi
         yrct=y1*180./pi
        end if
       end do
      end do
      if(counter.gt.0.) then
        tauavr=tavr/counter
      else
        tauavr=0.0
      end if
c
      if(mig.eq.1.and.migtau.eq.0) then
c--------- falk 08-19-03 use taumax instead of tauavr
c      if(tauavr.gt.tauavrp) then
c       awucon=tauavrp/tauavr
c       bwucon=(tauavr-tauavrp)/tauavr
       if(taumax.gt.taumaxp) then
        awucon=taumaxp/taumax
        bwucon=(taumax-taumaxp)/taumax
       else
        awucon=1.
        bwucon=0.
       end if
c-------------
       print *,' avrtau: first step in coupled model'
       print *,'migtau,mig=',migtau,mig
       write(6,101) tauavrp,tauavr,awucon,bwucon
 101   format(' tauavrp,tauavr,awucon,bwucon=',4(1PE10.2))
c-------------
       migtau=1
      end if
c
      if(mig.eq.1) then
       wucon=awucon+SIN(time*24./timesm*pi*0.5)*bwucon
       if(time*24..gt.timesm) wucon=1.
      else
       wucon=SIN(time*24./timesm4*pi*0.5)
       if(time*24..gt.timesm4) wucon=1.
      end if
c
      do j=1,jm
        do i=1,im
         wusurf(i,j)=wusurf(i,j)*wucon
         wvsurf(i,j)=wvsurf(i,j)*wucon
         taux(i,j)=taux(i,j)*wucon
         tauy(i,j)=tauy(i,j)*wucon
        end do
      end do
c-------------
c----------- falk 09-12-05 change output
      if(MOD(IINT,24).EQ.0) then
      if(mig.eq.1) then
       print *,'avrtau: couped model'
      else
       print *,'avrtau: phase4'
      end if
c-------------
      write(6,102) time*24,xlnc,yltc
 102  format('time*24,xlnc,yltc=',3f7.2)
      print *,'closest point to the center'
      write(6,204) xrct,yrct,RRCT,irct,jrct
 204  format('xrct,yrct,RRCT,irct,jrct=',2f7.2,f10.0,2i7)
      write(6,103) tauavrp,tauavr,taumaxp,taumax
 103  format('tauavrp,tauavr,taumaxp,taumax=',4(1PE10.2))
      write(6,202) awucon,bwucon,wucon
 202  format('  awucon,bwucon,wucon=',3f10.4)
c     write(6,203) timesm4,timesm
 203  format(   'timesm4,timesm=',2f7.2)
c-------------
      end if
      return
      end
c
c------------------- falk 06-21-05 use new SST assimilation procedure
c
	SUBROUTINE DATE2DAY(year,julday,date)
      integer*4 date
      integer dat2day(12),dat2dayl(12),day,month,year,hour
      real julday
      real*8 tmp
      data dat2day/31,28,31,30,31,30,31,31,30,31,30,31/
      data dat2dayl/31,29,31,30,31,30,31,31,30,31,30,31/

      year=int(date/1000000.)
      month=nint(100*(date/1000000.-int(date/1000000.)))
      julday=0
        if(mod(year,4).eq.0) then
      do n=1,month-1
      julday=julday+dat2dayl(n)
      end do
        else
      do n=1,month-1
      julday=julday+dat2day(n)
      end do
        end if
      julday=julday+nint(100*(date/10000.-int(date/10000.)))
      hour=date-nint(date/100.)*100
      julday=julday+float(hour)/24.

      return
      end
	SUBROUTINE DAY2DATE(year,julday,date)
      integer*4 date
      integer dat2day(12),dat2dayl(12),day,month,year,year1,hour
      real julday,julday1
      real*8 tmp
      data dat2day/31,28,31,30,31,30,31,31,30,31,30,31/
      data dat2dayl/31,29,31,30,31,30,31,31,30,31,30,31/

      if(int(julday).gt.365+int(1./(mod(year,4)*100.+1.))) then
      julday1=julday-365-int(1./(mod(year,4)*100.+1.))
      year1=year+1
      else
      julday1=julday
      year1=year
      end if
      day=0
      n=1
        if(mod(year1,4).eq.0) then  
      do while(day+dat2dayl(n).lt.int(julday1))
      day=day+dat2dayl(n)
      n=n+1   
      end do
        else
      do while(day+dat2day(n).lt.int(julday1))
      day=day+dat2day(n)
      n=n+1   
      end do
        end if
      month=n
      day=int(julday1-day)
      hour=nint((julday1-int(julday1))*24.)
      date=year1*1000000+month*10000+day*100+hour
c
      return
      end
      FUNCTION EXPWND(R,RMAX,Rref18,Rref26,WSMAX,WS18,WS26)
      real Rref18,Rref26,R,RMAX,WSMAX,WS18,WS26
      real b,expwnd
c 
       r1=0.5*(Rref18+Rref26)
       WS=0.5*(WS18+WS26)
       if(Rref18.le.0.) then
       r1=Rref26
       WS=WS26
       end if
       if(Rref26.le.0.) then
       r1=Rref18
       WS=WS18
       end if
c
         if(R.GE.RMAX) then
         b=(RMAX-r1)/log(WS/WSMAX)                     
         expwnd=WSMAX*exp((RMAX-R)/b)
         else
         expwnd=R*WSMAX/RMAX
         end if
c
      return
      end
      SUBROUTINE INTERP1d(mask,n,ni,x,y,xi,yi)
c--------------------------------------------------------------
c   This subroutine determines ni values yi at the points xi 
c   interpolating between the n values y at the points x
c   values equal to mask are ignored
c--------------------------------------------------------------
      real x(n),y(n),xi(ni),yi(ni),tmp(500)
      real cmp,mask
      integer n,ni,ii
c
      do i=1,ni
        if((xi(i).gt.x(n).and.xi(i).gt.x(1)).or.
     *  (xi(i).lt.x(1).and.xi(i).lt.x(n))) then
          if(xi(i).gt.x(n).and.xi(i).gt.x(1)) then
            if(x(n).gt.x(1)) then
              yi(i)=y(n)
            else
              yi(i)=y(1)
            end if
          else
            if(x(n).gt.x(1)) then
              yi(i)=y(1)
            else
              yi(i)=y(n)
            end if
          end if
        else
          do j=1,n-1
            tmp(j)=(xi(i)-x(j))*(xi(i)-x(j+1))
          end do
          do j=1,n-1
            if(tmp(j).le.0) ii=j
          end do
          if(y(ii).eq.mask.or.y(ii+1).eq.mask) then
            yi(i)=mask
          else
            yi(i)=(y(ii)*abs(x(ii+1)-xi(i))+y(ii+1)*abs(xi(i)-x(ii)))/
     1      abs(x(ii+1)-x(ii))
          end if
        end if
      end do
c
      return
      end
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

      SUBROUTINE OCEANSTEP(icoupling)
C----------------------------------------------------------------------
C     THIS SUBROUTINE DOES ONE STEP IN TIME OF THE OCEAN MODEL
C----------------------------------------------------------------------
C
      INCLUDE 'comblk.h'
      INCLUDE 'comblk1.h'
c
      TIME=DAYI*DTI*FLOAT(IINT)+TIME0
      RAMP=1.E0
C
        do j=2,jmm1
          do i=2,imm1
            wtsurf(i,j)=0.e0
            swrad(i,j)=0.e0
          end do
        end do
c
        do j=1,jm
          do i=1,im
            taux(i,j)=0.e0
            tauy(i,j)=0.e0
            wusurf(i,j)=0.e0
            wvsurf(i,j)=0.e0
          enddo
        enddo
c
        if((time*24.0).lt.float(inowindh)) then
          call wind(pathname,startdate)
          if(pathname1.ne.pathname) then
            call wind(pathname1,startdate)
          end if
        else
          do j=1,jm
            do i=1,im
              taux(i,j)=0.e0
              tauy(i,j)=0.e0
              wusurf(i,j)=0.e0
              wvsurf(i,j)=0.e0
            enddo
          enddo
        endif
c
C--------------------------------------------------------------------
C           BEGIN PRINT SECTION
C--------------------------------------------------------------------
c
      IF(MOD(IINT,IPRINT).NE.0.) GO TO 7000
C
      CALL OUTPUT
C
 7000 CONTINUE
C
C----------------------- END PRINT SECTION -----------------------------
C
      return
      END
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
      SUBROUTINE VERINTERP(n,ni,x,y,xi,yi)
c--------------------------------------------------------------
c   This subroutine determines ni values yi at the points xi 
c   interpolating between the n values y at the points x
c--------------------------------------------------------------
      real x(n),y(n),xi(ni),yi(ni),tmp(500)
      real cmp
      integer n,ni,ii

            do i=1,ni
      if((xi(i).gt.x(n).and.xi(i).gt.x(1)).or.
     *   (xi(i).lt.x(1).and.xi(i).lt.x(n))) then
        if(xi(i).gt.x(n).and.xi(i).gt.x(1)) then
            if(x(n).gt.x(1)) then
              yi(i)=y(n)
            else
              yi(i)=y(1)
            end if
        else
            if(x(n).gt.x(1)) then
              yi(i)=y(1)
            else
              yi(i)=y(n)
            end if
        end if
      else
        do j=1,n-1
        tmp(j)=(xi(i)-x(j))*(xi(i)-x(j+1))
        end do
c
         do j=1,n-1
           if(tmp(j).le.0) then
             ii=j
           end if
         end do
c
        yi(i)=(y(ii)*abs(x(ii+1)-xi(i))+y(ii+1)*abs(xi(i)-x(ii)))
     1        /abs(x(ii+1)-x(ii))
      end if
            end do

      return
      end
      SUBROUTINE WIND(filename,startdate)
      INCLUDE 'comblk.h'
      INTEGER PLN
      PARAMETER(PLN=100)
c
      integer hour, lat, long, mx, rmw
      integer*4 startdate,date
      integer day,month,year
      integer garb(5),Rd1(4),Rd2(4)
      character*19 name
      character*15 filename
      character*1 letter
c
      DIMENSION X(PLN),Y(PLN),TM(PLN),PRES(PLN),PRES0(PLN),
     *      RMAXa(PLN),WSPMAX(PLN),
     *      R18v(PLN,4),R26v(PLN,4),Rref18v(5),Rref26v(5),alphv(5)
      DIMENSION RAD(14),WS(14),RADM(14),WSM(14),ANGL(14)
      REAL CMP,T1,T2,F0,F1,L0,L1,REARTH,R,A7,B,E,DELP,x0,y0
      REAL LATMIN,LATMAX,LONGMIN,LONGMAX
      REAL DELTAX,DELTAX1,DELTAY,DELTAY1,DXDY,julday
      COMMON/sphere/REARTH,LATMIN,LATMAX,LONGMIN,LONGMAX

      DATA RM,R0/60.E3,480.E3/,RHO_0/1024.E0/
      DATA RAD/0.,.4,.7,.8,.95,1.,1.35,2.7,4.05,5.4,6.75
     * ,8.1,10.8,13.5/
      DATA WS/0.,.1,.5,.8,.95,1.,.97,.72,.54,.44,.4,.36
     * ,.27,.23/
      DATA ANGL/0.,2.,4.,6.,7.,7.,14.,23.,24.,22.,
     * 21.,21.,21.,21./
c     print*,'In subroutine WIND ...'
c
      WIND_SCALE=1.0
      ROA=1.28
      RMAX=50.E3
      PI=3.1415927
      E=exp(1.)
c
c----------------------- Reading message file --------------------
  17  format(A19,I6,1x,I4,1x,I3,1x,I4,1x,I3,1x,I3,3I5,
     * 1x,i2,1x,I3,1x,I4,1x,I4,1x,I4,1x,I4,1x,I4,
     * 1x,I4,1x,I4,1x,I4)
c     print*,'reading file ',filename
      open(15,file=filename,status='old')
      end=0.
      I=0
      do while(end.eq.0)
        read(15,17) name,date,hour,lat,long,garb,mx,rmw,Rd1,Rd2
c----------- falk 09-12-05 change output
       if(MOD(IINT,24).EQ.0) then
        print*,'reading file ',filename
        write(6,17) name,date,hour,lat,long,garb,mx,rmw,Rd1,Rd2
       end if
        if(date.eq.0) goto 20
        I=I+1
c
        date=date*100+hour/100
        call date2day(year,julday,date)
c
        TM(i)=julday
        X(i)=-long/10.
        Y(i)=lat/10.
        PRES(i)=float(garb(3))
        PRES0(i)=float(garb(4))
        WSPMAX(i)=float(mx)
        RMAXa(i)=float(rmw)
        do n=1,4
          n1=n+(1-mod(n,2))*sign(2,3-n)
          R18v(i,n)=Rd1(n1)
          R26v(i,n)=Rd2(n1)
          if(wspmax(i).le.26.or.R26v(i,n).le.RMAXa(i)) R26v(i,n)=-999
          if(wspmax(i).le.18.or.R18v(i,n).le.RMAXa(i)) R18v(i,n)=-999
          if(R26v(i,n).gt.R18v(i,n)) R26v(i,n)=-999
        end do
      end do
  20  end=1.
      ipmax=I
c----------- falk 09-12-05 change output
      if(MOD(IINT,24).EQ.0) then
      print*,'Number of hurricane path datapoints read: ',ipmax
c---------------- falk 01-05-04 change printing
c     print*,'tm=',(tm(i),i=1,ipmax)
      write(6,101) (tm(i),i=1,ipmax)
 101  format(10f7.2)
      print*,'year=',year
      print*,'startdate=',startdate
      end if
      close(15)
c
c--------------------- Calculating starting day -----------------
c
      call date2day(year,julday,startdate)
      do i=1,ipmax
        TM(i)=TM(i)-julday
      end do
C++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  INTERPOLATION OF HURRICANE PATH TO DETERMINE THE CURRENT POSITION

      CMP=TM(ipmax)
      if(cmp.lt.time.or.time.lt.tm(1)) then
c---------- falk 01-05-04 change printing
        print*,'NO HURRICANE PATH DATA FOR THIS TIME'
        print*,'   time=',time
        print*,'   tm(1)=',tm(1)
        print*,'   tm(ipmax)=',tm(ipmax)
        print*,'   day=',julday
        RETURN
      end if

      call verinterp(ipmax,1,TM,X,TIME,F0)
      call verinterp(ipmax,1,TM,Y,TIME,L0)
      xcen=F0
      ycen=L0
      call verinterp(ipmax,1,TM,PRES,TIME,PRES1)
      call verinterp(ipmax,1,TM,PRES0,TIME,PRES2)
      DELP=(PRES2-PRES1)*100.
      call verinterp(ipmax,1,TM,WSPMAX,TIME,WSMAX)
      prsmin=PRES1
      wndmax=WSMAX
      WSMAX=WSMAX*WIND_SCALE
      WS18=18*WIND_SCALE
      WS26=26*WIND_SCALE
      call verinterp(ipmax,1,TM,RMAXa,TIME,RMAX)
      RMAX=RMAX*1.e3
      do n=1,4
        call interp1d(-999.0,ipmax,1,TM,R18v(1,n),TIME,Rref18v(n))
        if(Rref18v(n).ne.-999) Rref18v(n) = Rref18v(n)*1.e3
        call interp1d(-999.0,ipmax,1,TM,R26v(1,n),TIME,Rref26v(n))
        if(Rref26v(n).ne.-999) Rref26v(n) = Rref26v(n)*1.e3
        alphv(n) = (n-1)*pi/2
      end do
      do n=2,6
        n1=mod(n-1,4)+1
        nm1=mod(n-2,4)+1
        np1=mod(n,4)+1
        if(Rref18v(n1).eq.-999) then
          if(Rref18v(nm1).ne.-999) then
            if(Rref18v(np1).ne.-999) then
              Rref18v(n1)=0.5*(Rref18v(nm1)+Rref18v(np1))
            else
              Rref18v(n1)=Rref18v(nm1)
            end if
          else
            if(Rref18v(np1).ne.-999) then
              Rref18v(n1)=Rref18v(np1)
            else
              Rref18v(n1)=-999
            end if
          end if
        end if
        if(Rref26v(n1).eq.-999) then
          if(Rref26v(nm1).ne.-999) then
            if(Rref26v(np1).ne.-999) then
              Rref26v(n1)=0.5*(Rref26v(nm1)+Rref26v(np1))
            else
              Rref26v(n1)=Rref26v(nm1)
            end if
          else
            if(Rref26v(np1).ne.-999) then
              Rref26v(n1)=Rref26v(np1)
            else
              Rref26v(n1)=-999
            end if
          end if
        end if
      end do
c
      Rref18v(5) = Rref18v(1)
      Rref26v(5) = Rref26v(1)
      alphv(5) = alphv(4)+pi/2
c
c----------- falk 09-12-05 change output
      if(MOD(IINT,24).EQ.0) then
      print*,'Time=',Time
      print*,'Current hurricane position (x,y): ',f0,l0
      print*,'WSMAX=',WSMAX,'; DELP=',DELP,'; RMAX=',RMAX
      print*,'Rref18v=',Rref18v
      print*,'Rref26v=',Rref26v
      end if
c
      x0=f0
      y0=l0
      F0=F0*2.*PI/360
      L0=L0*2.*PI/360
C--- F0,L0 ARE THE CARRENT (LONGITUDE,LATTITUDE) COORDINATES OF THE HURRICANE
C
C--- CALCULATING UTX AND UTY (HURRICANE SPEED)
c
      cmp=tm(ipmax)
      do i=1,ipmax
        if(abs(tm(i)-time).le.cmp) then
          cmp=abs(tm(i)-time)
          ii=i
        end if
      end do
c----------- falk 01-05-04 not to go out bnd
      if((tm(ii)-time).le.0.and.ii.ne.ipmax) then
        t1=tm(ii)
        t2=tm(ii+1)
        x1=x(ii)
        x2=x(ii+1)
        y1=y(ii)
        y2=y(ii+1)
      else
        t2=tm(ii)
        t1=tm(ii-1)
        x2=x(ii)
        x1=x(ii-1)
        y2=y(ii)
        y1=y(ii-1)
      end if
      deltax1=rearth*cos(l0)*(x2-x1)*2.*pi/360
      deltay1=rearth*(y2-y1)*2.*pi/360
      utx=deltax1/((t2-t1)*24.*3600.)
      uty=deltay1/((t2-t1)*24.*3600.)
c
c----------- falk 09-12-05 change output
      if(MOD(IINT,24).EQ.0) 
     *print*,'utx,uty: ',utx,uty
C
C--- CALCULATING PARAMETERS FOR WIND PROFILE FORMULA
c
      B=WSMAX**2*E*ROA/DELP
      A7=RMAX**B
c----------- falk 09-12-05 change output
      if(MOD(IINT,24).EQ.0) 
     *print*,'B= ',B
c
      do i=1,14
        RADM(I)=RMAX*RAD(I)
      end do
C
      DO 350 J=1,JM
      DO 351 I=1,IM

C  CALCULATING U-WIND STRESS FOR I,J POINT OF U-VELOCITY GRID
      F1=LONGMIN+(I-1)*(LONGMAX-LONGMIN)/((IM-0.5)-1)
      F1=F1*2.*PI/360
      L1=LATMIN+(J-1)*(LATMAX-LATMIN)/(JM-1)
      L1=L1*2.*PI/360
      DELTAX=REARTH*COS(L0)*(F1-F0)
      if(abs(DELTAX).lt.1.e-8) DELTAX=1.e-8
      DELTAY=REARTH*(L1-L0)
      DXDY=DELTAX*DELTAY
      R=SQRT(DELTAX**2+DELTAY**2)
      alpha=atan(abs(DELTAY/DELTAX))*sign(1.,DXDY) +
     1      (1 - sign(1.,DXDY))*pi/2 + (1 - sign(1.,DELTAY))*pi/2
      if(alpha.ge.pi/4) then
        alpha = alpha - pi/4
      else
        alpha = alpha - pi/4 + 2*pi
      end if
      call verinterp(5,1,alphv,Rref18v,alpha,Rref18)
      call verinterp(5,1,alphv,Rref26v,alpha,Rref26)

      call verinterp(14,1,radm,angl,R,RANGL)

c      if(R.GT.RADM(14).OR.R.EQ.0) then
c        k=14
c      else
c        k=1
c        do while(R.GE.RADM(k+1))
c          k=k+1
c        end do
c      end if

C    CALCULATING WIND SPEED 

      if(Rref18.le.0.and.Rref26.le.0) then
        WND=SQRT(A7*B*DELP*EXP(-A7/R**B)/(ROA*R**B)+R**2*
     1  COR(I,J)**2/4.)-R*COR(I,J)/2.
        UTXa=UTX/2.
        UTYa=UTY/2.
      else
        WND=EXPWND(R,RMAX,Rref18,Rref26,WSMAX,WS18,WS26)
        UTXa=0.
        UTYa=0.
      end if
c      RANGL=ANGL(K)*PI/180.
      RANGL=RANGL*PI/180.
c      WF=WND
c      WR=-WND*TAN(RANGL)
      WF=WND*COS(RANGL)
      WR=-WND*SIN(RANGL)
      WX=WR*(DELTAX)/R-WF*(DELTAY)/R+UTXa
      WY=WF*(DELTAX)/R+WR*(DELTAY)/R+UTYa
      WM=SQRT(WX**2+WY**2)
cRMY USE NEW STRESS PARAMETERIZATION: CD3 FROM TUNG'S THESIS
cRMY      IF(WM.LT.10.) CD=1.14*1.E-3
cRMY      IF(WM.GE.10.) CD=(0.49+.065*WM)*1.E-3
cRMY      IF(WM.LE.35.) then
cRMY      WUSURF(I,J)=WUSURF(I,J)-CD*ROA*WM*WX/RHO_0
cRMY      else
cRMY      WUSURF(I,J)=WUSURF(I,J)-(3.3368+
cRMY     1      (WM-34.0449)**0.3)*WX/(WM*RHO_0)
cRMY      end if
      IF(WM.LT.12.5) then
        Z0W=0.0185/GRAV*(0.001*WM**2.+0.028*WM)**2.
        CD=0.4**2./(log(10./Z0W))**2.
      else
        CD1=3.58e-9*WM**3.-9.88e-7*WM**2.+7.81e-5*WM+7.9107e-4
        DCDL=1.12e-7*WM**2.+2.49e-5*WM-3.32e-4
        CD=CD1-DCDL
      end if
      WUSURF(I,J)=WUSURF(I,J)-CD*ROA*WM*WX/RHO_0

C  CALCULATING V-WIND STRESS FOR I,J POINT OF V-VELOCITY GRID
      F1=LONGMIN+(I-1)*(LONGMAX-LONGMIN)/(IM-1)
      F1=F1*2.*PI/360
      L1=LATMIN+(J-1)*(LATMAX-LATMIN)/((JM-0.5)-1)
      L1=L1*2.*PI/360
      DELTAX=REARTH*COS(L0)*(F1-F0)
c------------ falk 01-05-04 check DELTAX
      if(abs(DELTAX).lt.1.e-8) DELTAX=1.e-8
      DELTAY=REARTH*(L1-L0)
      DXDY=DELTAX*DELTAY
      R=SQRT(DELTAX**2+DELTAY**2)
      alpha=atan(abs(DELTAY/DELTAX))*sign(1.,DXDY)+
     1      (1 - sign(1.,DXDY))*pi/2 + (1 - sign(1.,DELTAY))*pi/2
      if(alpha.ge.pi/4) then
        alpha = alpha - pi/4
      else
        alpha = alpha - pi/4 + 2*pi
      end if
      call verinterp(5,1,alphv,Rref18v,alpha,Rref18)
      call verinterp(5,1,alphv,Rref26v,alpha,Rref26)

      call verinterp(14,1,radm,angl,R,RANGL)

C    CALCULATING WIND SPEED 

      if(Rref18.le.0.and.Rref26.le.0) then
        WND=SQRT(A7*B*DELP*EXP(-A7/R**B)/(ROA*R**B)+R**2*
     1  COR(I,J)**2/4.)-R*COR(I,J)/2.
        UTXa=UTX/2.
        UTYa=UTY/2.
      else
        WND=EXPWND(R,RMAX,Rref18,Rref26,WSMAX,WS18,WS26)
        UTXa=0.
        UTYa=0.
      end if
c      RANGL=ANGL(K)*PI/180.
      RANGL=RANGL*PI/180.
      WF=WND*COS(RANGL)
      WR=-WND*SIN(RANGL)
      WX=WR*(DELTAX)/R-WF*(DELTAY)/R+UTXa
      WY=WF*(DELTAX)/R+WR*(DELTAY)/R+UTYa
      WM=SQRT(WX**2+WY**2)
cRMY USE NEW STRESS PARAMETERIZATION: CD3 FROM TUNG'S THESIS
cRMY      IF(WM.LT.10.) CD=1.14*1.E-3
cRMY      IF(WM.GE.10.) CD=(0.49+.065*WM)*1.E-3
cRMY      IF(WM.le.35.) then
cRMY      WVSURF(I,J)=WVSURF(I,J)-CD*ROA*WM*WY/RHO_0
cRMY      else
cRMY      WVSURF(I,J)=WVSURF(I,J)-(3.3368+
cRMY     1      (WM-34.0449)**0.3)*WY/(WM*RHO_0)
cRMY      end if
      IF(WM.LT.12.5) then
        Z0W=0.0185/GRAV*(0.001*WM**2.+0.028*WM)**2.
        CD=0.4**2./(log(10./Z0W))**2.
      else
        CD1=3.58e-9*WM**3.-9.88e-7*WM**2.+7.81e-5*WM+7.9107e-4
        DCDL=1.12e-7*WM**2.+2.49e-5*WM-3.32e-4
        CD=CD1-DCDL
      end if
      WVSURF(I,J)=WVSURF(I,J)-CD*ROA*WM*WY/RHO_0

C  CALCULATING WIND STRESS FOR I,J POINT OF DEPTH GRID
      F1=LONGMIN+(I-1)*(LONGMAX-LONGMIN)/(IM-1)
      F1=F1*2.*PI/360.
      L1=LATMIN+(J-1)*(LATMAX-LATMIN)/(JM-1)
      L1=L1*2.*PI/360.
      DELTAX=REARTH*COS(L0)*(F1-F0)
      if(abs(DELTAX).lt.1.e-8) DELTAX=1.e-8
      DELTAY=REARTH*(L1-L0)
      DXDY=DELTAX*DELTAY
      R=SQRT(DELTAX**2+DELTAY**2)
      alpha=atan(abs(DELTAY/DELTAX))*sign(1.,DXDY)+
     1      (1 - sign(1.,DXDY))*pi/2 + (1 - sign(1.,DELTAY))*pi/2
      if(alpha.ge.pi/4) then
        alpha = alpha - pi/4
      else
        alpha = alpha - pi/4 + 2*pi
      end if
      call verinterp(5,1,alphv,Rref18v,alpha,Rref18)
      call verinterp(5,1,alphv,Rref26v,alpha,Rref26)

      call verinterp(14,1,radm,angl,R,RANGL)

C    CALCULATING WIND SPEED 

      if(Rref18.le.0.and.Rref26.le.0) then
        WND=SQRT(A7*B*DELP*EXP(-A7/R**B)/(ROA*R**B)+R**2*
     1  COR(I,J)**2/4.)-R*COR(I,J)/2.
        UTXa=UTX/2.
        UTYa=UTY/2.
      else
        WND=EXPWND(R,RMAX,Rref18,Rref26,WSMAX,WS18,WS26)
        UTXa=0.
        UTYa=0.
      end if
c      RANGL=ANGL(K)*PI/180.
      RANGL=RANGL*PI/180.
      WF=WND*COS(RANGL)
      WR=-WND*SIN(RANGL)
      WX=WR*(DELTAX)/R-WF*(DELTAY)/R+UTXa
      WY=WF*(DELTAX)/R+WR*(DELTAY)/R+UTYa
c      WINDX(i,j)=WX
c      WINDY(i,j)=WY
      WM=SQRT(WX**2+WY**2)
cRMY USE NEW STRESS PARAMETERIZATION: CD3 FROM TUNG'S THESIS
cRMY      IF(WM.LT.10.) CD=1.14*1.E-3
cRMY      IF(WM.GE.10.) CD=(0.49+.065*WM)*1.E-3
cRMY      IF(WM.GT.35.) THEN
cRMY      TMAX=3.3368+(WM-34.0449)**0.3
cRMY      TAUX(I,J)=TAUX(I,J)+TMAX*WX/WM
cRMY      TAUY(I,J)=TAUY(I,J)+TMAX*WY/WM
cRMY      ELSE
cRMY      TAUY(I,J)=TAUY(I,J)+CD*ROA*WM*WY
cRMY      TAUX(I,J)=TAUX(I,J)+CD*ROA*WM*WX
cRMY
cRMYC          WTSURF(I,J)=30.*WM/(4000*1024.)
cRMY
cRMY     END IF
      IF(WM.LT.12.5) then
        Z0W=0.0185/GRAV*(0.001*WM**2.+0.028*WM)**2.
        CD=0.4**2./(log(10./Z0W))**2.
      else
        CD1=3.58e-9*WM**3.-9.88e-7*WM**2.+7.81e-5*WM+7.9107e-4
        DCDL=1.12e-7*WM**2.+2.49e-5*WM-3.32e-4
        CD=CD1-DCDL
      end if
      TAUY(I,J)=TAUY(I,J)+CD*ROA*WM*WY
      TAUX(I,J)=TAUX(I,J)+CD*ROA*WM*WX

 351  CONTINUE
 350  CONTINUE

 230  FORMAT(151e16.7)
 231  FORMAT(10F8.3)
c
c------------------------- falk 01-05-04 add call avrtau
      call avrtau(x0,y0,0)
c
c----------- falk 09-12-05 change output
      if(MOD(IINT,24).EQ.0) 
     *print*,'Exiting WIND ...'
      RETURN
      END
