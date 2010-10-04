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
