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
