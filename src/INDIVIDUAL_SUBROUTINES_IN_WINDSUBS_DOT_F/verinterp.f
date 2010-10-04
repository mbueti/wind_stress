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
