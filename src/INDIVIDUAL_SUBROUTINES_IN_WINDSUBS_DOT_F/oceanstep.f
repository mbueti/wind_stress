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
