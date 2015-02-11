      program fit_driv

C     f77 FITEXY-driver.f -L$PGPLOT_DIR -lpgplot -L/usr/X11R6 -lX11 bces_regress.f ran3.f bootspbec.f sort.f 

      implicit none
      integer choice,Mmax,regres
      real x(99),y(99),sigx(99),sigy(99),ylow,yhigh
      real epsilon,scat,top,astart,bstart,epsstart
      real chisq(400,1500),amin,bmin,chisqmin,Abscat
      real a,b,siga,sigb,denom
      real*8 OrthA,OrthB,OrthC,OrthD
       real OrthA4,OrthB4,OrthC4,OrthD4
      integer ndat,i,loop,epsdirect,k,l,i2,j2
      character thefile*16 
      common /DOG/ thefile,regres


      write(6,*) '1) Graham-Scott M-M    2) McConnell & Ma M-sigma'
      write(6,*) '3) Wandel (1999)       4) Savorgnan 2014'
      read (5,*) choice

      if (choice.eq.1) then 
         open (40,file='for-FITEXY.dat',status='old')
      elseif (choice.eq.2) then 
         open (40,file='McConnell-Ma-2013.dat',status='old')
      elseif (choice.eq.3) then 
         open (40,file='Wandel-1999.dat',status='old')
      elseif (choice.eq.4) then 
         open (40,file='Savorgnan_2014.dat',status='old')
      endif
      open (50, file='MassBH-Bersh.dat')
         read (40,*) ndat
         write (50,*) ndat
         Mmax=ndat
         do i=1,ndat
            read (40,*) X(i),sigx(i),Y(i),sigy(i)
            write (50,*) X(i),sigx(i),Y(i),sigy(i),0.0
         enddo
      close (unit=40)
      close (unit=50)

         if (choice.eq.1) then
            astart=5.3
            bstart=1.0
            epsstart=0.0
         elseif (choice.eq.2) then
            astart=8.2
            bstart=5.5 
            epsstart=0.0
         elseif (choice.eq.3) then
            astart=8.2
            bstart=1.0
            epsstart=0.0
         elseif (choice.eq.4) then
            astart=8.0
            bstart=5.0 
            epsstart=0.0
         endif


      do loop=1,500
         epsilon=epsstart + 0.001*loop

          chisqmin=1200.  ! some arbitrary high value
          do i2=1,200
            a=astart+i2/1000.
            do j2=1,1500
               b=bstart+j2/1000.

           chisq(i2,j2)=0.0
       do l=1,ndat
          top=y(l) - ( a + b*x(l) )
c          denom=(sigy(l)*sigy(l) + b*b*sigx(l)*sigx(l) +      ! XXX residual in x
c     :           b*b*epsilon*epsilon)                         ! XXX residual in x   
          denom=(sigy(l)*sigy(l) + b*b*sigx(l)*sigx(l) +     ! XXX residual in y
     :           epsilon*epsilon)                            ! XXX residual in y
          if (choice.eq.2) then 
          if (l.eq.16.OR.l.eq.17) then 
             denom=denom*sqrt(2.0)     ! NGC 1399 
          endif
          endif         
          chisq(i2,j2) = chisq(i2,j2) + top*top/denom         
       enddo

       if (chisq(i2,j2).le.(chisqmin)) then
         chisqmin=chisq(i2,j2)
         amin=a
         bmin=b
       endif
          enddo  ! i2 loop
         enddo   !  j2 loop

      if (chisqmin.le.(1.*(ndat-2))) then
         goto 20
      endif

      enddo  !  return to the epsilon loop.
20    continue

       write(6,*) 'FINAL: a,b and intrinsic scatter:',
     :             amin,bmin,epsilon

       abscat=0.0
       do l=1,ndat
         abscat=abscat+ (y(l)- (amin+bmin*x(l)))**2
       enddo
       abscat=sqrt(abscat/ndat)
       write(6,*) 'Abs. scatter in log(M_bh) = ',abscat
       write(6,*) ' '

C    ===========================================

      END
