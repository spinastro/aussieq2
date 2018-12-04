C >>> YCK
C Fall 2000 Written for the standard set of the YY
C Dec 2002 Modfied for the introduction of pseudo-eep
C Dec 2015 Modified for the Yapsi isochrone set
C Jan 2017 put a switch, Lremap, for re-mapping EEPs or not
C <<< YCK
C WARNING     THIS ROUTINE ASSUMES numage=41 AT EACH ISOCHRONE
C             and the set of the 41 ages should be the same.
C numfile = the number of the isochrone files + 1, for each Y set
C
C The names of the two input files are fixed to be 'YaPsI.nml' and 'YaPsI.age'
C
C This program reads data from YaPsI.nml file.
C when the read-in Age=0., the age is read-in from YaPsI.age 
C The name of the source isochrone tables (yapsi*.*)
C should be assigned at YYiso(*).
CCCCCCCCCC   YaPsI.nml   CCCCCCCCCCCCC
C $INPUT
C targtY=0.27    ! Y  between 0.25 and 0.37
C targtFe=-0.03  ! [Fe/H] between -1.5 and +0.3
C Age= 9.        ! age in Gyr: if negative, interactive input
C                !             if 0., read 'YaPsI.age'
C YYout='yy.out' ! output file name: if ' ', interactive input
C                ! if 'each', the output file name will be asked 
C                !   for each age, and stored in seperate files.
C Lremap=.FALSE. ! 
C                ! if .TRUE., re-map the EEPs to set the MSTO at the same grid
C $END
C $ISET1         ! the set for Y=0.25
C NYYiso=5
C YYiso(1)='LCB/yapsi_l_X0p749455_Z0p000545.dat'! the file for [Fe/H]=-1.5
C YYiso(2)='LCB/yapsi_l_X0p748279_Z0p001721.dat'! the file for [Fe/H]=-1.0
C YYiso(3)='LCB/yapsi_l_X0p744584_Z0p005416.dat'! the file for [Fe/H]=-0.5
C YYiso(4)='LCB/yapsi_l_X0p733138_Z0p016862.dat'! the file for [Fe/H]= 0.0
C YYiso(5)='LCB/yapsi_l_X0p717092_Z0p032908.dat'! the file for [Fe/H]=+0.3
C $END
C $ISET2         ! the set for Y=0.28
C NYYiso=5
C YYiso(1)='LCB/yapsi_l_X0p719477_Z0p000523.dat'! the file for [Fe/H]=-1.5
C YYiso(2)='LCB/yapsi_l_X0p718348_Z0p001652.dat'! the file for [Fe/H]=-1.0
C YYiso(3)='LCB/yapsi_l_X0p714801_Z0p005199.dat'! the file for [Fe/H]=-0.5
C YYiso(4)='LCB/yapsi_l_X0p703812_Z0p016188.dat'! the file for [Fe/H]= 0.0
C YYiso(5)='LCB/yapsi_l_X0p688408_Z0p031592.dat'! the file for [Fe/H]=+0.3
C $END
C $ISET3         ! the set for Y=0.31
C NYYiso=5
C YYiso(1)='LCB/yapsi_l_X0p689499_Z0p000501.dat'! the file for [Fe/H]=-1.5
C YYiso(2)='LCB/yapsi_l_X0p688417_Z0p001583.dat'! the file for [Fe/H]=-1.0
C YYiso(3)='LCB/yapsi_l_X0p685018_Z0p004982.dat'! the file for [Fe/H]=-0.5
C YYiso(4)='LCB/yapsi_l_X0p674487_Z0p015513.dat'! the file for [Fe/H]= 0.0
C YYiso(5)='LCB/yapsi_l_X0p659725_Z0p030275.dat'! the file for [Fe/H]=+0.3
C $END
C $ISET4         ! the set for Y=0.34
C NYYiso=5
C YYiso(1)='LCB/yapsi_l_X0p659520_Z0p000480.dat'! the file for [Fe/H]=-1.5
C YYiso(2)='LCB/yapsi_l_X0p658485_Z0p001515.dat'! the file for [Fe/H]=-1.0
C YYiso(3)='LCB/yapsi_l_X0p655234_Z0p004766.dat'! the file for [Fe/H]=-0.5
C YYiso(4)='LCB/yapsi_l_X0p645161_Z0p014839.dat'! the file for [Fe/H]= 0.0
C YYiso(5)='LCB/yapsi_l_X0p631041_Z0p028959.dat'! the file for [Fe/H]=+0.3
C $END
C $ISET5         ! the set for Y=0.37
C NYYiso=5
C YYiso(1)='LCB/yapsi_l_X0p629542_Z0p000458.dat'! the file for [Fe/H]=-1.5
C YYiso(2)='LCB/yapsi_l_X0p628554_Z0p001446.dat'! the file for [Fe/H]=-1.0
C YYiso(3)='LCB/yapsi_l_X0p625451_Z0p004549.dat'! the file for [Fe/H]=-0.5
C YYiso(4)='LCB/yapsi_l_X0p615836_Z0p014164.dat'! the file for [Fe/H]= 0.0
C YYiso(5)='LCB/yapsi_l_X0p602357_Z0p027643.dat'! the file for [Fe/H]=+0.3
C $END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Warning!
C (Z/X)_sun value used for the conversion from [Fe/H] to Z
C      parameter (ZXsun=0.023D0 )
C Warning!!
C When re-mapped, the MSTO is at kipgrd
C      parameter (kipgrd=130)

      program YYinbetween
      IMPLICIT REAL*8 (A-H,O-Z)
      integer UNIN,UNOUT
      parameter (IYY=4)
      parameter (UNIN=7)
      parameter (UNOUT=9)
      parameter (numfile=6)
      parameter (numYset=5)
      parameter (numage=41)
      parameter (numpoint=200)
      parameter (numparam=12)
C when the differencce between a target value and grid is less than epsval,
C use the grid value.
      parameter (epsval=1.1d-5)
C
      character*50 YYiso(numfile),YYout,YYout2
      data YYiso/numfile*' '/
      data YYout/' '/
      real*8 YYage(numage,numfile),YYzarr(numfile),YYyarr(numYset)
      real*8 YYarray(numparam,numpoint,numage,numfile)
      integer iYYmass(numage,numfile),ofile
      common /workspace/YYarray,YYage,YYzarr,YYyarr,iYYmass
      data NYYiso,targtFe,Age/1,0.0,12.0/
C YYzarr is the [Fe/H] grid, 
      data YYzarr/-1.5d0,-1.0d0,-0.5d0,0.0d0,0.3d0,0.0d0/
C YYyarr is the Y grid, 
C      data YYyarr/0.25d0,0.265d0,0.278d0,0.295d0,0.31d0/
      data YYyarr/0.25d0,0.28d0,0.31d0,0.34d0,0.37d0/
c a switch for re-mapping EEPs
      LOGICAL*4 Lremap
      DATA Lremap/.FALSE./
      common /switch/Lremap
c
      NAMELIST /INPUT/targtY,targtFe,Age,YYout,Lremap
      NAMELIST /ISET1/NYYiso,YYiso
      NAMELIST /ISET2/NYYiso,YYiso
      NAMELIST /ISET3/NYYiso,YYiso
      NAMELIST /ISET4/NYYiso,YYiso
      NAMELIST /ISET5/NYYiso,YYiso
C
      open(UNIT=IYY,file='YaPSI.nml',status='old')
      READ(UNIT=IYY, NML=INPUT)
      close(IYY)
         FeH0=targtFe
c       write(6,NML=INPUT)
c       call  feh2z(targtFe,targtY,targtX,targtZ)
c      write(6,831) targtX,targtY,targtZ,targtFe,Age
c  831 format(18H# interpolated for,
c     +   ' X= ',f7.5,' Y= ',f7.5,' Z=',f7.5,' [Fe/H] = ',f6.3,
c     +   ' Age (Gyr)=',f6.3)

      nforyi=0
      if(YYout.ne."each")then
      nforyi=1
      YYout2=YYout
      if(YYout2.eq." ")then
      write(6,'(1x,A,$)')' output file name ? '
      read(5,'(a)')YYout2
      endif
      open(UNIT=UNOUT,file=YYout2,status='new')
      endif

      if(Age.gt.0.0d0)then
        agenew=Age
      elseif(Age.eq.0.0d0)then
        open(UNIT=UNIN,file='YaPSI.age',status='old')
        read(UNIN,*,end=900,err=900)agenew
      else
        write(6,'(A,$)')' New Age (input a negative value to stop) '
        read(5,*,end=900,err=900)agenew
        if(agenew.le.0.0e0)go to 900
      endif
      goto 401
  600 continue
      if(Age.lt.0.0d0)then
        write(6,'(A,$)')' New Age (input a negative value to stop) '
        read(5,*,end=900,err=900)agenew
        if(agenew.le.0.0e0)go to 900
      else
        read(UNIN,*,end=900,err=900)agenew
      endif
  401 continue
C
      if(YYout.eq."each")then
      write(6,'(1x,A,$)')' output file name ? '
      read(5,'(a)')YYout2
      open(UNIT=UNOUT,file=YYout2,status='new')
      endif
C
        iyfile=numYset*0.5
       call findz(YYyarr,numYset,targtY,iyfile)
cD     write(6,*)YYyarr(iyfile),targtY,iyfile


      if(DABS((YYyarr(iyfile)-targtY)).lt.epsval)then

c      write(6,*)' only 1 Y set is used '
      open(UNIT=IYY,file='YaPSI.nml',status='old')
        if(iyfile.eq.1) READ(UNIT=IYY, NML=ISET1)
        if(iyfile.eq.2) READ(UNIT=IYY, NML=ISET2)
        if(iyfile.eq.3) READ(UNIT=IYY, NML=ISET3)
        if(iyfile.eq.4) READ(UNIT=IYY, NML=ISET4)
        if(iyfile.eq.5) READ(UNIT=IYY, NML=ISET5)
      close(IYY)
cD     write(6,'(11H Y set ISET,i1,8H read in)')iyfile
         Nalp=1
         call Y_each(targtFe,YYiso,agenew,NYYiso,Nalp,iyfile)
         call fin_write(Nalp,nforyi)

      else
cD     write(6,*)' 3 Y sets are used'
       iyfile=iyfile-1
       if(iyfile.le.1)iyfile=1
       if(iyfile.ge.(numYset-1))iyfile=iyfile-1

      Nalp=0
       do iyset=iyfile, iyfile+2
        open(UNIT=IYY,file='YaPSI.nml',status='old')
          if(iyset.eq.1) READ(UNIT=IYY, NML=ISET1)
          if(iyset.eq.2) READ(UNIT=IYY, NML=ISET2)
          if(iyset.eq.3) READ(UNIT=IYY, NML=ISET3)
          if(iyset.eq.4) READ(UNIT=IYY, NML=ISET4)
          if(iyset.eq.5) READ(UNIT=IYY, NML=ISET5)
        close(IYY)
cD       write(6,'(11H Y set ISET,i1,8H read in)')iyset
        Nalp=Nalp+1
        targtFe=FeH0
        call Y_each(targtFe,YYiso,agenew,NYYiso,Nalp,iyset)
       enddo
         
         targtFe=FeH0
         call intpol_Y(targtFe,targtY,iyfile)
         Nalp=1
         call fin_write(Nalp,nforyi)

      endif
C

      if(YYout.eq."each")then
      close(UNOUT)
      write(6,*)' File created ==> ',YYout2
      endif

      if(Age.gt.0.0d0)goto 900
      go to 600
  900 continue
      if(Age.eq.0.0d0)close(UNIN)
      if(YYout.ne."each")then
      close(UNOUT)
      write(6,*)' File created ==> ',YYout2
      endif
      stop 'All done...'
      end

      subroutine Y_each(targtFe,YYiso,Age,NYYiso,Nalp,iyfile)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer UNIN,UNOUT
      parameter (IYY=14)
      parameter (UNIN=7)
      parameter (UNOUT=9)
      parameter (numfile=6)
      parameter (numYset=5)
      parameter (numage=41)
      parameter (numpoint=200)
      parameter (numparam=12)
      parameter (epsval=1.1d-5)
C
      character*50 YYiso(numfile)
      real*8 YYage(numage,numfile),YYzarr(numfile),YYyarr(numYset)
      real*8 YYarray(numparam,numpoint,numage,numfile)
      integer iYYmass(numage,numfile),ofile
c a switch for re-mapping EEPs or not
      LOGICAL*4 Lremap
      common /switch/Lremap
      common /workspace/YYarray,YYage,YYzarr,YYyarr,iYYmass

C Read tables
      do 10 ifile=1,NYYiso
      call yyread(YYiso(ifile), ifile,iyfile)
C
      if ( Lremap )call resetGRID(ifile)
C
   10 continue
      if(NYYiso.le.1)go to 100
C Z interpolation
      Znew=targtFe
      ifile=NYYiso*0.5
      call findz(YYzarr,NYYiso,Znew,ifile)
C in case
      if( DABS((YYzarr(ifile)-targtFe)).lt.epsval)then
       goto 100
      endif
C
      ifile=ifile-1
      if(ifile.ge.(NYYiso-3))ifile=NYYiso-3
      if(ifile.le.1)ifile=1
CD     write(6,'(A)')'# [Fe/H] interpolation    V'
CD     write(6,'(5f10.3)')YYzarr(ifile),YYzarr(ifile+1),Znew,
CD    +          YYzarr(ifile+2),YYzarr(ifile+3)
      ofile=NYYiso+1

      YYzarr(ofile)=Znew

      call intpol_z(ifile,ofile)
C the output stored in (NYYiso+1)
      ifile=ofile
C
  100 continue
C Age interpolation
      if(NYYiso.le.1)then
         ifile=1
      endif
CD     write(6,*)ifile,' find age index',Age
CD     write(6,*) (YYage(i1,ifile),i1=1,numage)

        agenew=Age
      call findz(YYage(1,ifile),numage,agenew,iage)
CD     write(6,*)iage,' found'
C for linear
      if(iage.ge.(numage-1))iage=numage-1
      if(iage.le.1)iage=1
cD     write(6,'(A)')'# Age interpolation       V'
cD     write(6,'(10x,5f10.5)')YYage(iage,ifile),agenew,
cD    +                       YYage(iage+1,ifile)
C
      call intpol_age(ifile,iage,agenew,Nalp,iyfile)
  900 continue
      return
      end
C
      subroutine yyread(YYisofl,ifile,iyfile)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer UNIN,UNOUT
      parameter (IYY=24)
      parameter (UNIN=7)
      parameter (UNOUT=9)
      parameter (numfile=6)
      parameter (numYset=5)
      parameter (numage=41)
      parameter (numpoint=200)
      parameter (numparam=12)
      parameter (epsval=1.1d-5)
      character*50 YYisofl
      character*138 header,chfmt(1)
c
      common /head/YYos,YYalp,YYfeh,YYalpfe,nump,header,chfmt
      real*8 YYarray(numparam,numpoint,numage,numfile)
      real*8 YYage(numage,numfile),YYz(numfile),YYy(numYset)
      integer iYYmass(numage,numfile),ofile
      common /workspace/YYarray,YYage,YYz,YYy,iYYmass
c for nomalization
      real*8 eep(numpoint),xeep(numpoint),dteep,tlkip,blkip,slope
      real*8 temar(numparam,numpoint),xxeep
cD write(6,*)YYisofl,ifile,iyfile
      open(IYY,file=YYisofl,status='old')
  810 format(2x,f8.6,3x,f8.6,4x,f4.2,6x,f8.6,8x,f9.6,12x,f5.2)
!      chfmt(1)='(f10.7,3f8.4,5f7.3)'
! new format >>>
      chfmt(1)='(f9.4,f10.7,3f8.4,8f7.3)'
! <<<
      nump=numparam
      nmass=numpoint
      read(IYY,'(A)')header
c
c >>> to check the consistency
!      read(header,'(9x,i5,11x,i5)')iiiix,iiiiz
      read(header,'(12x,i6,11x,i6)')iiiix,iiiiz
      xisoch=float(iiiix)*1.0e-6
      zisoch=float(iiiiz)*1.0e-6
      yisoch=1.0e0-xisoch-zisoch
      call feh2z( YYz(ifile), YYy(iyfile),Xgrid,Zgrid)
      if( abs(zisoch-Zgrid) .gt. epsval .or.
     +     abs(yisoch-YYy(iyfile)). gt. epsval)then
       write(6,*)' Warning...  isochrone file mismatch '
       write(6,'(1x, 3f10.6)')
     +    xisoch,yisoch,zisoch,Xgrid,YYy(iyfile),Zgrid 
      endif
c <<<

      read(IYY,'(A)')header
      iage=0
 901  continue
      read(IYY,'(12x,f6.3)',end=910)age
      iage=iage+1
      YYage(iage,ifile)=age
      iYYmass(iage,ifile)=nmass
      do 100 imass=1,nmass
C#       m,t,l,g,mv,ub,bv,vr,vi
C#  Age    M/Msun     logT  logL/Ls   logg    Mv     
C   U-B    B-V    V-R    V-I    J-K    H-K    V-K
      read(IYY,chfmt,err=901)Agetmp,
     +      (YYarray(i,imass,iage,ifile),i=1,nump)
      if(Agetmp.ne.age)then
       write(6,*)' Age is wrong in the file'
       write(6,*)YYisofl,ifile,iyfile,Age
       close(IYY)
       stop
      endif
  100 continue

      read(IYY,'()')
      goto 901
 910  continue
      close(IYY)
      nage=iage
      if(nage .ne. numage)then
      write(6,*)' the number of isochrone is differ from numage=41'
      endif
      return
      end
C

Cc a switch for re-mapping EEPs or not
C      LOGICAL*4 Lremap
C      common /switch/Lremap
CC
C      if ( Lremap )call resetGRID(ifile)
CC
      subroutine resetGRID(ifile)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer UNIN,UNOUT
      parameter (IYY=24)
      parameter (UNIN=7)
      parameter (UNOUT=9)
      parameter (numfile=6)
      parameter (numYset=5)
      parameter (numage=41)
      parameter (numpoint=200)
      parameter (numparam=12)
      parameter (epsval=1.1d-5)
      parameter (kipgrd=130)
      character*50 YYisofl
      character*138 header,chfmt(1)
c
      common /head/YYos,YYalp,YYfeh,YYalpfe,nump,header,chfmt
      real*8 YYarray(numparam,numpoint,numage,numfile)
      real*8 YYage(numage,numfile),YYz(numfile),YYy(numYset)
      integer iYYmass(numage,numfile),ofile
      common /workspace/YYarray,YYage,YYz,YYy,iYYmass
c for nomalization
      integer meep,nmass,igrd
      real*8 eep(numpoint),xeep(numpoint),dteep,tlkip,blkip,TOpoint
      real*8 temar(numparam,numpoint),xxeep
      pollin(x1,y1,x2,y2,x)=
     +  (x-x2)*y1/(x1-x2) +(x-x1)*y2/(x2-x1)

      nump=numparam
      nmass=numpoint
c to re-map TO at kipgrd
c When re-mapped, the MSTO will be at the grid number of kipgrd
c        kipgrd=130

c Because the younger ones do not have TO
c To use the older one's TO for younger ones
      iage=numage+1
 901  continue
      if(iage.le.1)goto 910
      iage=iage-1
c--------------yi----------------use length
       imass=1
        tlkip=YYarray(2,imass,iage,ifile)
        blkip=YYarray(3,imass,iage,ifile)
        eep(imass)=1.0d0
      do 100 imass=2,nmass
        eep(imass)=eep(imass-1)
     +           + dsqrt(1.0d2*(YYarray(2,imass,iage,ifile)-tlkip)
     +                  *(YYarray(2,imass,iage,ifile)-tlkip)
     +                  +(YYarray(3,imass,iage,ifile)-blkip)
     +                  *(YYarray(3,imass,iage,ifile)-blkip))
        tlkip=YYarray(2,imass,iage,ifile)
        blkip=YYarray(3,imass,iage,ifile)
c--------------yi----------------
c--------------kim---------------use mass
c      do 100 imass=1,nmass
c        eep(imass)=YYarray(1,imass,iage,ifile)
c--------------kim---------------
        xeep(imass)=YYarray(2,imass,iage,ifile)
  100 continue

C >>> YCK modified for the eep
c to find the turnoff mass to utilize it as an anchor point
         call eepset(eep,xeep,nmass,igrd,TOpoint) 
C        if(igrd.eq.nmass)then
c the isochrone does not have a TOpoint.  
c use the previous TO
c        else
c do have a TO
c        TOkip=TOpoint
C        endif
        if(igrd.ne.nmass) TOkip=TOpoint
c EEP spacing before the TO
        delbf=(TOkip-eep(1))/dfloat(kipgrd-1)
c EEP spacing after the TO
        delaf=(eep(nmass)-TOkip)/dfloat(nmass-kipgrd)

         do im=1,kipgrd
         xeep(im)= dfloat(im-1)*delbf  +eep(1)
         enddo
         do im=kipgrd+1,nmass
         xeep(im)= dfloat(im-kipgrd)*delaf  +TOkip
C         xeep(im)= dfloat(im-kipgrd)*delaf  +xeep(kipgrd)
         enddo

c
       meep=1
       do 2003 im=1,nmass
       xxeep=xeep(im)
       call findz(eep,nmass,xxeep,meep)
       if(meep.ge.(nmass-1))meep=nmass-1
       if(meep.le.1)meep=1
       do  20031 ii=1,nump
       temar(ii,im)= pollin(
     +  eep(meep),YYarray(ii,meep,iage,ifile),
     +  eep(meep+1),YYarray(ii,meep+1,iage,ifile),
     +  xxeep )
20031  continue
 2003  continue
c normalized
       do 20032 im=1,nmass
       do 20033 ii=1,nump
        YYarray(ii,im,iage,ifile)=temar(ii,im)
20033  continue
20032  continue
cc <<< YCK  Done!

      goto 901
 910  continue
c Re-mapping of EEPs done!
      return
      end
C
C
      subroutine findz(AX,NX,X,M)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AX(NX)
C FIND THE 'M'
      IF(M.LT.1.OR.M.GT.NX)M=1
      KIP=M
CCC      DO 21 K=1,NX
        IF(X.LT.AX(KIP))THEN
           DO 211 IYC=KIP-1,1,-1
             IF(AX(IYC).LE.X)THEN
               KIP=IYC
               GOTO 213
             ENDIF
  211      CONTINUE
c           WRITE(6,*)' SMALLER THAN THE RANGE : extrapolation '
C           WRITE(ISHORT,*)' ERROR findex: SMALLER THAN THE RANGE '
C           WRITE(6,*) 1,AX(1),X
           KIP=1
        ELSE
           DO 212 IYC=KIP,NX-1
             IF(AX(IYC+1).GT.X)THEN
               KIP=IYC
               GOTO 213
             ENDIF
  212      CONTINUE
               KIP = NX 
c         WRITE(6,*)' LARGER THAN THE RANGE : extrapolation '
C         WRITE(ISHORT,*)' ERROR findex: LARGER THAN THE RANGE '
C         WRITE(6,*)NX,AX(NX),X
        ENDIF
  213   M=KIP
      RETURN
      END
C
      subroutine intpol_z(ifile,ofile)
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter (numfile=6)
      parameter (numYset=5)
      parameter (numage=41)
      parameter (numpoint=200)
      parameter (numparam=12)
      character*50 YYiso(numfile)
      real*8 YYarray(numparam,numpoint,numage,numfile)
      real*8 YYage(numage,numfile),YYzarr(numfile),YYyarr(numYset)
      integer iYYmass(numage,numfile),ofile
      character*138 header,chfmt(1)
      common /head/YYos,YYalp,YYfeh,YYalpfe,nump,header,chfmt
      common /workspace/YYarray,YYage,YYzarr,YYyarr,iYYmass
      cubeint(x1,y1,x2,y2,x3,y3,x4,y4,x)=
     +  (x-x2)*(x-x3)*(x-x4)*y1/((x1-x2)*(x1-x3)*(x1-x4))
     + +(x-x1)*(x-x3)*(x-x4)*y2/((x2-x1)*(x2-x3)*(x2-x4))
     + +(x-x1)*(x-x2)*(x-x4)*y3/((x3-x1)*(x3-x2)*(x3-x4))
     + +(x-x1)*(x-x2)*(x-x3)*y4/((x4-x1)*(x4-x2)*(x4-x3))

      do 300 iage=1,numage
        iYYm=iYYmass(iage,ifile+1)
       if((YYage(iage,ifile).ne.YYage(iage,ifile+1)) .or.
     +    (YYage(iage,ifile).ne.YYage(iage,ifile+2)) .or.
     +    (YYage(iage,ifile).ne.YYage(iage,ifile+3)))then
        write(6,*)'Warning. Need age interpolation'
       write(6,*)YYage(iage,ifile),YYage(iage,ifile+1),
     +          YYage(iage,ifile+2),YYage(iage,ifile+3)
       stop 'error'
       endif
      YYage(iage,ofile)=YYage(iage,ifile)
      iYYmass(iage,ofile)=iYYm
      do 200 ipoint=1,iYYm
      do 100 iparam=1,nump
      YYarray(iparam,ipoint,iage,ofile)=
     +cubeint(YYzarr(ifile),YYarray(iparam,ipoint,iage,ifile),
     +        YYzarr(ifile+1),YYarray(iparam,ipoint,iage,ifile+1),
     +        YYzarr(ifile+2),YYarray(iparam,ipoint,iage,ifile+2),
     +        YYzarr(ifile+3),YYarray(iparam,ipoint,iage,ifile+3),
     +        YYzarr(ofile))
  100 continue
  200 continue
  300 continue

      return
      end
C
      subroutine intpol_age(ifile,iage,agenew,Nalp,iyfile)
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter (numfile=6)
      parameter (numYset=5)
      parameter (numage=41)
      parameter (numpoint=200)
      parameter (numparam=12)
      real*8 YYnew(numparam)
      character*50 YYiso(numfile)
      real*8 YYarray(numparam,numpoint,numage,numfile)
      real*8 YYage(numage,numfile),YYzarr(numfile),YYyarr(numYset)
      integer iYYmass(numage,numfile),ofile
      character*138 header,chfmt(1)
      common /head/YYos,YYalp,YYfeh,YYalpfe,nump,header,chfmt
      common /workspace/YYarray,YYage,YYzarr,YYyarr,iYYmass
      common /final/AYYfin(8,numYset),YYfin(numparam,numpoint,numYset)

      pollin(x1,y1,x2,y2,x)=
     +  (x-x2)*y1/(x1-x2) +(x-x1)*y2/(x2-x1)

      iYYm=iYYmass(iage+1,ifile)
      if(iYYmass(iage,ifile).lt.iYYm)iage=iage+1
c
       AYYfin(1,Nalp)=YYzarr(ifile)
       AYYfin(2,Nalp)=YYyarr(iyfile)
       AYYfin(7,Nalp)=agenew
       AYYfin(8,Nalp)=iYYm
c
      do 200 ipoint=1,iYYm
      do 100 iparam=1,nump
C output
      YYfin(iparam,ipoint,Nalp)=
     +pollin(YYage(iage,ifile),YYarray(iparam,ipoint,iage,ifile),
     +        YYage(iage+1,ifile),YYarray(iparam,ipoint,iage+1,ifile),
     +        agenew)
  100 continue
  200 continue

CD      iparam=2
CD      do ipoint=1,iYYm
CD       write(6,'(1x,6f10.5)')
CD     + 0.5*(YYarray(iparam,ipoint,iage,ifile)+
CD     +  YYarray(iparam,ipoint,iage+1,ifile)),
CD     + 0.5*(YYarray(iparam+1,ipoint,iage,ifile)+
CD     +  YYarray(iparam+1,ipoint,iage+1,ifile)),
CD     +  YYfin(iparam,ipoint,Nalp),
CD     +  YYfin(iparam+1,ipoint,Nalp)
CD      enddo

      return
      end

      subroutine feh2z(FeH,Y,X,Z)
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter (ZXsun=0.0230D0 )
c Compute X,Y,Z when the input is in [Fe/H]
      FeH0=FeH
      X=(1.0d0-Y)/(1.0d0+ZXsun*(10.0d0**FeH0))
      Z=1.0d0-Y-X
c      write(6,9)X,Y,Z,FeH0
c    9 format(1x,'New mixture',
c     +       ' X= ',f7.5,' Y= ',f7.5,' Z=',f7.5,' [Fe/H]= ',f6.3)
      return
      end
c
      subroutine fin_write(Nalp,nforyi)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer UNIN,UNOUT
      parameter (UNOUT=9)
      parameter (numfile=6)
      parameter (numYset=5)
      parameter (numage=41)
      parameter (numpoint=200)
      parameter (numparam=12)
      parameter (epsval=1.1d-5)
      real*8 YYnew(numparam)
      character*50 YYiso(numfile)
      real*8 YYarray(numparam,numpoint,numage,numfile)
      real*8 YYage(numage,numfile),YYzarr(numfile),YYyarr(numYset)
      integer iYYmass(numage,numfile),ofile
      real*8 mass(numpoint), logt(numpoint),logl(numpoint)
      real*8 massup, massdo,num(3,numpoint)
      character*138 header,chfmt(1)
      common /head/YYos,YYalp,YYfeh,YYalpfe,nump,header,chfmt
      common /workspace/YYarray,YYage,YYzarr,YYyarr,iYYmass
      common /final/AYYfin(8,numYset),YYfin(numparam,numpoint,numYset)
      data icounter/0/
      save icounter
      icounter=icounter+1


CCCCCCCCCCCCCCCC  
       iYYm=AYYfin(8,Nalp)
       targtFe=AYYfin(1,Nalp)
       targtY=AYYfin(2,Nalp)
       call  feh2z(targtFe,targtY,targtX,targtZ)
C screen out
      write(6,9)targtX,targtY,targtZ,targtFe,AYYfin(7,Nalp)
    9 format(1x,'New mixture',
     +       ' X= ',f7.5,' Y= ',f7.5,' Z= ',f7.5,' [Fe/H]= ',f6.3,
     +       ' Age= ',f7.4,' Gyr')

      if(nforyi.le.0)then
      write(UNOUT,830) targtX,targtY,targtZ,targtFe
  830 format(19H#  interpolated for,
     +   ' X= ',f7.5,' Y= ',f7.5,' Z=',f7.5,' [Fe/H]= ',f6.3)
      write(UNOUT,'(1h#,A97)')header
c      write(UNOUT,'(1h#,A79,/)')header
      write(UNOUT,840)AYYfin(7,Nalp)
  840 format(12H#  age(Gyr)=,f6.3)
      else
C-------
      if(icounter.le.1)then
      write(UNOUT,831) targtX,targtY,targtZ,targtFe
  831 format(19H#  interpolated for,
     +   ' X= ',f7.5,' Y= ',f7.5,' Z=',f7.5,' [Fe/H]= ',f6.3)
      write(UNOUT,'(A97)')header
c      write(UNOUT,'(A79,/)')header
      else
      write(UNOUT,'(1h )')
      endif
      write(UNOUT,841)AYYfin(7,Nalp)
  841 format(12H#  age(Gyr)=,f6.3)
C------
      endif
CCCCCCCCCCCCCCCC
      do 3000 ipoint=1,iYYm
       mass(ipoint)=YYfin(1,ipoint,Nalp)
       logt(ipoint)=YYfin(2,ipoint,Nalp)
       logl(ipoint)=YYfin(3,ipoint,Nalp)
c check the logg value
       complogg = dlog10(mass(ipoint))-logl(ipoint)
     +          +4.0d0*logt(ipoint) -10.612D0
c      if(DABS((YYfin(4,ipoint,Nalp)-complogg))
c     +    .gt.epsval)then
c
cc       write(6,*)' the gravity has been recomputed'
cc       write(6,*)YYfin(4,ipoint,Nalp),complogg,
cc     + DABS((YYfin(4,ipoint,Nalp)-complogg))
       YYfin(4,ipoint,Nalp)=complogg
c      endif

 3000 continue
c The output
      do 200 ipoint=1,iYYm
c

      write(UNOUT,chfmt)AYYfin(7,Nalp),(YYfin(ii,ipoint,Nalp),ii=1,nump)
  200 continue
c <<< Apr 2003
      return
      end
      
      subroutine intpol_Y(targtFe,targtY,iyfile)
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter (numfile=6)
      parameter (numYset=5)
      parameter (numage=41)
      parameter (numpoint=200)
      parameter (numparam=12)

      real*8 YYnew(numparam)
      character*50 YYiso(numfile)
      real*8 YYarray(numparam,numpoint,numage,numfile)
      real*8 YYage(numage,numfile),YYzarr(numfile),YYyarr(numYset)
      integer iYYmass(numage,numfile),ofile
      character*138 header,chfmt(1)
      common /head/YYos,YYalp,YYfeh,YYalpfe,nump,header,chfmt
      common /workspace/YYarray,YYage,YYzarr,YYyarr,iYYmass
      common /final/AYYfin(8,numYset),YYfin(numparam,numpoint,numYset)

      quad(x1,y1,x2,y2,x3,y3,x)=y1*(x2-x)*(x3-x)/((x2-x1)*(x3-x1))
     +                         +y2*(x1-x)*(x3-x)/((x1-x2)*(x3-x2))
     +                         +y3*(x1-x)*(x2-x)/((x1-x3)*(x2-x3))

cD     write(6,'(A)')'# Y interpolation         V '
cD     write(6,'(1x,4f10.5)')YYyarr(iyfile), YYyarr(iyfile+1),
cD    +                     targtY, YYyarr(iyfile+2)

      AYYfin(1,1)=targtFe
      AYYfin(2,1)=targtY
      iYYm=AYYfin(8,1)
c
      do 200 ipoint=1,iYYm
      do 100 iparam=1,nump
      YYfin(iparam,ipoint,1)=
     +       quad(YYyarr(iyfile)    ,YYfin(iparam,ipoint,1),
     +            YYyarr(iyfile+1)  ,YYfin(iparam,ipoint,2),
     +            YYyarr(iyfile+2)  ,YYfin(iparam,ipoint,3),targtY)
  100 continue
  200 continue
      return
      end
C
      subroutine eepset(x,y,nstep,igrd,tomass)
c      parameter (nstep=200)
      real*8 x(nstep),y(nstep),xp,xpmean,ybar,ymax
      real*8 tomass, totemp,anchor,xhi,xlo,slope,eeptag
      integer nstep,igrd
c To find the turnoff point
      ikip=1
      ymax=0.0
      do i=1,nstep
c      write(1,*)x(i),y(i)
      if(y(i).ge.ymax)then
       ikip=i
       ymax=y(i)
      endif
      enddo
cD
      if(ikip.ge.nstep)then
c this is when the isochrone does not have MSTO.
       xp=x(ikip)
       ybar=y(ikip)
       goto 88
      endif
c this is to find the MSTO
      call ToffM(x(ikip-1),y(ikip-1),xp,ybar,1,*99)
      if(xp.lt.x(ikip))call ToffM(x(ikip-2),y(ikip-2),xp,ybar,2,*99)
      if(xp.lt.x(ikip).and.xp.gt.x(ikip+1))then
        write(6,*)' possibly incorrect '
      endif
      go to 88
  99  continue
      xp=x(ikip)
  88  tomass=xp
c      totemp=ybar
      igrd=ikip
cD
c Got the turnoff point
c returns two value; igrd and tomass
c tomass is the turn off point
c igrd is the original grid number near by
      return
      end
      subroutine ToffM(x,yy,xp,ybar,iorder,*)
      parameter (n=4)
      real*8 x(n), y(n),yy(n),a,b,c,s,xp,xm,xbar,ybar
      integer k,l,np,m
      do k=1,n
       y(k)=yy(k)
c      write(6,*)x(k),y(k)
      enddo
      do 10 k=1,n-1
      do 20 l=1,(n-k)
        y(l)=(y(l+1)-y(l))/(x(l+k)-x(l))
   20 continue
   10 continue
      a=3.0d0*y(1)
      b=-2.0d0*y(1)*(x(4)+x(3)+x(2))+2.0d0*y(2)
      c=y(1)*(x(3)*x(4)+x(2)*x(4)+x(2)*x(3))-y(2)*(x(4)+x(3))+y(3)
      s=sqrt(b*b-4.0d0*a*c)
c      write(6,*)a,b,c,s
      xp=(-b+s)/(2.d0*a)
      xm=(-b-s)/(2.d0*a)
c      write(6,*)xp,xm
      if(xp.ge.x(iorder).and.xp.le.x(3))then
        xbar=xp
      else if(xm.ge.x(iorder).and.xm.le.x(3))then
       xbar=xm
      else
       if (iorder.ge.2)then
         return 1
       else
         write(6,*)' failed to find the turnoff mass'
         write(6,*)xp,xm
         stop
       endif
      endif
c      xp=(-2.0d0*c) /(b+s)
c      xm=(-2.0d0*c) /(b-s)
      ybar=y(1)
      do 100 m=2,n
       ybar=ybar*(xbar-x(m))+y(m)
  100 continue  
c      write(6,*)ybar
      xp=xbar
      return
      end
