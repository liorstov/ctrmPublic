	dimension a(2500),t(2500,4 )
         character *80 filerad,filef
        integer*2 bb(1000)
         print *,'ENTER THE NAME FOR THE OUTPUT  FILE'
        read *,filerad
********************************
         print *,'ENTER THE NUMBER OF TRACES'
        read *,ntot
       ! ntot=184*92
         print *,'ENTER THE NUMBER OF SAMPLES'
        read *,nsamp
        dt=0.0005
        open(20,file=filerad,form='unformatted',
     *    access='direct'
     *    ,recl=nsamp*4)
********************************
        print *,'ENTER THE LEVEL OF NOISE 1 HIGH 10 LOW'
        read *,ilev
        irec=0
        do 77 kk=1,ntot
        do i=1,nsamp
        a(i)=0.
        enddo
        call randtr(irec,a,nsamp,dt,2,80.,1,1,ilev)
        irec=irec+1
        do 17 j=1,nsamp
17      bb(j)=a(j)*10000.
100       write(20,rec=irec)(a(j),j=1,nsamp)
        print *,'SHOT # ',ns,irec
77      continue
        print *,'total number of records - ',irec
        stop
        end
       subroutine randtr(irec,trace,nsamp,dt,iwave,hcut,jump,ist,ilev)
********************************************************
**
**    trace  -  input trace
**    nsamp  -  # of samples
**    dt     -  time interval (sec)
**    iwave  -  1 - gaussian
**              2 - ricker
**              3 - min. phase
**    hcut   -  highcut frequency (Hz)
**    jump   -  interval for random samples
**    ist    -  first sample to start noise
**
**    niose level:  0.- 1.
********************************************************
        save
        dimension trace(1),wavelt(500)
        data isss/0/
        if(isss.eq.0)then
        isss=1
	lx=256+irec
        endif
	do 3 j=ist,nsamp,jump
	rrr=random(lx)
3       trace(j)=trace(j)+rrr/ilev
        call wave(wavelt,nwavlt,iwave,hcut,dt)
        call convol(trace,wavelt,nsamp,nwavlt)
	return
	end
************************************************	
      SUBROUTINE CONVOL(A,B,NA,NB)
      dIMENSION A(1),B(1),D(7000)
      ishift=nb/2
      ishift=0
      MM=NA+NB-1
      DO 1 I=1,MM
1     D(I)=0.
      DO 2 I=1,NA
      DO 2 J=1,NB
 2    D(I+J-1)=D(I+J-1)+A(I)*B(J)
      DO 5 I=1,NA
 5    A(I)=D(I+ishift)
      RETURN
      eND
***********************************************
        SUBROUTINE WAVE(wavelt,nwavlt,iwave,hcut,dt)
***********************************************
**
**     iwave   :  1 - gaussian
**                2 - ricker
**                3 - min. phase
**
***********************************************
        DIMENSION WAVELT(1),wmin(16)
c        data wmin/0.01,0.7,1.,0.47,-0.45,-0.9,-0.62,-0.25,0.3,
c     .            0.4,0.37,0.24,0.18,0.12,0.08,0.01/
        data wmin/1.,1.368,0.6765,-0.6542,-1.221,-0.8854,-0.3306,
     .            0.4044,0.5619,0.5086,0.3384,0.2478,0.17,0.1096,
     .            0.0148,0.00019/
        PI=3.14159265
        TWPI=2.*PI
        if(iwave.eq.1)then
        agauss=hcut*2.
        tcut=3./agauss
        else
        agauss=hcut*0.5
        tcut=1.5/agauss
        endif
        NWAVLT=TCUT*2/DT
        IF(IWAVE.EQ.1) GO TO 50
        IF(IWAVE.EQ.3) GO TO 60
C       ****************
C       ...RICKER.......
C       ****************
        NW=NWAVLT/2
        DO 10 I=1,NWAVLT
        S=(-NW+I-1)*DT*AGAUSS
10      WAVELT(I)=EXP(-2.*S*S)*COS(TWPI*S)
        RETURN
C       ****************
C       ...GAUSSIAN.....
C       ****************
50      DO 11 I=1,NWAVLT
        X=(I-NWAVLT/2)*DT*AGAUSS
11      WAVELT(I)=EXP(-X*X*0.5)
        RETURN
C       ****************
C       ...MIN. PHASE...
C       ****************
60      nwavlt=16
        do 12 i=1,nwavlt
12      wavelt(i)=wmin(i)/1.368
        return
        END
*********************************************
	function random(lx)
	lx=mod(164907*lx,2147483647)
	random=float(lx)*4.6566128752458e-10
	return
	end
