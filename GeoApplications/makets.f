c              
	  real a(30000)
		character *80   filin
		   
        open(10,file='fort.15',form='formatted')
c*******************************
c         
c              
        dt=0.0005
c        
c*******************************
       print *,'number of samples'
       read *,nsamp
       print *,'number of shots'
       read *, nsh
       print *,'number of receivers'
       read *, nrec
       print *,'number of layers'
       read *, nlay
		print *,'file input'
		read * ,filin	
		open(10,file=filin,form='formatted')		
		print *,'file out'
		read * ,filin	
       open(11,file=filin,form='unformatted',
     .       access='direct',recl=nsamp*4)
     
ch      xs=20.
ch      ys=30.
ch      p=-5.
ch      h=10.
ch      v=200.
ch      do dik=-100,100,0.1
ch      dxh=sqrt((dik-xs)**2+(p-ys)**2+h**2)
c
c

ch      do i=1,96
ch      do j=1,48
ch      xr=(i-1.)*10.
ch      yr=(j-1.)*10.
ch      dlr=sqrt((xr-p)**2+h**2)
ch      ttt=(1./v)*sqrt((dlr+dls)**2+(yr-ys)**2)
c 
c
c      t1=10000
c      do k=1,4004
c      q=-1001.+k/2.
c       print *,q
c      s1=sqrt((xs-q)**2+(ys-p)**2+h**2)
c      s2=sqrt((xr-q)**2+(yr-p)**2+h**2)
c      c=(1./v)*(s1+s2)
c      if (c.lt.t1) then
c      s1b=s1
c      s2b=s2
c      qb=q
c                   endif
c      t1=c
c      enddo
ch      print*, ttt-t ,i,j
ch      enddo
ch      enddo
ch      stop
ch      end

c*******************************
        irec=0
        irecc=0
        do 77  ns=1,nsh
        do 77 i=1,nrec
            do 1 j=1,nsamp
1       a(j)=0.
            do 88 il=1,nlay
c       if(i.eq.1) irecc=irecc+1
         irecc=irecc+1
        read(10,*), xs,xr,xss,xrr,tas,taa
		print *, xs,xr,xss,xrr,tas,taa
c       do 100 i=1,4608
        
      if(il.eq.1) it1=taa/dt +0.5
      if(il.eq.2)  it2=taa/dt +0.5
      
      if(il.eq.3)  it3=taa/dt +0.5
c       if(il.eq.1) a(it1)=1.0/sqrt((xs-xss)**2+(xr-xrr)**2)**2*100000.
       if(il.eq.1) a(it1)=1.0
       if(il.eq.2) a(it2)=1.0
       if(il.eq.3) a(it3)=1.0
 88      continue
         iwave=4
         hcut=40.
        call synttr(a,nsamp,dt,iwave,hcut)
c       call randtr(a,nsamp,dt,2,25.,8,ist)
         
c       if(i.eq.1) write(12,rec=irecc) (a(j),j=1,nsamp)
         write(11,rec=irecc)(a(j),j=1,nsamp)
        print *,'SHOT # ',ns,irecc
77      continue
        print *,'total number of records - ',irec
        close(11)
        stop
        end
                 
          SUBROUTINE WAVE(wavelt,nwavlt,iwave,hcut,dt)
c***********************************************
c**     iwave   :  1 - gaussian
c**                2 - ricker
c**                3 - min. phase
c**                4- vladi signal
c***********************************************
         DIMENSION WAVELT(1),wmin(16)
        data wmin/1.,1.368,0.6765,-0.6542,-1.221,-0.8854,-0.3306, 
     *            0.4044,0.5619,0.5086,0.3384,0.2478,0.17,0.1096, 
     *            0.0148,0.00019/
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
          IF(IWAVE.EQ.2) GO TO 40
          IF(IWAVE.EQ.3) GO TO 60
          if(IWAVE.EQ.4) then
c***********
c...vladi...
c***********    NF=100
          NT=50
          nf=100
          nwavlt=nf
          do 111 I=1,NF
          dd=i
          wavelt(i)=((1.-dd/(NF-1))**2)*sin(TWPI*dd/NT)
111           continue
          endif
          return
 
c************
c...RICKER...
c************
40       NW=NWAVLT/2
       DO 10 I=1,NWAVLT
       S=(-NW+I-1)*DT*AGAUSS
10         WAVELT(I)=EXP(-2.*S*S)*COS(TWPI*S)
       RETURN
c**************
c...GAUSSIAN...
c**************
50      DO 11 I=1,NWAVLT
        X=(I-NWAVLT/2)*DT*AGAUSS
11      WAVELT(I)=EXP(-X*X*0.5)
        RETURN
c****************
c...MIN. PHASE...
c****************
60       nwavlt=25
         do 12 i=1,nwavlt
12           wavelt(i)=wmin(i)/1.368
         return
        END
       subroutine synttr(trace,nsamp,dt,iwave,hcut)
c*******************************************************
c*    convolution of a trace with synthetic wavelet.
c*
c*    trace  -  input trace
c*    nsamp  -  # of samples
c*    dt     -  time interval (sec)
c*    iwave  -  1 - gaussian
c*              2 - ricker
c*              3 - min. phase
c*    hcut   -  highcut frequency (Hz)
c*
c*******************************************************
        dimension trace(3000),wavelt(500)
        call wave(wavelt,nwavlt,iwave,hcut,dt)
        call convol(trace,wavelt,nsamp,nwavlt)
	    return
	    end
c***********************************************	
      SUBROUTINE CONVOL(A,B,NA,NB)
C     CONVOLUTION OF A AND B. RESULT IN A. NA,NB LENGTH OF A AND B.
      dIMENSION A(1),B(1),D(70000)
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
       