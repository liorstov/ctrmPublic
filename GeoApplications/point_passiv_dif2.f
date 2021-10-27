c reflection receivers vertical streamers
       real cosphi, costhet,sinthet
       !h=20
       !v=1000.
       !xdif=14
       !ydif=55
       zdif=20
       !xdif=20
       xdif=100
       ydif=140
       !ydif=70
       v=2000
       vs=v*1./(3**0.5)
       print *,vs
       iflag=1
        !open(10,file='Ssave.prn')
        open(11,file='coord.txt')
        !!! if iflag =1  1 SYSPZ  )Y hit S P waves Z component  2-SYSPY 3-SYSPX
      
       xs1=0
       ys1=0
       zs1=0
c      print * ,ii, xs1,ys1,zs1
          
           ! read (10,*) ,ii, xs1,ys1,zs1
       !zs1=0
              t01=0.07
              t02=0.09
              t03=0.1
              t04=0.14
              t05=0.18
              t06=0.20
             
             ! ATT=  rand(20)*500
             ! print *,ATT,rand(20)
             ! read *
       do j=1,44
        read (11,*) ii, xr1,yr1,zr1
         
c      print * ,ii, xr1,yr1,zr1,i,j
c projection of distance to receiver on xy plane
c
           !do kkk=1,20
            
      !ds=sqrt((xs1-xdif)**2+( ys1-ydiff)**2+(zs1-h)**2)
        dr=sqrt((xr1-xdif)**2+( yr1-ydif)**2+(zr1-zdif)**2)
              AT= dr
              !find radiation pattern hit in direction of y S source  
        cosphi=(ydif-yr1)/
     *  ((xdif-xr1)**2+(ydif-yr1)**2)**0.5
        sinthet=(zdif-zr1)/((ydif-yr1)**2+(zdif-zr1)**2)**0.5
        if((ydif-yr1).le.0) sinthet=-sinthet
       costhet=(ydif-yr1)/((ydif-yr1)**2+(zdif-zr1)**2)**0.5 
       if(ydif.eq.yr1.and.xdif.eq.xr1) cosphi=0.
       sinphi=(xdif-xr1)/((ydif-yr1)**2+(xdif-xr1)**2)**0.5
       
       !!!! Source hit y S P wave Z component
       if (iflag.eq.1) then 
       uspp=sinthet*costhet*cosphi
	   !uspp = 1
       amplp=-uspp/dr
       !if S
       ampls=-amplp
       print *,'j',j
                      endif
       
       
      if (iflag.eq.3) then 
      !!! Source hit y P wave x component
        ! if((xdif-xr1).le.0) sinphi=-sinphi
         ! sinphi=-sinphi
          if((ydif-yr1).le.0) cosphi=-cosphi
        uspp=costhet*cosphi*sinphi
       amplp=-uspp/dr
       !!!! Source hit y S wave x component
       !if((xdif-xr1).le.0) sinphi=-sinphi
         if((ydif-yr1).le.0) sinthet=-sinthet
      uspp=sinthet*cosphi*sinphi
       !ampls=uspp/dr
       ampls=-amplp
       print *,'j',j
                      endif
                      if(iflag.eq.2) then
              !!!! Source hit y S wave y component
              if((ydif-yr1).le.0) sinthet=-sinthet
       uspp=sinthet*cosphi**2
       !uspp=sinthet**2*cosphi**2
       !print *,'j,sinthet, yr1,xr1',j,sinthet, cosphi,yr1,xr1
       ampls=-uspp/dr
      
      !!!! Source hit y P wave y component
       if((ydif-yr1).le.0) costhet=-costhet
       uspp=costhet*cosphi**2
       amplp=-uspp/dr
       print *,'j',j
                                      endif
      if(uspp.eq.0.or.dr.eq.0) then
      !amplp=0.0000001
       !ampls=0.000001
      endif
       tt=AT/v
       tts=AT/vs
       t1=tt+t01
       t2=tt+t02
       t2=tt+t02
       t3=tt+t03
       t4=tt+t04
       t5=tt+t05
       t6=tt+t06
       
          !  print *,'j at tt urpp', at,tt,urpp
       !enddo
       !t=AT/v
      !print *,t,iiii
c       print *,geom
c      read *
c      t=sqrt(geom)/v
      ! write(13,*)  ampl, tt,xr1,yr1
      write (13,*) amplp,ampls,tt,tts
        enddo
      
        
        
       
       
       stop
       end
       
