c reflection receivers vertical streamers
       h=12
       v=2000.
       xdif=50
       ydif=25
        open(10,file='Ssave.prn')
        open(11,file='Rsave.prn')
      
       
c      print * ,ii, xs1,ys1,zs1
        do iiii=1,210
            read (10,*) ,ii, xs1,ys1,zs1
       zs1=0
             
             ! ATT=  rand(20)*500
             ! print *,ATT,rand(20)
             ! read *
       do j=1,432
        read (11,*) ii, xr1,yr1,zr1
        zr1=0
c      print * ,ii, xr1,yr1,zr1,i,j
c projection of distance to receiver on xy plane
c
           !do kkk=1,20
             ydiff=(kkk-1)+ydif
      ds=sqrt((xs1-xdif)**2+( ys1-ydiff)**2+(zs1-h)**2)
        dr=sqrt((xr1-xdif)**2+( yr1-ydiff)**2+(zr1-h)**2)
              AT= dr+ds
              
          
       t=AT/v
          print *,'kkk  t j iii', kkk,t,j,iiii
       !enddo
       !t=AT/v
      !print *,t,iiii
c       print *,geom
c      read *
c      t=sqrt(geom)/v
       write(12,*),xs1,ys1,xr1,yr1,zr1,t, t
        enddo
        
       rewind(11)
       
       enddo
       stop
       end
       
