c reflection receivers vertical streamers
       !h=20
       !v=1000.
       !xdif=14
       !ydif=55
       h=20
       !xdif=20
       xdif=44
       ydif=-18
       !ydif=70
       v=2000
        !open(10,file='Ssave.prn')
        open(11,file='coord.txt')
      
        
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
       do
			read (11,*,end = 1) dii, xr1,yr1,zr1

			c zr1=0
			print * ,ii, xr1,yr1,zr1,i
			
			!do kkk=1,20

			!ds=sqrt((xs1-xdif)**2+( ys1-ydiff)**2+(zs1-h)**2)
			dr=sqrt((xr1-xdif)**2+( yr1-ydif)**2+(zr1-h)**2)
			AT= dr


			!tt=AT/v-h/v
			tt=AT/v 
			!t1=tt+t01
			!t2=tt+t02
			!t3=tt+t03
			! t4=tt+t04
			!t5=tt+t05
			!t6=tt+t06

			print *,'kkk  t  iii', kkk,t,iiii
			!enddo
			!t=AT/v
			!print *,t,iiii
c       print *,geom
c      read *
c      t=sqrt(geom)/v
			write(15,*) xs1,ys1,xr1,yr1, tt, tt
		enddo
1		print *,'eof'
       
       stop
       end
       
