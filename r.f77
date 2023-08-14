cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                                                                                    c
c          the model is described as the following																				     c
c          1. the model is two dimension																							 c
c  	       2. the length in x direction is 40 meters,the length in y direction is the same								             c
c     	   3. the square is 10 m * 10 m, so the number o grid is 5 *5															     c
c          4. the interval between two shots is 3 meters,the interval between two geophones is 1 meters						      	 c
c          5. the first shot is at 2 meters in y direction,and last one is at 38 ,so the numberof shots is 13						 c
c          6. tge first geophone is at 1 meter in y direction ,and the last one is at 39 meter ,so the number of geophone is 39  	 c
c																																     c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      PROGRAM first_program
cccccccc     1.  vaviables description     ccccccccccccccccccccccccccccccccccccccccccccc

      parameter(num_s=5,num_g=5,igrid_r=21,igrid_v=21)
      parameter(d_x=20,d_y=40)
      real,parameter::pi=3.1415926
      real t0(num_s,num_g),t1(num_s,num_g),t2(num_s,num_g)
      real dt(num_s,num_g)
      real v(igrid_r,igrid_v),et(2*num_s*num_g,igrid_r*igrid_v)
      real aa(igrid_r*igrid_v,igrid_r*igrid_v),b(igrid_r*igrid_v)
      real aet(igrid_r*igrid_v,2*num_s*num_g),xc(igrid_r,igrid_v)
      real f(2*num_s*num_g),af(igrid_r*igrid_v)
      real vv(igrid_r,igrid_v),vt(igrid_r,igrid_v)
      real xcoe1(10000),ycoe1(10000),xcoe2(10000),ycoe2(10000)
      real tcoe1(10000),tcoe2(10000),dis(10000),ang(10000)
      real pcoe1(10000),pcoe2(10000)
      real p0(num_s,num_g),p1(num_s,num_g),p2(num_s,num_g)
      real ang0(num_s,num_g),ang1(num_s,num_g), ang2(num_s,num_g)
      real zn(1000),zni(1000)
      real tio(igrid_r*igrid_v),tia(igrid_r*igrid_v)
      real w1(10000),eet(igrid_r*igrid_v,num_s*num_g*2)
      real ww(num_s*num_g*2,num_s*num_g*2)
      real eet1(igrid_r*igrid_v,num_s*num_g*2)
      real vvxx(igrid_r,igrid_v),vvzz(igrid_r,igrid_v)
      real vxx(igrid_r,igrid_v),  vzz(igrid_r,igrid_v)

      grid_interval_x=d_x/((igrid_v-1)*1.)
      grid_interval_y=d_y/((igrid_r-1)*1.)
      interval_s=(d_y-30)/((num_s-1)*1.)
      interval_g=(d_y-20)/((num_g-1)*1.)
      y_first_s=20
      y_first_g=10
      source_x=0
      geophone_x=d_x
      velocity_true=2000
      velocity_initial=0
      x_cord_begin=0
      x_cord_end=d_x
      y_cord_begin=0
      y_cord_end=d_y
      idt_scan=20
      a_c=1/1.
      itt=5
      ijk=1
      cof=1.
      dv=0.01
      dvz=dv
      dvx=dv*0
      ddv=200
      open(12,file='result311.dat')
      open(123,file='see-p-3.dat')
      open(15,file='result-2.dat')
      open(16,file='true-2.dat')
      open(7,file='c_saa-2-1-for-10.dat')
      open(18,file='curve-test1.dat')
      open(19,file='curve-test2.dat')
      open(31,file='c1.dat')
      open(32,file='c2.dat')
      open(33,file='cc1.dat')
      open(50,file='dis.dat')
      idd=0
      ttmax=0
      ppmax=0
      ttmin=100000000.
      ppmin=100000000.
      st=0
      sp=0
      ix=2
      rant=0.
      ranp=0.
      beta=0
cccccccc     2. define true velcoity model  cccccccccccccccccccccccccccccccccccccccccccc
      do i=1,igrid_r
         do j=1,igrid_v
            c_x=(j-1)*grid_interval_x
            c_y=(i-1)*grid_interval_y
            v(i,j)=velocity_true
            if(c_y.ge.120..and.c_y.le.280) v(i,j)=v(i,j)+100

            write(123,*) 'true',i,j,v(i,j)
         end do
      end do

      do i=1,igrid_r
         do j=1,igrid_v
         write(16,*) (j-1)*grid_interval_x,(i-1)*grid_interval_y,v(i,j)
         write(16,*) 'qwqwqw',i,j,v(i,j)
         end do
      end do
cccccccc     2. get the true traveltime in geophone cccccccccccccccccccccccccccccccccccc
      write(*,*) 'begin to get the true traveltime'
      iio=0
      do i=num_s,1,-1
         do j=1,num_g
            source_y=y_first_s+(i-1)*interval_s
            geophone_y=y_first_g+(j-1)*interval_g
            s=1.0E+10
            if((geophone_y-source_y).ne.0) then
               a1=atan((x_cord_end-x_cord_begin)/(geophone_y-source_y))
     *                /pi*180.
               ia1=int(a1)
               if(ia1.lt.0) ia1=ia1+180
            else
               ia1=90
            end if
            jk1=(ia1-idt_scan)*10
            jk2=(ia1+idt_scan)*10
            if(jk1.lt.100) jk1=100
            if(ik2.gt.1700) jk2=1700
            i_a=0
            do k=jk1,jk2,itt
               i_a=i_a+1
               angle=pi/180.*(k/10.)
c	            write(*,*) 'pass1'
            call ray11(x_cord_begin,y_cord_begin,x_cord_end,y_cord_end,
     *		            igrid_r,igrid_v,
     *			         grid_interval_x,grid_interval_y,
     *					   source_x,source_y,
     *	               angle,
     *	               v,
     *					   x,y,t,p,pz,xx,yy,tt,pp,ppz,idd)

               if(x-geophone_x.lt.0) then
                  dis(i_a)=1000000000000000.
                  cycle
               end if
               s1=sqrt((x-geophone_x)**2+(y-geophone_y)**2)
     *		      /velocity_true+t
               dis(i_a)=s1
               xcoe1(i_a)=xx
               ycoe1(i_a)=yy
               tcoe1(i_a)=tt
               if(ppz.gt.pp) ppz=pp
               if(ppz.lt.pp.and.abs(ppz).gt.pp) ppz=-pp
               pcoe1(i_a)=acos(ppz/pp)/pi*180
               xcoe2(i_a)=x
               ycoe2(i_a)=y
               tcoe2(i_a)=t
               if(pz.gt.pp) pz=pp
               if(pz.lt.pp.and.abs(pz).gt.pp) pz=-pp
               pcoe2(i_a)=acos(pz/pp)/pi*180
               ang(i_a)=k/10.
            end do
c	         write(*,*) 'pass3'
            do m=1,i_a
               if(dis(m).lt.s) then
                  k_min1=m
                  s=dis(m)
               end if
            end do
            ss=1.0E+10
            do m=1,i_a
               if(m.ne.k_min1) then
                  if(dis(m).lt.ss) then
                     k_min2=m
                     ss=dis(m)
                  end if
               end if
            end do
c	         write(*,*) 'pass4'
            sas1=(tcoe2(k_min1)-tcoe2(k_min1-1))**2
     *	       +(pcoe2(k_min1)-pcoe2(k_min1-1))**2
            sas2=(tcoe2(k_min1)-tcoe2(k_min1+1))**2
     *           +(pcoe2(k_min1)-pcoe2(k_min1+1))**2
            if(sas1.gt.sas2) then
               k_min2=k_min1+1
            else
               k_min2=k_min1-1
            end if
            ttt2_1=tcoe1(k_min2)
            xxx2_1=xcoe1(k_min2)
            yyy2_1=ycoe1(k_min2)
            ppp2_1=pcoe1(k_min2)
            ttt2_2=tcoe2(k_min2)
            xxx2_2=xcoe2(k_min2)
            yyy2_2=ycoe2(k_min2)
            ppp2_2=pcoe2(k_min2)
            ttt1_1=tcoe1(k_min1)
            xxx1_1=xcoe1(k_min1)
            yyy1_1=ycoe1(k_min1)
            ppp1_1=pcoe1(k_min1)
            ttt1_2=tcoe2(k_min1)
            xxx1_2=xcoe2(k_min1)
            yyy1_2=ycoe2(k_min1)
            ppp1_2=pcoe2(k_min1)
c	         write(*,*) 'pass5',xxx1_2,xxx1_1

         tt1=ttt1_1+(ttt1_2-ttt1_1)/(xxx1_2-xxx1_1)*(geophone_x-xxx1_1)
         pp1=ppp1_1+(ppp1_2-ppp1_1)/(xxx1_2-xxx1_1)*(geophone_x-xxx1_1)
         yy1=yyy1_1+(yyy1_2-yyy1_1)/(xxx1_2-xxx1_1)*(geophone_x-xxx1_1)
c	         write(*,*) 'pass5-1',xxx1_2,xxx1_1
         tt2=ttt2_1+(ttt2_2-ttt2_1)/(xxx2_2-xxx2_1)*(geophone_x-xxx2_1)
         pp2=ppp2_1+(ppp2_2-ppp2_1)/(xxx2_2-xxx2_1)*(geophone_x-xxx2_1)
         yy2=yyy2_1+(yyy2_2-yyy2_1)/(xxx2_2-xxx2_1)*(geophone_x-xxx2_1)
c 	         write(*,*) 'pass5-2'
            tt=tt1+(tt2-tt1)/(yy2-yy1)*(geophone_y-yy1)
            pp=pp1+(pp2-pp1)/(yy2-yy1)*(geophone_y-yy1)
c  	      write(*,*) 'pass5-3'
            if(abs(tt).gt.ttmax) ttmax=abs(tt)
            if(abs(pp).gt.ppmax) ppmax=abs(pp)
            if(abs(tt).lt.ttmin) ttmin=abs(tt)
            if(abs(pp).lt.ppmin) ppmin=abs(pp)
c	         write(*,*) 'passs6'
            st=st+abs(tt)
            sp=sp+abs(pp)
            t0(i,j)=tt
            p0(i,j)=pp/cof
            ang0(i,j)=ang(k_min1)+(ang(k_min2)-ang(k_min1))
     *		/(yy2-yy1)*(geophone_y-yy1)
            ang0(i,j)=ang0(i,j)/cof
            iio=iio+1
            if(j.eq.8) then
               write(33,*) i,p0(i,j),ang0(i,j)
            end if
         end do
      end do

      do i=1,num_s
         do j=1,num_g
            call rann(ix,yfl)
            ran_t=(yfl-0.5)*2*rant
            t0(i,j)=t0(i,j)+ran_t
            call rann(ix,yfl)

            ran_p=(yfl-0.5)*2*ranp
            p0(i,j)=p0(i,j)+ran_p
	      end do
      end do


cccccccccc    initial model cccccccccccccccccccccccccccccccccccccccccc

      write(*,*) 'begin to calculate the traveltime for initial model'
      ii_k=0
c     open(186,file='qdf1.dat')
      do ii=1,igrid_r
         do jj=1,igrid_v
            call rann(ix,yfl)
            ran_v=(yfl-0.5)*0
            v(ii,jj)=velocity_initial+ran_v
            c_x=(jj-1)*grid_interval_x
            c_y=(ii-1)*grid_interval_y
            vvxx(ii,jj)=velocity_initial
            vvzz(ii,jj)=velocity_initial
            write(123,*) 'imitial',ii,jj,v(ii,jj)
         end do
      end do

      ks=0
      sss_min=0
      ign=0
198   ks=ks+1
      if(ign.eq.-1) a_c=a_c*0.1
      write(7,*) 'aaaa',ks,it_min,sss_min,ign,a_c
      do i=1,igrid_r
	      do j=1,igrid_v
         write(7,*) (j-1)*grid_interval_x,(i-1)*grid_interval_y,v(i,j)
	      end do
      end do

      write(*,*) 'the number of iteration is ',ks
      ii_p=0
      do i=1,num_s
         do j=1,num_g
            ii_p=ii_p+1
            source_y=y_first_s+(i-1)*interval_s
            geophone_y=y_first_g+(j-1)*interval_g
            s=1.0E+10
            if((geophone_y-source_y).ne.0) then
               a1=atan((x_cord_end-x_cord_begin)/(geophone_y-source_y))
     *            /pi*180.
               ia1=int(a1)
               if(ia1.lt.0) ia1=ia1+180
            else
               ia1=90
            end if
            jk1=(ia1-idt_scan)*10
	         jk2=(ia1+idt_scan)*10
	         if(jk1.lt.100) jk1=100
	         if(ik2.gt.1700) jk2=1700
            i_a=0
               do k=jk1,jk2,itt
                  i_a=i_a+1
                  angle=pi/180.*(k/10.)
                  if(i.eq.3.and.j.eq.1) then
                     idd=-90
                  end if
                  call ray1(x_cord_begin,y_cord_begin,
     *                     x_cord_end,y_cord_end,
     *			            igrid_r,igrid_v,
     *			            grid_interval_x,grid_interval_y,
     *					      source_x,source_y,
     *	                  angle,
     *	                  vvxx,vvzz,
     *					      x,y,t,p,pz,xx,yy,tt,pp,ppz,idd)
                  if(x-geophone_x.lt.0) then
                     dis(i_a)=1000000000000000.
                     cycle
                  end if
            s1=sqrt((x-geophone_x)**2+(y-geophone_y)**2)/velocity_true
     *		         +t
                  dis(i_a)=s1
                  xcoe1(i_a)=xx
                  ycoe1(i_a)=yy
                  tcoe1(i_a)=tt
                  if(ppz.gt.pp) ppz=pp
                  if(ppz.lt.pp.and.abs(ppz).gt.pp) ppz=-pp
                  pcoe1(i_a)=acos(ppz/pp)/pi*180
                  xcoe2(i_a)=x
                  ycoe2(i_a)=y
                  tcoe2(i_a)=t
                  if(pz.gt.pp) pz=pp
                  if(pz.lt.pp.and.abs(pz).gt.pp) pz=-pp
                  pcoe2(i_a)=acos(pz/pp)/pi*180
                  ang(i_a)=k/10.
            end do
            do m=1,i_a
               if(dis(m).lt.s) then
                  k_min1=m
                  s=dis(m)
               end if
            end do
            ss=1.0E+10
            do m=1,i_a
               if(m.ne.k_min1) then
                  if(dis(m).lt.ss) then
                     k_min2=m
                     ss=dis(m)
                  end if
               end if
            end do
            sas1=(pcoe2(k_min1)-pcoe2(k_min1-1))**2
            sas2=(pcoe2(k_min1)-pcoe2(k_min1+1))**2
            if(sas1.gt.sas2) then
               k_min2=k_min1+1
            else
               k_min2=k_min1-1
            end if

            ttt2_1=tcoe1(k_min2)
            xxx2_1=xcoe1(k_min2)
            yyy2_1=ycoe1(k_min2)
            ppp2_1=pcoe1(k_min2)
            ttt2_2=tcoe2(k_min2)
            xxx2_2=xcoe2(k_min2)
            yyy2_2=ycoe2(k_min2)
            ppp2_2=pcoe2(k_min2)
            ttt1_1=tcoe1(k_min1)
            xxx1_1=xcoe1(k_min1)
            yyy1_1=ycoe1(k_min1)
            ppp1_1=pcoe1(k_min1)
            ttt1_2=tcoe2(k_min1)
            xxx1_2=xcoe2(k_min1)
            yyy1_2=ycoe2(k_min1)
            ppp1_2=pcoe2(k_min1)

      tt1=ttt1_1+(ttt1_2-ttt1_1)/(xxx1_2-xxx1_1)*(geophone_x-xxx1_1)
      pp1=ppp1_1+(ppp1_2-ppp1_1)/(xxx1_2-xxx1_1)*(geophone_x-xxx1_1)
      yy1=yyy1_1+(yyy1_2-yyy1_1)/(xxx1_2-xxx1_1)*(geophone_x-xxx1_1)
      tt2=ttt2_1+(ttt2_2-ttt2_1)/(xxx2_2-xxx2_1)*(geophone_x-xxx2_1)
      pp2=ppp2_1+(ppp2_2-ppp2_1)/(xxx2_2-xxx2_1)*(geophone_x-xxx2_1)
      yy2=yyy2_1+(yyy2_2-yyy2_1)/(xxx2_2-xxx2_1)*(geophone_x-xxx2_1)
      tt=tt1+(tt2-tt1)/(yy2-yy1)*(geophone_y-yy1)
      pp=pp1+(pp2-pp1)/(yy2-yy1)*(geophone_y-yy1)
            t1(i,j)=tt
            p1(i,j)=pp/cof
            ang1(i,j)=ang(k_min1)+(ang(k_min2)-ang(k_min1))
     *		      /(yy2-yy1)*(geophone_y-yy1)
            ang1(i,j)=ang1(i,j)/cof
            bbb=10000
            if(abs(p1(i,j)).eq.90) then
               ukk=1
            else
               ukk=1./abs(p1(i,j)-90)
            end if
            if(abs(ang1(i,j)).eq.90) then
               wkk=1
            else
               wkk=1./abs(ang1(i,j)-90)
            end if
            w1(ii_p)=1
            w1(num_s*num_g+ii_p)=1
            if(ks.eq.1.and.ijk.eq.2) then
            sss_min=sss_min+(t1(i,j)-t0(i,j))**2+(p1(i,j)-p0(i,j))**2
            end if
            if(ks.eq.1.and.ijk.eq.1) then
               sss_min=sss_min+(p1(i,j)-p0(i,j))**2
            end if
         end do
      end do

cccccccccccccccc  calculating the derivation         ccccccccccccccccccccccccccccccccccccccccccccccc
      write(*,*) 'begin to calculate the deviation'
*$acc data copy(vxx,vzz,vvxx,vvzz,vvxx1,vvzz1,x_ray,z_ray,t,pp,ppx,ppz)
      ii_k=0
      ia=0
      do iijj=2,2,4
         do ii=1,igrid_r
            do jj=1,igrid_v
               ia=ia+1
               qsc=1+(jj-1)*beta
               ii_k=ii_k+1
               xvv=0
               ixvv=0
               do i=1,igrid_r
                  do j=1,igrid_v
                     vxx(i,j)=vvxx(i,j)
                     vzz(i,j)=vvzz(i,j)
                     if(iijj.eq.1) then
                        if(i.eq.ii.and.j.eq.jj) then
                           vxx(i,j)=vxx(i,j)+dvx
                        end if
                     end if
                     if(iijj.eq.2) then
                        if(i.eq.ii.and.j.eq.jj) then
                           vzz(i,j)=vzz(i,j)+dvz
                        end if
                     end if
                  end do
               end do
               i_k=0
c	            write(*,*) 'pass1'
               do i=1,num_s
                  do j=1,num_g
                     i_k=i_k+1
                     source_y=y_first_s+(i-1)*interval_s
                     geophone_y=y_first_g+(j-1)*interval_g
                     s=1.0E+10
                     if((geophone_y-source_y).ne.0) then
                        a1=atan((x_cord_end-x_cord_begin)
     *		               /(geophone_y-source_y))
                        a1=a1/pi*180.
                        ia1=int(a1)
                        if(ia1.lt.0) ia1=ia1+180
                     else
                        ia1=90
                     end if
                     jk1=(ia1-idt_scan)*10
                     jk2=(ia1+idt_scan)*10
                     if(jk1.lt.100) jk1=100
                     if(ik2.gt.1700) jk2=1700
                     i_a=0
*$acc region
*$acc loop private(i_a)
                     do k=jk1,jk2,itt
                        i_a=i_a+1
                        angle=pi/180.*(k/10.)
c                       write(*,*) 'pass2',iijj,ii,jj,i,j
                        call ray1(x_cord_begin,y_cord_begin,
     *				         x_cord_end,y_cord_end,
     *			            igrid_r,igrid_v,
     *			            grid_interval_x,grid_interval_y,
     *					      source_x,source_y,
     *	                  angle,
     *	                  vxx,vzz,
     *					      x,y,t,p,pz,xx,yy,tt,pp,ppz,idd)
                        if(x-geophone_x.lt.0) then
                           dis(i_a)=1000000000000000.
                           cycle
                        end if
                        s1=sqrt((x-geophone_x)**2+(y-geophone_y)**2)
     *				            /velocity_true+t
                        dis(i_a)=s1
                        xcoe1(i_a)=xx
                        ycoe1(i_a)=yy
                        tcoe1(i_a)=tt
                        if(ppz.gt.pp) ppz=pp
                        if(ppz.lt.pp.and.abs(ppz).gt.pp) ppz=-pp
                        pcoe1(i_a)=acos(ppz/pp)/pi*180
                        xcoe2(i_a)=x
                        ycoe2(i_a)=y
                        tcoe2(i_a)=t
                        if(pz.gt.pp) pz=pp
                        if(pz.lt.pp.and.abs(pz).gt.pp) pz=-pp
                        pcoe2(i_a)=acos(pz/pp)/pi*180
                        ang(i_a)=k/10.
c                    write(*,*) 'pass4'
                  end do
*$acc end region
                  do m=1,i_a
                     if(dis(m).lt.s) then
                        k_min1=m
                        s=dis(m)
                     end if
                  end do
                  ss=1.0E+10
                  do m=1,i_a
                     if(m.ne.k_min1) then
                        if(dis(m).lt.ss) then
                           k_min2=m
                           ss=dis(m)
                        end if
                     end if
                  end do
c	               write(*,*) 'pass5'
                  sas1=(tcoe2(k_min1)-tcoe2(k_min1-1))**2
     *	             +(pcoe2(k_min1)-pcoe2(k_min1-1))**2
                  sas2=(tcoe2(k_min1)-tcoe2(k_min1+1))**2
     *                +(pcoe2(k_min1)-pcoe2(k_min1+1))**2
                  if(sas1.gt.sas2) then
                     k_min2=k_min1+1
                  else
                     k_min2=k_min1-1
                  end if
c	               write(*,*) 'pass6'
                  ttt2_1=tcoe1(k_min2)
                  xxx2_1=xcoe1(k_min2)
                  yyy2_1=ycoe1(k_min2)
                  ppp2_1=pcoe1(k_min2)
                  ttt2_2=tcoe2(k_min2)
                  xxx2_2=xcoe2(k_min2)
                  yyy2_2=ycoe2(k_min2)
                  ppp2_2=pcoe2(k_min2)
                  ttt1_1=tcoe1(k_min1)
                  xxx1_1=xcoe1(k_min1)
                  yyy1_1=ycoe1(k_min1)
                  ppp1_1=pcoe1(k_min1)
                  ttt1_2=tcoe2(k_min1)
                  xxx1_2=xcoe2(k_min1)
                  yyy1_2=ycoe2(k_min1)
                  ppp1_2=pcoe2(k_min1)
c	               write(*,*) 'pass7'
                  tt1=ttt1_1+(ttt1_2-ttt1_1)/(xxx1_2-xxx1_1)
     *				    *(geophone_x-xxx1_1)
                  pp1=ppp1_1+(ppp1_2-ppp1_1)/(xxx1_2-xxx1_1)
     *				    *(geophone_x-xxx1_1)
                  yy1=yyy1_1+(yyy1_2-yyy1_1)/(xxx1_2-xxx1_1)
     *				    *(geophone_x-xxx1_1)
                  tt2=ttt2_1+(ttt2_2-ttt2_1)/(xxx2_2-xxx2_1)
     *				    *(geophone_x-xxx2_1)
                  pp2=ppp2_1+(ppp2_2-ppp2_1)/(xxx2_2-xxx2_1)
     *				    *(geophone_x-xxx2_1)
                  yy2=yyy2_1+(yyy2_2-yyy2_1)/(xxx2_2-xxx2_1)
     *				    *(geophone_x-xxx2_1)
                  tt=tt1+(tt2-tt1)/(yy2-yy1)*(geophone_y-yy1)
                  pp=pp1+(pp2-pp1)/(yy2-yy1)*(geophone_y-yy1)
                  t2(i,j)=tt
                  p2(i,j)=pp/cof
                  ang2(i,j)=ang(k_min1)+(ang(k_min2)-ang(k_min1))
     *		            /(yy2-yy1)*(geophone_y-yy1)
                  ang2(i,j)=ang2(i,j)/cof

                  f(i_k)=p1(i,j)-p0(i,j)
                  f(num_s*num_g+i_k)=ang1(i,j)-ang0(i,j)
                  et(i_k,ii_k)=(p2(i,j)-p1(i,j))/dv*qsc
               et(num_s*num_g+i_k,ii_k)=(ang2(i,j)-ang1(i,j))/dv*qsc

                     if(jj.eq.1.or.jj.eq.11) then
 				   et(i_k,ii_k)=(p2(i,j)-p1(i,j))/dv*qsc
               et(num_s*num_g+i_k,ii_k)=(ang2(i,j)-ang1(i,j))
     *			   /dv*qsc
                     end if
                     if(jj.eq.2.or.jj.eq.10) then
 				   et(i_k,ii_k)=(p2(i,j)-p1(i,j))/dv*qsc
               et(num_s*num_g+i_k,ii_k)=(ang2(i,j)-ang1(i,j))
     *			   /dv*qsc
                     end if
                     if(p2(i,j)-p1(i,j).ne.0) then
                        ixvv=ixvv+1
                        xvv=xvv+(p2(i,j)-p1(i,j))**2
                     end if
                     if(ang2(i,j)-ang1(i,j).ne.0) then
                        ixvv=ixvv+1
                        xvv=xvv+(ang2(i,j)-ang1(i,j))**2
                     end if
c	               write(*,*) 'pass8'
                  end do
               end do
c	            write(*,*) 'pass9'
               tio(ia)=sqrt(xvv)
               tia(ia)=ixvv
               if(jj.eq.1.or.jj.ge.9) then
                  tio(ia)=50
               else
                  tio(ia)=1
               end if
               write(31,*) ia,xvv
               write(32,*) ia,ixvv
               if(ixvv.ne.0) then
                  write(50,*) ii,jj,ixvv
                  write(18,*) ia,xvv*1.5*1000000
               end if
               zni(ii_k)=ixvv
               zn(ii_k)=sqrt(xvv)
c	            write(*,*) 'pass10'
            end do
         end do
      end do

*$acc end data
c	   write(*,*) 'pass11'

      xmax=0
      do i=1,igrid_r*igrid_v
         if(zni(i).gt.xmax) then
            xmax=zni(i)
         end if
      end do
c	   write(*,*) 'pass12'
      xmax=0
      do i=1,igrid_r*igrid_v
         if(tia(i).gt.xmax) then
            xmax=tia(i)
         end if
      end do
c	   write(*,*) 'pass13'

ccccccccccccccccccccccc    inversion   cccccccccccccccccccccccccccccccccccccccccccccccccccc
      write(*,*) 'begin to inverse '

      do i=1,num_s*num_g*2
         do j=1,igrid_r*igrid_v
            eet(j,i)=et(i,j)
         end do
      end do

      do i=1,num_s*num_g*2
         do j=1,num_s*num_g*2
            ww(i,j)=0
            if(i.eq.j) ww(i,j)=w1(i)
         end do
      end do

      do i=1,igrid_r*igrid_v
         do j=1,num_s*num_g*2
            eet1(i,j)=0
            do ii=1,num_s*num_g*2
               eet1(i,j)=eet1(i,j)+eet(i,ii)*ww(ii,j)
            end do
         end do
      end do

      tm_max=0
      do i=1,igrid_r*igrid_v
         do j=1,igrid_r*igrid_v
            aa(i,j)=0
            do ii=1,2*num_s*num_g
               aa(i,j)=aa(i,j)+eet1(i,ii)*et(ii,j)
            end do
            if(abs(aa(i,j)).gt.tm_max) tm_max=abs(aa(i,j))
         end do
      end do
      write(*,*) 'tm_max',tm_max
      ia=0
      do i=1,igrid_r*igrid_v
         ia=ia+1
         do j=1,igrid_r*igrid_v
            aaa=10
            bbb=aaa*tm_max
            if(i.eq.j) aa(i,j)=aa(i,j)+bbb
         end do
      end do


      do 13 i=1,igrid_r*igrid_v
         b(i)=1
13    continue
c	   write(*,*) 'pass1'

      call gaussj(aa,igrid_r*igrid_v,igrid_r*igrid_v,b,1,1)
c	   write(*,*) 'pass'
      do i=1,igrid_r*igrid_v
         do j=1,2*num_s*num_g
            aet(i,j)=0
            do k=1,igrid_r*igrid_v
               aet(i,j)=aet(i,j)+aa(i,k)*et(j,k)
            end do
         end do
      end do

      do i=1,igrid_r*igrid_v
         do j=1,num_s*num_g*2
            eet1(i,j)=0
            do ii=1,num_s*num_g*2
               eet1(i,j)=eet1(i,j)+aet(i,ii)*ww(ii,j)
            end do
         end do
      end do

      do i=1,igrid_r*igrid_v
         af(i)=0
         do j=1,2*num_s*num_g
            af(i)=af(i)+eet1(i,j)*f(j)
         end do
      end do

cccccccccccc one dimension searching cccccccccccccccccccccccccccccccc
      write(*,*) 'begin to one dimension searching'

      ign=-1
      it_min=0.
      sss_min=10000000000.
      do it=1,100,2
         write(16,*) 'it=',it
         sss=0
         i_k=0
         do iijj=2,2,4
            do i=1,igrid_r
               do j=1,igrid_v
                  i_k=i_k+1
                  xvf=((it-1)*1./10)*af(i_k)
                  if(iijj.eq.1) then
                     vxx(i,j)=vvxx(i,j)-xvf
                  end if
                  if(iijj.eq.2) then
                     vzz(i,j)=vvzz(i,j)-xvf
                  end if
                  if(abs(xvf).gt.0.1) goto 197
               end do
            end do
         end do
         do i=1,num_s
            do j=1,num_g
               source_y=y_first_s+(i-1)*interval_s
               geophone_y=y_first_g+(j-1)*interval_g
               s=1.0E+10
               if((geophone_y-source_y).ne.0) then
         a1=atan((x_cord_end-x_cord_begin)/(geophone_y-source_y))
     *            /pi*180.
                  ia1=int(a1)
                  if(ia1.lt.0) ia1=ia1+180
               else
                  ia1=90
               end if
               jk1=(ia1-idt_scan)*10
               jk2=(ia1+idt_scan)*10
               if(jk1.lt.100) jk1=100
               if(ik2.gt.1700) jk2=1700
               i_a=0
               do k=jk1,jk2,itt
                  i_a=i_a+1
                  angle=pi/180.*(k/10.)
                  call ray1(x_cord_begin,y_cord_begin,
     *				         x_cord_end,y_cord_end,
     *			            igrid_r,igrid_v,
     *			            grid_interval_x,grid_interval_y,
     *					      source_x,source_y,
     *	                  angle,
     *	                  vxx,vzz,
     *					      x,y,t,p,pz,xx,yy,tt,pp,ppz,idd)
                  if(x-geophone_x.lt.0) then
                     dis(i_a)=1000000000000000.
                     cycle
                  end if
               s1=sqrt((x-geophone_x)**2+(y-geophone_y)**2)
     *				   /velocity_true+t
                  dis(i_a)=s1
                  xcoe1(i_a)=xx
                  ycoe1(i_a)=yy
                  tcoe1(i_a)=tt
                  if(ppz.gt.pp) ppz=pp
                  if(ppz.lt.pp.and.abs(ppz).gt.pp) ppz=-pp
                  pcoe1(i_a)=acos(ppz/pp)/pi*180
                  xcoe2(i_a)=x
                  ycoe2(i_a)=y
                  tcoe2(i_a)=t
                  if(pz.gt.pp) pz=pp
                  if(pz.lt.pp.and.abs(pz).gt.pp) pz=-pp
                  pcoe2(i_a)=acos(pz/pp)/pi*180
                  ang(i_a)=k/10.
               end do
               do m=1,i_a
                  if(dis(m).lt.s) then
                     k_min1=m
                     s=dis(m)
                  end if
               end do
               ss=1.0E+10
               do m=1,i_a
                  if(m.ne.k_min1) then
                     if(dis(m).lt.ss) then
                        k_min2=m
                        ss=dis(m)
                     end if
                  end if
               end do
               sas1=(tcoe2(k_min1)-tcoe2(k_min1-1))**2
     *	          +(pcoe2(k_min1)-pcoe2(k_min1-1))**2
               sas2=(tcoe2(k_min1)-tcoe2(k_min1+1))**2
     *             +(pcoe2(k_min1)-pcoe2(k_min1+1))**2
               if(sas1.gt.sas2) then
                  k_min2=k_min1+1
               else
                  k_min2=k_min1-1
               end if
               ttt2_1=tcoe1(k_min2)
               xxx2_1=xcoe1(k_min2)
               yyy2_1=ycoe1(k_min2)
               ppp2_1=pcoe1(k_min2)
               ttt2_2=tcoe2(k_min2)
               xxx2_2=xcoe2(k_min2)
               yyy2_2=ycoe2(k_min2)
               ppp2_2=pcoe2(k_min2)
               ttt1_1=tcoe1(k_min1)
               xxx1_1=xcoe1(k_min1)
               yyy1_1=ycoe1(k_min1)
               ppp1_1=pcoe1(k_min1)
               ttt1_2=tcoe2(k_min1)
               xxx1_2=xcoe2(k_min1)
               yyy1_2=ycoe2(k_min1)
               ppp1_2=pcoe2(k_min1)
               tt1=ttt1_1+(ttt1_2-ttt1_1)/(xxx1_2-xxx1_1)
     *			   *(geophone_x-xxx1_1)
               pp1=ppp1_1+(ppp1_2-ppp1_1)/(xxx1_2-xxx1_1)
     *			   *(geophone_x-xxx1_1)
               yy1=yyy1_1+(yyy1_2-yyy1_1)/(xxx1_2-xxx1_1)
     *			   *(geophone_x-xxx1_1)
               tt2=ttt2_1+(ttt2_2-ttt2_1)/(xxx2_2-xxx2_1)
     *			   *(geophone_x-xxx2_1)
               pp2=ppp2_1+(ppp2_2-ppp2_1)/(xxx2_2-xxx2_1)
     *			   *(geophone_x-xxx2_1)
               yy2=yyy2_1+(yyy2_2-yyy2_1)/(xxx2_2-xxx2_1)
     *			   *(geophone_x-xxx2_1)
               tt=tt1+(tt2-tt1)/(yy2-yy1)*(geophone_y-yy1)
               pp=pp1+(pp2-pp1)/(yy2-yy1)*(geophone_y-yy1)
               t1(i,j)=tt
               p1(i,j)=pp/cof
               ang1(i,j)=ang(k_min1)+(ang(k_min2)-ang(k_min1))
     *		           /(yy2-yy1)*(geophone_y-yy1)
               ang1(i,j)=ang1(i,j)/cof
               if(ijk.ne.1) then
                  sss=sss+(t1(i,j)-t0(i,j))**2+(p1(i,j)-p0(i,j))**2
               else
                  sss=sss+(p1(i,j)-p0(i,j))**2
                  sss=sss+(ang1(i,j)-ang0(i,j))**2
               end if
            end do
         end do
         if(sss.lt.sss_min) then
            ign=1
            sss_min=sss
            it_min=it
         end if
         write(16,*) it,sss
      end do
197   continue
      write(*,*) 'it_min',it_min,sss_min,ign,a_c
167   continue
      if(ign.eq.1) then
         i_k=0
         do iijj=2,2,4
            do i=1,igrid_r
               do j=1,igrid_v
                  qsc=1+(j-1)*beta
                  i_k=i_k+1
                  if(iijj.eq.1) then
                     vxx(i,j)=vvxx(i,j)
            vvxx(i,j)=vvxx(i,j)-((it_min-1)*1./10.)*af(i_k)
                  end if
                  if(iijj.eq.2) then
                     vzz(i,j)=vvzz(i,j)
            vvzz(i,j)=vvzz(i,j)-((it_min-1)*1./10.)*af(i_k)
                  end if
               end do
            end do
         end do
         do i=1,igrid_r
            do j=1,igrid_v
               write(16,*) ks,i,j,vvxx(i,j),vvzz(i,j)
               write(123,*) i,j,vvzz(i,j)
            end do
         end do
c        write(*,*) 'one iteration'
      else
         do i=1,igrid_r
            do j=1,igrid_v
               write(16,*) ks,i,j,vvxx(i,j),vvzz(i,j)
               write(123,*) i,j,vvzz(i,j)
            end do
         end do
c        write(*,*) 'one iteration'
      end if
c     write(*,*) 'ok,ok'
      if(ks .LT. 10) goto 198
      END PROGRAM first_program

ccccccccccccc    subroutine cccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ray1(x0,z0,x_e,z_e,nz,nx,dx,dz,xu,zu,qq,
     *   	      vvxx1,vvzz1,xxx,yyy,ttt,ppp,pppz,xxx2,yyy2,ttt2
     *           ,ppp2,pppz2,idd)

c 			call ray1(x_cord_begin,y_cord_begin,x_cord_end,y_cord_end,
c     *			             igrid_r,igrid_v,
c     *			             grid_interval_x,grid_interval_y,
c     *					     source_x,source_y,
c     *	                     angle,
c     *	                     vv,
c     *					     x,y,t,p,pz,xx,yy,tt,pp,ppz,idd)

      real vvxx1(nz,nx),vvzz1(nz,nx)
      real vvxx(nx,nz),vvzz(nx,nz)
      real x_ray(90000),z_ray(90000)
      real t(9000),pp(9000),ppx(9000),ppz(9000)
      double precision x,z,dt,p,px,pz,v,p1,x1,x2,z1,z2,tx,tz
     *       ,tv1,tv2,vv1,vv2,dvx,dvx1,dvx2,dvz,dvz1,dvz2
      do i=1,nx
         do j=1,nz
            vvxx(i,j)=vvxx1(j,i)
            vvzz(i,j)=vvzz1(j,i)
         end do
      end do
      num=1
      ds=1
      p=1000
      px=p*sin(qq)
      pz=p*cos(qq)
      xxu=xu
      zzu=zu
      x_ray(num)=xu
      z_ray(num)=zu
      pp(num)=p
      ppx(num)=px
      ppz(num)=pz

      do while((xxu.gt.x0.and.xxu.lt.x0+(nx-1)*dx).and.
     *   (zzu.gt.z0.and.zzu.lt.z0+(nz-1)*dz))
      num=num+1
      ix=(xxu/dx)+1
      iz=(zzu/dz)+1
      ks=ks+1
      vx1=vvxx(ix,iz)+(vvxx(ix+1,iz)-vvxx(ix,iz))/dx*(xxu-(ix-1)*dx)
      vx2=vvxx(ix,iz+1)+(vvxx(ix+1,iz+1)-vvxx(ix,iz+1))/dx
     * *(xxu-(ix-1)*dx)
      vx=vx1+(vx2-vx1)/dz*(zzu-(iz-1)*dz)
      vz1=vvzz(ix,iz)+(vvzz(ix+1,iz)-vvzz(ix,iz))/dx*(xxu-(ix-1)*dx)
      vz2=vvzz(ix,iz+1)+(vvzz(ix+1,iz+1)-vvzz(ix,iz+1))/dx
     * *(xxu-(ix-1)*dx)
      vz=vz1+(vz2-vz1)/dz*(zzu-(iz-1)*dz)
      dpx=-vx*ds
      dpz=-vz*ds
      px=px+dpx
      pz=pz+dpz
      if(pz.ne.0) qq=atan (px/pz)
      if(pz.eq.0) qq=pi/2.
      if(qq.lt.0) qq=qq+pi
      xxu=xxu+ds*sin(qq)
      zzu=zzu+ds*cos(qq)
      x_ray(num)=xxu
      z_ray(num)=zzu
      pp(num)=p
      ppx(num)=px
      ppz(num)=pz
      t(num)=1
c      if(xxu.le.x0.or.xxu.ge.x0+(nx-1)*dx) goto 11
c      if(zzu.le.z0.or.zzu.ge.z0+(nz-1)*dz) goto 11
c      goto 1
c11    continue
      end do
      xxx=x_ray(num)
      yyy=z_ray(num)
      ttt=t(num)
      ppp=pp(num)
      pppx=ppx(num)
      pppz=ppz(num)
      xxx2=x_ray(num-1)
      yyy2=z_ray(num-1)
      ttt2=t(num-1)
      ppp2=pp(num-1)
      pppx2=ppx(num-1)
      pppz2=ppz(num-1)
c      continue
      return
      end


      subroutine gaussj(a,n,np,b,m,mp)
      INTEGER m,mp,n,np,NMAX
      REAL a(np,np),b(np,mp)
      PARAMETER (NMAX=2200)
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      REAL big,dum,pivinv
      ipiv(1:n)=0
      do i=1,n
         big=0.
         do j=1,n
            if(ipiv(j).ne.1)then
               do k=1,n
                  if (ipiv(k).eq.0) then
                     if (abs(a(j,k)).ge.big)then
                        big=abs(a(j,k))
                        irow=j
                        icol=k
                     endif
                  else if (ipiv(k).gt.1) then
                     pause 'singular matrix in gaussj'
                  endif
                end do
            endif
         end do
         ipiv(icol)=ipiv(icol)+1
         if (irow.ne.icol) then
            do l=1,n
               dum=a(irow,l)
               a(irow,l)=a(icol,l)
               a(icol,l)=dum
            end do
            do l=1,m
               dum=b(irow,l)
               b(irow,l)=b(icol,l)
               b(icol,l)=dum
            end do
         endif
         indxr(i)=irow
         indxc(i)=icol
         if (a(icol,icol).eq.0.) pause 'singular matrix in gaussj'
         pivinv=1./a(icol,icol)
         a(icol,icol)=1.
         do l=1,n
            a(icol,l)=a(icol,l)*pivinv
         end do
         do l=1,m
            b(icol,l)=b(icol,l)*pivinv
         end do
         do ll=1,n
            if(ll.ne.icol)then
               dum=a(ll,icol)
               a(ll,icol)=0.
               do l=1,n
                  a(ll,l)=a(ll,l)-a(icol,l)*dum
               end do
               do l=1,m
                  b(ll,l)=b(ll,l)-b(icol,l)*dum
               end do
            endif
         end do
      end do
      do l=n,1,-1
         if(indxr(l).ne.indxc(l))then
            do k=1,n
               dum=a(k,indxr(l))
               a(k,indxr(l))=a(k,indxc(l))
               a(k,indxc(l))=dum
            end do
         endif
      end do
      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine smooth(v,nx,nz)
	   real v(nx,nz),vel(nx,nz)
	   c_sm=0.2
	   r=1
      do i=1,nx
	      do j=1,nz
	         s=0
	         k=0
	         do m=1,nx
	            do n=1,nz
	               if(abs(m-i).le.r.and.abs(n-j).le.r.and.(abs(m-i)+
     *               abs(n-j)).ne.0) then
				         s=s+c_sm*v(m,n)
	                  k=k+1
	               end if
	            end do
	         end do
	         vel(i,j)=(v(i,j)+s)/(k*c_sm+1)
	      end do
      end do
      do i=1,nx
	      do j=1,nz
	         v(i,j)=vel(i,j)
	      end do
      end do
	   return
      end

      subroutine rann(ix,yfl)
	   if(ix.eq.0) ix=67107
	   ix=125*ix
	   ix=ix-ix/2796203*2796203
	   yfl=float(ix)
	   yfl=yfl/2796203.
	   return
      end

ccccccccccccc    subroutine cccccccccccccccccccccccccccccccccccccccccccccc
     	subroutine ray11(x0,z0,x_e,z_e,nz,nx,dx,dz,xu,zu,qq,
     *   	   va,xxx,yyy,ttt,ppp,pppz,xxx2,yyy2,ttt2,ppp2,pppz2,idd)
	   real vel(nx,nz),x_ray(90000),z_ray(90000)
	   real t(9000),pp(9000),ppx(9000),ppz(9000),va(nz,nx)
	   double precision x,z,dt,p,px,pz,v,p1,x1,x2,z1,z2,tx,tz
     *       ,tv1,tv2,vv1,vv2,dvx,dvx1,dvx2,dvz,dvz1,dvz2

	   kkk=0
	   sign=-1
	   x=xu
	   z=zu

	   do i=1,nx
	      do j=1,nz
	         vel(i,j)=va(j,i)
	      end do
	   end do

	   if(idd.eq.-9) write(*,*) 'sub 1'
	   if(x.lt.x0.or.x.gt.x_e) goto 1
	   do i=2,nx
	      xn1=x0+(i-2)*dx
	      xn2=x0+(i-1)*dx
	      if(x.ge.xn1.and.x.lt.xn2) then
	         vn=vel(i-1,1)+(vel(i,1)-vel(i-1,1))/dx*(x-xn1)
         end if
	   end do
	   p=1./vn
	   px=p*sin(qq)
	   pz=p*cos(qq)
	   v=vn
	   dt=0.0005
      ix1=int((x-x0)/dx)+1
      izz=int(z/dz)+1
	   x1=x0+(ix1-1)*dx
	   x2=x1+dx
	   z1=z0+(izz-1)*dz
	   z2=z1+dz
      tx=(x-x1)/dx
	   tz=(z-z1)/dz
	   tv1=vel(ix1+1,izz)-vel(ix1,izz)
	   tv2=vel(ix1+1,izz+1)-vel(ix1,izz+1)
      vv1=vel(ix1,izz)+tx*tv1
	   vv2=vel(ix1,izz+1)+tx*tv2
	   v=vv1+tz*(vv2-vv1)
cccccccccccccccccccccccccccccccccccc
	   p=1./v
	   px=p*sin(qq)
	   pz=p*cos(qq)
ccccccccccccccccccccccccccccccccccc
	   dvx1=tv1/dx
	   dvx2=tv2/dx
	   dvx=dvx1+tz*(dvx2-dvx1)
      dvz1=(vel(ix1,izz+1)-vel(ix1,izz))/dz
	   dvz2=(vel(ix1+1,izz+1)-vel(ix1+1,izz))/dz
	   dvz=dvz1+tx*(dvz2-dvz1)
101	  nray_point=1
      x_ray(nray_point)=x
	   z_ray(nray_point)=z
	   t(nray_point)=0
	   pp(nray_point)=p
	   ppx(nray_point)=px
	   if(idd.eq.-9) write(*,*) 'sub2'
	   do k=1,1000000
	      sign=-1
	      x=x+dt*v*v*px
	      z=z+dt*v*v*pz
	      px=px-dt*p*p*v*dvx
	      pz=pz-dt*p*p*v*dvz
	      p1=px**2+pz**2
	      p1=sqrt(p1)
         ix1=int((x-x0)/dx)+1
         izz=int((z-z0)/dz)+1
	      x1=x0+(ix1-1)*dx
	      x2=x1+dx
	      z1=z0+(izz-1)*dz
	      z2=z1+dz
	      tx=(x-x1)/dx
	      tz=(z-z1)/dz
	      if(ix1.le.1) ix=1
	      if(ix1.ge.nx) ix1=nx-1
	      if(izz.le.1) izz=1
	      if(izz.ge.nz) izz=nz-1
	      tv1=vel(ix1+1,izz)-vel(ix1,izz)
c 	      write(*,*) 'pass sss4',ix1,izz,x,z

	      tv2=vel(ix1+1,izz+1)-vel(ix1,izz+1)
c	      write(*,*) 'pass sss3'

         vv1=vel(ix1,izz)+tx*tv1
	      vv2=vel(ix1,izz+1)+tx*tv2
	      v=vv1+tz*(vv2-vv1)
c	      write(*,*) 'pass sss1'

	      dvx1=tv1/dx
	      dvx2=tv2/dx
	      dvx=dvx1+tz*(dvx2-dvx1)
	      dvz1=(vel(ix1,izz+1)-vel(ix1,izz))/dz
	      dvz2=(vel(ix1+1,izz+1)-vel(ix1+1,izz))/dz
	      dvz=dvz1+tx*(dvz2-dvz1)
  	      if(v.lt.100) write(*,*) 'pass22-1',k,v,vv1,vv2,ix1,izz

	      p=1/v

         a=1.
	      px=px*a
	      pz=pz*a

	      nray_point=nray_point+1
	      x_ray(nray_point)=x
	      z_ray(nray_point)=z
	      t(nray_point)=k*dt
		   pp(nray_point)=p
	      ppx(nray_point)=px
	      ppz(nray_point)=pz
	      sign=1

		   if(x.le.x0.or.x.ge.x0+(nx-1)*dx) goto 1
	      if(z.le.z0.or.z.ge.z0+(nz-1)*dz) goto 1
	   end do
1     continue
      xxx=x_ray(nray_point)
      yyy=z_ray(nray_point)
	   ttt=t(nray_point)
	   ppp=pp(nray_point)
	   pppx=ppx(nray_point)
	   pppz=ppz(nray_point)
      xxx2=x_ray(nray_point-1)
      yyy2=z_ray(nray_point-1)
	   ttt2=t(nray_point-1)
	   ppp2=pp(nray_point-1)
	   pppx2=ppx(nray_point-1)
	   pppz2=ppz(nray_point-1)
      continue
      if(idd.eq.-9) write(*,*) 'sub 4'
      return
      end

