!main.f90

program main
	use subfunc
	implicit real*4(a-h,o-z)
	implicit integer(i-n)

	d_x=200
	d_y=400
	grid_interval_x=200./(igrid_v-1)    !200./(11-1)=20.
	grid_interval_y=400./(igrid_r-1)    !400./(11-1)=40.
	grid_interval_x_1=200./(igrid_v_1-1)!200./(11-1)=20.
	grid_interval_y_1=400./(igrid_r_1-1)!400./(11-1)=40.
	interval_s=370./(num_s-1)   !370./(11-1)=37.
	interval_g=380./(num_g-1)   !380./(11-1)=38.
	y_first_s=20
	y_first_g=10
	source_x=0 !Ô´
	geophone_x=200 !½ÓÊÕÆ÷
	velocity_true=2000
	velocity_initial=2000
	x_cord_begin=0
	x_cord_end=200
	y_cord_begin=0
	y_cord_end=400
	idt_scan=20
	a_c=1/1.
	itt=5
	ijk=1
	cof=1.
	!cccccccccccccccc 'seet-t.dat' is the file of result  ccccccccccccccccc
	open(123,file='see-t.dat')
	!open(16,file='true-2-t.dat')
	!open(17,file='c_saa-2-1-for-10-t.dat')
	!open(18,file='curve-test1.dat')
	!open(19,file='curve-test2.dat')

	idd=0
	ttmax=0
	ppmax=0
	ttmin=100000000.
	ppmin=100000000.
	st=0
	sp=0
	ix=2
	rant=0.00
	ranp=0.
	!      cccccc  2. define true velcoity model  ccccccccccccccc

    v_1(:,:)=velocity_true
    v_1((igrid_r_1+1)/2:,:)=velocity_true+100

	do i=1,igrid_r_1
		do j=1,igrid_v_1
			write(123,*) i,j,v_1(i,j)
			write(*,*) i,j,v_1(i,j)
		end do
	end do

	!cccccccc     3. get the true traveltime in geophone ccccccccccc
	write(*,*) 'begin to get the true traveltime'
	do i=1,num_s; do j=1,num_g
		source_y=y_first_s+(i-1)*interval_s
		geophone_y=y_first_g+(j-1)*interval_g
		s=1.0E+10
		if((geophone_y-source_y).ne.0) then
			a1=atan((x_cord_end-x_cord_begin)/(geophone_y-source_y))
			a1=a1/pi*180.
			ia1=int(a1)
			if(ia1.lt.0) ia1=ia1+180
		else
			ia1=90
		end if
		jk1 = max(100,(ia1-idt_scan)*10)
		jk2 = min(1700,(ia1+idt_scan)*10)
		i_a=0
		do k=jk1,jk2,itt
			i_a=i_a+1
			angle=pi/180.*(k*0.1)

			call ray1(x_cord_begin,y_cord_begin,x_cord_end,y_cord_end, &
			igrid_r_1,igrid_v_1, &
			grid_interval_x_1,grid_interval_y_1, &
			source_x,source_y, angle, v_1, &
			x,y,t,p,pz,xx,yy,tt,pp,ppz)

			if(x-geophone_x.lt.0) then
				dis(i_a)=1000000000000000.
			else
				s1=sqrt((x-geophone_x)**2+(y-geophone_y)**2)/velocity_true+t
				dis(i_a)=s1
				xcoe1(i_a)=xx
				ycoe1(i_a)=yy
				tcoe1(i_a)=tt
				ppz=min(pp,ppz)
				if(ppz.lt.pp.and.abs(ppz).gt.pp) ppz=-pp
				pcoe1(i_a)=acos(ppz/pp)/pi*180
				xcoe2(i_a)=x
				ycoe2(i_a)=y
				tcoe2(i_a)=t
				pz=min(pp,pz)
				if(pz.lt.pp.and.abs(pz).gt.pp) pz=-pp
				pcoe2(i_a)=acos(pz/pp)/pi*180
				ang(i_a)=k/10.
			end if
		end do

		k_min1=minloc(dis(1:i_a),1)
		s=minval(dis(1:i_a),1)

		ss=1.0E+10
		do m=1,i_a
			if(m.ne.k_min1 .and. dis(m).lt.ss) then
				k_min2=m
				ss=dis(m)
			end if
		end do

		sas1=(tcoe2(k_min1)-tcoe2(k_min1-1))**2+(pcoe2(k_min1)-pcoe2(k_min1-1))**2  !???
		sas2=(tcoe2(k_min1)-tcoe2(k_min1+1))**2+(pcoe2(k_min1)-pcoe2(k_min1+1))**2
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

		ttmax=max(abs(tt),ttmax)
		ppmax=max(abs(pp),ppmax)
		ttmin=min(abs(tt),ttmin)
		ppmin=min(abs(pp),ppmin)

		st=st+abs(tt)
		sp=sp+abs(pp)
		t0(i,j)=tt
		p0(i,j)=pp/cof
	end do; end do

	write(*,*) 'qqq',sp/st,ttmax,ppmax,ttmin,ppmin,ppmax/ttmax
	do i=1,num_s; do j=1,num_g
		call rann(ix,yfl)
		ran_t=(yfl-0.5)*2*rant
		t0(i,j)=t0(i,j)+ran_t
		call rann(ix,yfl)

		ran_p=(yfl-0.5)*2*ranp
		p0(i,j)=p0(i,j)+ran_p
	end do; end do


	!cccccccccc    4. define initial model cccccccccccccc

	write(*,*) 'begin to calculate the traveltime for initial model'
	ii_k=0
	do ii=1,igrid_r; do jj=1,igrid_v
		call rann(ix,yfl)
		ran_v=(yfl-0.5)*0
		c_x=(jj-1)*grid_interval_x
		c_y=(ii-1)*grid_interval_y

		v(ii,jj)=velocity_initial+ran_v*0   !???
		write(123,*) 'initial',ii,jj,v(ii,jj)
	end do; end do

	!ccccccccccccccccccc  5.  get time and path in initial model     cccccccccccccccccccc
	ks=0
	sss_min=0
	ign=0
	do while (sss_min.eq.0 .or. sss_min.gt.2*ep)
		ks = ks+1
		if(ign.eq.-1) a_c=a_c*0.1

		!write(17,*) 'aaaa',ks,it_min,sss_min,ign,a_c
		!do i=1,igrid_r; do j=1,igrid_v
		!    write(17,*) (j-1)*grid_interval_x,(i-1)*grid_interval_y,v(i,j)
		!end do; end do

		write(*,*) 'the number of iteration is ',ks
		ii_p=0
		do i=1,num_s; do j=1,num_g
			ii_p=ii_p+1
			source_y=y_first_s+(i-1)*interval_s
			geophone_y=y_first_g+(j-1)*interval_g
			s=1.0E+10
			if((geophone_y-source_y).ne.0) then
				a1=atan((x_cord_end-x_cord_begin)/(geophone_y-source_y))
				a1=a1/pi*180.
				ia1=int(a1)
				if(ia1.lt.0) ia1=ia1+180
			else
				ia1=90
			end if
			jk1 = max(100,(ia1-idt_scan)*10)
			jk2 = min(1700,(ia1+idt_scan)*10)
			i_a=0
			do k=jk1,jk2,itt
				i_a=i_a+1
				angle=pi/180.*(k/10.)

				call ray1(x_cord_begin,y_cord_begin,x_cord_end,y_cord_end, &
				igrid_r,igrid_v, grid_interval_x,grid_interval_y, &
				source_x,source_y, angle, v, &
				x,y,t,p,pz,xx,yy,tt,pp,ppz)

				if(x-geophone_x.lt.0) then
					dis(i_a)=1000000000000000.
				else
					s1=sqrt((x-geophone_x)**2+(y-geophone_y)**2)/velocity_true+t
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
				end if
			end do

            k_min1=minloc(dis(1:i_a),1)
            s=minval(dis(1:i_a),1)
			ss=1.0E+10
			do m=1,i_a
				if(m.ne.k_min1) then
					if(dis(m).lt.ss) then
						k_min2=m
						ss=dis(m)
					end if
				end if
			end do
			sas1=(tcoe2(k_min1)-tcoe2(k_min1-1))**2 &
			+(pcoe2(k_min1)-pcoe2(k_min1-1))**2
			sas2=(tcoe2(k_min1)-tcoe2(k_min1+1))**2 &
			+(pcoe2(k_min1)-pcoe2(k_min1+1))**2
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

			if(ks.eq.1)then
                if(ijk.eq.1) sss_min=sss_min+(t1(i,j)-t0(i,j))**2
                if(ijk.eq.2) &
                    sss_min=sss_min+(t1(i,j)-t0(i,j))**2+(p1(i,j)-p0(i,j))**2
            end if

			if(abs(yy2-yy1).lt.0.00001) then
				ang1(i,j)=ang(k_min1)
			else
				ang1(i,j)=ang(k_min1)+(ang(k_min2)-ang(k_min1)) &
				/(yy2-yy1)*(geophone_y-yy1)
			end if
			ang1(i,j)=ang1(i,j)/cof
			!cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
			gg=ang1(i,j)
			gg=gg/180.*pi
			call ray18(x_cord_begin,y_cord_begin,x_cord_end,y_cord_end, &
			igrid_r,igrid_v, grid_interval_x,grid_interval_y, &
			source_x,source_y, gg, v, x_ray,z_ray,ionum)

			do io=1,ionum+1
				if(io.eq.1) then
					ray_x(ii_p,io)=ionum
					ray_z(ii_p,io)=ionum
				else
					ray_x(ii_p,io)=x_ray(io-1)
					ray_z(ii_p,io)=z_ray(io-1)
				end if
			end do
		end do; end do

		!cccccccccccccccc  6. calculating the derivation    cccccccccccccccccccc
		write(*,*) 'begin to calculate the deviation'

		ii_k=0
		do ii=1,igrid_r; do jj=1,igrid_v
			ia=ia+1
			ii_k=ii_k+1
			i_k=0
			do i=1,num_s; do j=1,num_g
				i_k=i_k+1
				x0=(jj-1)*grid_interval_x
				y0=(ii-1)*grid_interval_y
				x_x(1)=x0
				y_y(1)=y0-grid_interval_y
				x_x(2)=x0-grid_interval_x
				y_y(2)=y0-grid_interval_y
				x_x(3)=x0-grid_interval_x
				y_y(3)=y0
				x_x(4)=x0
				y_y(4)=y0

				do m=1,4
					s_s(m)=1
					if(x_x(m).lt.0.or.y_y(m).lt.0 .or. x_x(m).gt.d_x-grid_interval_x .or. &
					y_y(m).gt.d_y-grid_interval_y) s_s(m)=-1
				end do

				s=0
				do m=1,4
					if(s_s(m).eq.1) then
						x1=x_x(m)
						y1=y_y(m)
						x2=x1+grid_interval_x
						y2=y1
						x3=x1
						y3=y1+grid_interval_y
						x4=x3+grid_interval_x
						y4=y3
						ionum=ray_x(i_k,1)

						i_num=0
						do n=2,ionum+1
							x=ray_x(i_k,n)
							y=ray_z(i_k,n)
							if(x.gt.x1.and.x.le.x2 .and. y.gt.y1.and.y.le.y3) then
								select case(m)
								case(1)
									s=s+(x4-x)*(y-y1)/(y3-y1)/(x4-x3)
								case(2)
									s=s+(x-x3)*(y-y1)/(y3-y1)/(x4-x3)
								case(3)
									s=s+(x-x1)/(x2-x1)+(x1-x)*(y-y1) &
									*(x4-x1)/(y3-y1)/(x2-x1)/(x4-x3)
								case(4)
									s=s+1+(x1-x)/(x2-x1)+(y1-y)/(y3-y1) &
									+(x-x1)*(y-y1)/(y3-y1)/(x2-x1)
								endselect
							end if
						end do
					end if
				end do
				f(i_k)=t1(i,j)-t0(i,j)
				et(i_k,ii_k)=s*(-1./v(ii,jj)/v(ii,jj))
			end do; end do
		end do; end do

		!ccccccccccccccccccccccc    7. inversion   cccccccccccccccccccc
		write(*,*) 'begin to inverse '
		tm_max=0
		do i=1,igrid_r*igrid_v; do j=1,igrid_r*igrid_v
			aa(i,j)=0
			do ii=1,ijk*num_s*num_g
				aa(i,j)=aa(i,j)+et(ii,i)*et(ii,j)
			end do
			if(abs(aa(i,j)).gt.tm_max) tm_max=abs(aa(i,j))
		end do; end do

		write(*,*) tm_max
		do i=1,igrid_r*igrid_v; do j=1,igrid_r*igrid_v
			if(i.eq.j) aa(i,j)=aa(i,j)+10*tm_max
		end do; end do

		do i=1,igrid_r*igrid_v
			b(i)=1
		end do
		call gaussj(aa,igrid_r*igrid_v,igrid_r*igrid_v,b,1,1)

		do i=1,igrid_r*igrid_v; do j=1,ijk*num_s*num_g
			aet(i,j)=0
			do k=1,igrid_r*igrid_v
				aet(i,j)=aet(i,j)+aa(i,k)*et(j,k)
			end do
		end do; end do

		do i=1,igrid_r*igrid_v
			af(i)=0
			do j=1,ijk*num_s*num_g
				af(i)=af(i)+aet(i,j)*f(j)
			end do
		end do

		!cccccccccccc  8.one dimension searching cccccccccccccccccccccccccccccccc
		write(*,*) 'begin to one dimension searching'

		ign=-1
		it_min=0.
		sss_min=10000000000.
		do it=2,100,2
			write(*,*) 'it=',it
			sss=0
			i_k=0
			do i=1,igrid_r; do j=1,igrid_v
				i_k=i_k+1
				!xzf=((it-1)*1./100)*af(i_k)
				!vt(i,j)=v(i,j)-xzf
				!if(vt(i,j).gt.4000) then
				!	vt(i,j)=4000
				!end if
				!if(vt(i,j).lt.1000) then
				!	vt(i,j)=1000
				!end if
				vt(i,j)=max(1000.,min(v(i,j)-((it-1)*1./100)*af(i_k),4000.))
			end do; end do

			do i=1,num_s; do j=1,num_g
				source_y=y_first_s+(i-1)*interval_s
				geophone_y=y_first_g+(j-1)*interval_g
				s=1.0E+10
				if((geophone_y-source_y).ne.0) then
					a1=atan((x_cord_end-x_cord_begin)/(geophone_y-source_y))
					a1=a1/pi*180.
					ia1=int(a1)
					if(ia1.lt.0) ia1=ia1+180
				else
					ia1=90
				end if
				jk1=(ia1-idt_scan)*10
				jk2=(ia1+idt_scan)*10
				if(jk1.lt.100) jk1=100
				if(jk2.gt.1700) jk2=1700
				i_a=0
				do k=jk1,jk2,itt
					i_a=i_a+1

					angle=pi/180.*(k/10.)
					call ray1(x_cord_begin,y_cord_begin,x_cord_end,y_cord_end, &
					igrid_r,igrid_v, &
					grid_interval_x,grid_interval_y, &
					source_x,source_y, angle, vt, &
					x,y,t,p,pz,xx,yy,tt,pp,ppz)

					if(x-geophone_x.lt.0) then
						dis(i_a)=1000000000000000.
					else
						s1=sqrt((x-geophone_x)**2+(y-geophone_y)**2)/velocity_true+t
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
					end if
				end do
				k_min1=minloc(dis(1:i_a),1)
				s=minval(dis(1:i_a),1)
				ss=1.0E+10
				do m=1,i_a
					if(m.ne.k_min1 .and. dis(m).lt.ss) then
						k_min2=m
						ss=dis(m)
					end if
				end do
				sas1=(tcoe2(k_min1)-tcoe2(k_min1-1))**2+(pcoe2(k_min1)-pcoe2(k_min1-1))**2
				sas2=(tcoe2(k_min1)-tcoe2(k_min1+1))**2+(pcoe2(k_min1)-pcoe2(k_min1+1))**2
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
				if(ijk.ne.1) then
					sss=sss+(t1(i,j)-t0(i,j))**2+(p1(i,j)-p0(i,j))**2
				else
					sss=sss+(t1(i,j)-t0(i,j))**2
				end if
			end do; end do

			if(sss.lt.sss_min) then
				ign=1
				sss_min=sss
				it_min=it
			end if
			write(*,*) it,sss
		end do
		write(*,*) 'it_min',ks,it_min,sss_min,ign,a_c

		!cccccccccccccc 9. data output  ccccccccccccccc
		if(ign.eq.1) then
			i_k=0
			do i=1,igrid_r; do j=1,igrid_v
				i_k=i_k+1
				xc(i,j)=v(i,j)
				v(i,j)=v(i,j)-((it_min-1)*1./100.)*af(i_k)
				!write(19,*) i_k,v(i,j)
			end do; end do

			do i=1,igrid_r; do j=1,igrid_v
				write(*,*) 'vel',ks,i,j,xc(i,j),v(i,j)
				write(123,*) i,j,v(i,j)
			end do;end do
			write(*,*) 'one iteration'
		else
			i_k=0
			it_min=10.
			do i=1,igrid_r; do j=1,igrid_v
				i_k=i_k+1
				xc(i,j)=v(i,j)
				v(i,j)=v(i,j)-(0./10.)*af(i_k)
			end do; end do

			do i=1,igrid_r; do j=1,igrid_v
				write(*,*) 'vel',ks,i,j,xc(i,j),v(i,j)
				write(123,*) i,j,v(i,j)
			end do; end do
			write(*,*) 'one iteration'
		end if
		write(*,*) 'ok,ok'
	end do
end program
