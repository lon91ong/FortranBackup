! module

module subfunc
	implicit real*4(a-h,o-z)
	implicit integer(i-n)
	!save
	real(4),parameter::pi=3.14159265
	real(4),parameter :: ep = epsilon(pi)

	parameter(num_s=11,num_g=11,igrid_r=11,igrid_v=11)
	parameter(igrid_r_1=11,igrid_v_1=11)
	real t0(num_s,num_g),t1(num_s,num_g),t2(num_s,num_g)
	real dt(num_s,num_g)
	real v(igrid_r,igrid_v),et(2*num_s*num_g,igrid_r*igrid_v)
	real aa(igrid_r*igrid_v,igrid_r*igrid_v),b(igrid_r*igrid_v)
	real aet(igrid_r*igrid_v,2*num_s*num_g),xc(igrid_r,igrid_v)
	real f(2*num_s*num_g),af(igrid_r*igrid_v)
	real vv(igrid_r,igrid_v),vt(igrid_r,igrid_v)
	real xcoe1(10000),ycoe1(10000),xcoe2(10000),ycoe2(10000)
	real tcoe1(10000),tcoe2(10000),dis(10000),ang(10000)
	real pcoe1(10000),pcoe2(10000),p0(num_s,num_g)
	real p1(num_s,num_g),p2(num_s,num_g)
	real zn(1000),zni(1000),ang1(num_s,num_g)
	real v_1(igrid_r_1,igrid_v_1)
	real x_ray(90000),z_ray(90000)
	real ray_x(num_s*num_g,90000),ray_z(num_s*num_g,90000)
	real x_x(4),y_y(4),s_s(4)

	contains

	subroutine ray1(x0,z0,x_e,z_e,nz,nx,dx,dz,xu,zu,qq,&
		va,xxx,yyy,ttt,ppp,pppz,xxx2,yyy2,ttt2,ppp2,pppz2)
		real vel(nx,nz),x_ray(90000),z_ray(90000)
		real t(9000),pp(9000),ppx(9000),ppz(9000),va(nz,nx)
		double precision x,z,dt,p,px,pz,v,p1,x1,x2,z1,z2,tx,tz,&
		tv1,tv2,vv1,vv2,dvx,dvx1,dvx2,dvz,dvz1,dvz2

		dt=0.0005
		kkk=0
		sign=-1
		x=xu
		z=zu

		do i=1,nx; do j=1,nz
			vel(i,j)=va(j,i)
		end do; end do

		if(x.ge.x0.and.x.le.x_e) then
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
			p=1./v
			px=p*sin(qq)
			pz=p*cos(qq)
			dvx1=tv1/dx
			dvx2=tv2/dx
			dvx=dvx1+tz*(dvx2-dvx1)
			dvz1=(vel(ix1,izz+1)-vel(ix1,izz))/dz
			dvz2=(vel(ix1+1,izz+1)-vel(ix1+1,izz))/dz
			dvz=dvz1+tx*(dvz2-dvz1)
			nray_point=1
			x_ray(nray_point)=x
			z_ray(nray_point)=z
			t(nray_point)=0
			pp(nray_point)=p
			ppx(nray_point)=px

			k=1
			do while((x.ge.x0.and.x.le.x0+(nx-1)*dx).and.(z.ge.z0.and.z.le.z0+(nz-1)*dz))
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
				tv1=vel(ix1+1,izz)-vel(ix1,izz)
				tv2=vel(ix1+1,izz+1)-vel(ix1,izz+1)
				vv1=vel(ix1,izz)+tx*tv1
				vv2=vel(ix1,izz+1)+tx*tv2
				v=vv1+tz*(vv2-vv1)
				dvx1=tv1/dx
				dvx2=tv2/dx
				dvx=dvx1+tz*(dvx2-dvx1)
				dvz1=(vel(ix1,izz+1)-vel(ix1,izz))/dz
				dvz2=(vel(ix1+1,izz+1)-vel(ix1+1,izz))/dz
				dvz=dvz1+tx*(dvz2-dvz1)
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
				k = k+1
			end do
		end if
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

		return
	end subroutine ray1

	subroutine gaussj(a,n,np,b,m,mp)
		integer m,mp,n,np,NMAX
		real a(np,np),b(np,mp)
		parameter (NMAX=2200)
		integer i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
		real big,dum,pivinv
		do j=1,n
			ipiv(j)=0
		end do
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
	end subroutine gaussj

	subroutine rann(ix,yfl)
		if(ix.eq.0) ix=67107
		ix=125*ix
		ix=ix-ix/2796203*2796203
		yfl=float(ix)/2796203.
		!yfl=yfl/2796203.
		return
	end subroutine rann

	subroutine ray18(x0,z0,x_e,z_e,nz,nx,dx,dz,xu,zu,qq,va,x_ray,z_ray,nray_point)
		real vel(nx,nz),x_ray(90000),z_ray(90000)
		real t(9000),pp(9000),ppx(9000),ppz(9000),va(nz,nx)
		double precision x,z,dt,p,px,pz,v,p1,x1,x2,z1,z2,tx,tz,&
		tv1,tv2,vv1,vv2,dvx,dvx1,dvx2,dvz,dvz1,dvz2

		dt=0.0005
		kkk=0
		sign=-1
		x=xu
		z=zu

		do i=1,nx; do j=1,nz
			vel(i,j)=va(j,i)
		end do;	end do

		do while(x.ge.x0 .and. x.le.x_e)
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
			p=1./v
			px=p*sin(qq)
			pz=p*cos(qq)
			dvx1=tv1/dx
			dvx2=tv2/dx
			dvx=dvx1+tz*(dvx2-dvx1)
			dvz1=(vel(ix1,izz+1)-vel(ix1,izz))/dz
			dvz2=(vel(ix1+1,izz+1)-vel(ix1+1,izz))/dz
			dvz=dvz1+tx*(dvz2-dvz1)
			nray_point=1
			x_ray(nray_point)=x
			z_ray(nray_point)=z
			t(nray_point)=0
			pp(nray_point)=p
			ppx(nray_point)=px

			k=1
			do while((x.ge.x0.and.x.le.x0+(nx-1)*dx).and.(z.ge.z0.and.z.le.z0+(nz-1)*dz))
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
				tv1=vel(ix1+1,izz)-vel(ix1,izz)
				tv2=vel(ix1+1,izz+1)-vel(ix1,izz+1)
				vv1=vel(ix1,izz)+tx*tv1
				vv2=vel(ix1,izz+1)+tx*tv2
				v=vv1+tz*(vv2-vv1)
				dvx1=tv1/dx
				dvx2=tv2/dx
				dvx=dvx1+tz*(dvx2-dvx1)
				dvz1=(vel(ix1,izz+1)-vel(ix1,izz))/dz
				dvz2=(vel(ix1+1,izz+1)-vel(ix1+1,izz))/dz
				dvz=dvz1+tx*(dvz2-dvz1)
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
				k = k+1
			end do
		end do

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

		return
	end subroutine ray18
end module
