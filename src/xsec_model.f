      subroutine xsec_model(ispec,ebeam,pcent,thcent,delta,yptar,xptar,
     > xbj,xsecvert,xsecrad)

      implicit none

      real*8 ebeam,pcent,thcent
      real*8 delta,yptar,xptar
      real*8 plabx,plaby,plabz,pz,p
      real*8 costheta,sintheta
      real*8 thetalab,raddeg
      real*8 nu,wsq,q2,xbj,mp,qabs
      real*8 A, Z, m_tgt
      real*8 sig_dis, sig_qe,sigtot
      real*8 xsecvert,xsecrad,xdum,q2dum,edum
      real*8 sig_rad_func, sig_vert_func
      integer ispec,filelen,i,j,ncomment
      character*80 line,radfile,xfile
      logical first
      data first/.true./

      include 'radcor.cmn'

      radfile = 'rctables/lad_rc_d2cryo_lad.out'
      mp=0.938272
      raddeg = 180.0/3.14159

      p=(pcent/1000.0)*(1.0+delta/100.0)
      pz=p/sqrt(1.0+xptar**2+yptar**2)

      sintheta=sin(thcent)
      costheta=cos(thcent)
      if(ispec.eq.2) then
         plabx=pz*xptar
         plaby=pz*(yptar*costheta+sintheta)
         plabz=pz*(-yptar*sintheta+costheta)
      elseif(ispec.eq.1) then
         plabx=pz*xptar
         plaby=pz*(yptar*costheta-sintheta)
         plabz=pz*(yptar*sintheta+costheta)
      endif

      thetalab=acos(plabz/p)

      qabs=sqrt((plabx**2+plaby**2+(ebeam-plabz)**2))
      nu=Ebeam-p
      q2=qabs**2-nu**2
      wsq=-q2+mp**2+2.*mp*nu
      xbj=q2/2./mp/nu

      xsecrad=0.0
      xsecvert=0.0
      if(first) then ! read in table
         do i=32,1,-1
            if (radfile(i:i).eq.' ') filelen=i-1
         enddo

         if(filelen.gt.0) then
            write(6,*) 'Radiative corrections file in input file'
            xfile = radfile(1:filelen)
            write(6,*) "opening ",xfile
            open (unit=10,file=xfile,status='old')
            ncomment=1
            do i=1,ncomment
               read(10,'(a)',end=99) line
               if(line(1:3).ne.'***') then
                  write(6,*) 'Warning - expected comment line!!'
               endif
            enddo               !loop over comment lines
            kmax=275
            write(6,*) 'number of momentum bins - RAD1',kmax
            ntheta = 69         ! my tables
            do i = 1,ntheta	!theta bins
               do j=1,kmax      !loop over eprime dbins
                  read(10,1005,end=99) edum,eprad(j),thrad(i),xdum,q2dum,sigv(i,j),
     >  	    sigvdis(i,j),sigvqe(i,j),sigrad(i,j),sigrel(i,j),
     >              sigrqe(i,j),sigrdis(i,j),coulcor(i,j)
                  if(sigrad(i,j) .gt.0.0) then
                     radcor(i,j) = sigv(i,j)/sigrad(i,j)
                  else
                     radcor(i,j) =0.0
                  endif
                  write(9,*) 'radcor read in:',i,j,eprad(j),thrad(i),radcor(i,j),sigrad(i,j),sigv(i,j),coulcor(i,j)
               enddo		! loop over eprime
            enddo               !loop over theta
 99         continue
            close(10)
         endif
         first=.false.                     
      endif


      xsecrad=sig_rad_func(p,thetalab)
      xsecvert=sig_vert_func(p,thetalab)
      
 1005 format(1x,5f9.4,9e13.5)
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function sig_vert_func(ep,theta)
	
      implicit none

      real*8 ep,theta,thdeg,thtmp,eptmp
      real*8 thhi,thlo,delta_th
      real*8 ephi,eplo,delta_ep
      real*8 A1,A2,A3,A4,A12,A34,A1234
      real*8 B1,B2,B3,B4,B12,B34,B1234
      real*8 sig_vert_qe
      
      integer thcount,epcount

      include 'radcor.cmn'


      thdeg = theta*180.0/3.141592654
      thtmp = thdeg
      
      if(thdeg.lt.thrad(1)) then
         thtmp=thrad(1)
         write(6,*) 'WARNING, theta .lt. theta_min: sig_vert'
      endif

      if(thdeg.gt.thrad(ntheta)) then
         thtmp=thrad(ntheta)
         write(6,*) 'WARNING, theta .gt. theta_max: sig_vert',thdeg,thrad(ntheta)
      endif

      eptmp=ep
      
      if(ep.lt.eprad(1)) then
         eptmp = eprad(1)
         write(6,*) 'WARNING, EP .lt. ep_min: sig_vert',ep,eprad(1)
      endif

      if(ep.gt.eprad(kmax)) then
         eptmp = eprad(kmax)
         write(6,*) 'WARNING, EP .gt. ep_max: sig_vert'
      endif
      
      
      do thcount=1,ntheta-1
         thhi=0.0
         thlo=0.0
         ephi=0.0
         eplo=0.0
         if( (thtmp.gt.thrad(thcount)) .and. (thtmp.le.thrad(thcount+1)) ) then
         thhi=thrad(thcount+1)
         thlo=thrad(thcount)
         delta_th=thhi-thlo
         do epcount=1,kmax-1
            if( (eptmp.gt.eprad(epcount)) .and. (eptmp.le.eprad(epcount+1)) )then
            ephi=eprad(epcount+1)
            eplo=eprad(epcount)
            delta_ep=ephi-eplo

            A1 = sigv(thcount,epcount)
            A2 = sigv(thcount+1,epcount)
            A3 = sigv(thcount,epcount+1)
            A4 = sigv(thcount+1,epcount+1)
            
            A12 = (A1*(thhi-thtmp) + A2*(thtmp-thlo))/delta_th
            A34 = (A3*(thhi-thtmp) + A4*(thtmp-thlo))/delta_th
            
            A1234 = (A12*(ephi-eptmp) + A34*(eptmp-eplo))/delta_ep

            B1 = sigvqe(thcount,epcount)
            B2 = sigvqe(thcount+1,epcount)
            B3 = sigvqe(thcount,epcount+1)
            B4 = sigvqe(thcount+1,epcount+1)
            
            B12 = (B1*(thhi-thtmp) + B2*(thtmp-thlo))/delta_th
            B34 = (B3*(thhi-thtmp) + B4*(thtmp-thlo))/delta_th
            
            B1234 = (B12*(ephi-eptmp) + B34*(eptmp-eplo))/delta_ep
         endif                  !ep check
      enddo                     !loop over ep
      endif                     !theta check
      enddo                     !loop over theta

      sig_vert_func=A1234
      sig_vert_func=sig_vert_func/1000.0 ! convert ot ub I think
c      sig_vert_qe=B1234/1000.0

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC This function calculates the radiated cross section from the "main" 
CC 2D radcor file
	real*8 function sig_rad_func(ep,theta)
	
	implicit none

	real*8 ep,theta,thdeg,thtmp,eptmp
	real*8 thhi,thlo,delta_th
	real*8 ephi,eplo,delta_ep
	real*8 A1,A2,A3,A4,A12,A34,A1234

	integer thcount,epcount

	include 'radcor.cmn'


	thdeg = theta*180.0/3.141592654

	thtmp=thdeg

	if(thdeg.lt.thrad(1)) then
	   thtmp=thrad(1)
	   write(6,*) 'WARNING, theta .lt. theta_min: sig_rad',thdeg,thrad(1)
	   return
	endif

	if(thdeg.gt.thrad(ntheta)) then
	   thtmp=thrad(ntheta)
	   write(6,*) 'WARNING, theta .gt. theta_max: sig_rad',thdeg,thrad(ntheta)
	   return
	endif

	eptmp=ep

	if(ep.lt.eprad(1)) then
	   eptmp = eprad(1)
c	   write(6,*) 'WARNING, EP .lt. ep_min: sig_rad'
	endif

	if(ep.gt.eprad(kmax)) then
c	   write(6,*) 'WARNING, EP .gt. ep_max: sig_rad',eptmp,eprad(kmax)
	   eptmp = eprad(kmax)
	endif


	do thcount=1,ntheta-1
	   thhi=0.0
	   thlo=0.0
	   ephi=0.0
	   eplo=0.0
	   if( (thtmp.gt.thrad(thcount)) .and. (thtmp.le.thrad(thcount+1)) ) then
	      thhi=thrad(thcount+1)
	      thlo=thrad(thcount)
	      delta_th=thhi-thlo
	      do epcount=1,kmax-1
		 if( (eptmp.gt.eprad(epcount)) .and. (eptmp.le.eprad(epcount+1)) )then
		    ephi=eprad(epcount+1)
		    eplo=eprad(epcount)
		    delta_ep=ephi-eplo

		    A1 = sigrad(thcount,epcount)
		    A2 = sigrad(thcount+1,epcount)
		    A3 = sigrad(thcount,epcount+1)
		    A4 = sigrad(thcount+1,epcount+1)

C This gives the radiated model with "Coulombic effects" - this should match the
C uncorrected data better and give better bin centering corrections

c		    A1 = sigrad(thcount,epcount)/coulcor(thcount,epcount)
c		    A2 = sigrad(thcount+1,epcount)/coulcor(thcount+1,epcount)
c		    A3 = sigrad(thcount,epcount+1)/coulcor(thcount,epcount+1)
c		    A4 = sigrad(thcount+1,epcount+1)/coulcor(thcount+1,epcount+1)

		    A12 = (A1*(thhi-thtmp) + A2*(thtmp-thlo))/delta_th
		    A34 = (A3*(thhi-thtmp) + A4*(thtmp-thlo))/delta_th

		    A1234 = (A12*(ephi-eptmp) + A34*(eptmp-eplo))/delta_ep
		 endif   !ep check
	      enddo      !loop over ep
	   endif    !theta check
	enddo       !loop over theta

	sig_rad_func=A1234
c	write(6,*) 'sig_rad_func 1', sig_rad_func
	sig_rad_func=sig_rad_func/1000.0
	if(sig_rad_func.lt.0.0) then
	   sig_rad_func=0.0
	endif
c	sig_rad_func=1.0
c	write(6,*) 'sig_rad_func 2', sig_rad_func

	return
	end
