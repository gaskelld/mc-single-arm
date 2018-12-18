	subroutine bdisnew4he3(eic1,ep1,theta1,aa1,zz1,ann1,esep1,pmax1,
     &	innt1,innp1,f01,bigB1,ag1,bg1,alpha1,sigdeep)
        
!	implicit real*8 (a-z)
	implicit none
	real*8 wp,ww1,ww2, w11, w22, nu,p
	real*8 eic1,ep1,theta1,aa1,zz1,ann1,esep1,pmax1
        real*8  innt1,innp1,ag1,bg1,cg1,dg1,sigdeep
	real*8 f01,bigB1,alpha1,aa,ag,agd,akf
	real*8 alpha,ama,amn,amp,amp2,ann,anu,anup
	real*8 arg1,arg2,arg21,arg22,ama_1,bg,bgd,BigB,cgd
	real*8 dgd,dp,du,ei,eic,ep,esep,f0,fackf,innp,innt
	real*8 pi,pmax,q2,q3,radical,s2,sigm,thec,theta,tt2,u,w1
        real*8 w1a,w1an,w1an1,w1ap,w1ap1,w2,w2a,w2an,w2an1
        real*8 w2ap,w2ap1,weight1,weight2,wp_max,wpmin2,x,zz
!	modified version of BWF's routine bdis
!	ag,bg,cg,dg are coeff of a 2 gaussian fit to the f(y)
!	which was converted to n(k)
	real*8 rk,rho
	real*8 rq2		!to pass to w1w2

	integer ip,inp,iu
!	common/interp/npts,nterms,rk(200),rho(200)
	
!	eic = real(eic1)
!	ep = real(ep1)
!	theta = real(theta1)
!	aa = real(aa1)
!	zz = real(zz1)
!	ann = real(ann1)
!	esep = real(esep1)
!	pmax = real(pmax1)
!	innt = real(innt1)
!	innp = real(innp1)
!	ag = real(ag1)
!	bg = real(bg1)
!	f0 = real (f01)
!	BigB = real (bigB1)
!	alpha=real(alpha1)
	eic = eic1
	ep = ep1
	theta = theta1
	aa = aa1
	zz = zz1
	ann = ann1
	esep = esep1
	pmax = pmax1
	innt = innt1
	innp = innp1
	ag = ag1
	bg = bg1
	f0 = f01
	BigB = bigB1
	alpha=alpha1
	
	f0=f0*1000
	bigB=bigB*1000
	alpha=alpha/1000
	ag=ag*1000
	bg=bg*1000

	
        if (aa .eq. 3)then       !he3
	  agd =    0.39167E+01
	  bgd =    0.78468E+02
	  cgd =    0.17924E+00
	  dgd =    0.66511E+01
        elseif (aa .eq. 4)then	!he4
	  agd =   0.35160E+01
	  bgd =   0.60049E+02
	  cgd =   0.38390E+00
	  dgd =   0.13421E+02
        elseif ((aa .eq. 12).or.(aa.eq.9)) then      !carbon
	  agd =   0.28757E+01
	  bgd =   0.41922E+02
	  cgd =   0.33801E+00
	  dgd =   0.13824E+02
        elseif (aa .eq. 27)then      !alum
	  agd =   0.25660E+01
	  bgd =   0.36962E+02
	  cgd =   0.41675E+00
	  dgd =   0.13772E+02
        elseif ((aa .eq. 56).or.(aa.eq.64))then      !iron
	  agd =   0.28316E+01
	  bgd =    0.44624E+02
	  cgd =    0.37850E+00
	  dgd =    0.12003E+02
        elseif(aa .eq. 197)then	!gold, Jerry, Gold!
	  agd =   0.24947E+01
	  bgd =   0.30614E+02
	  cgd =   0.42398E+00
	  dgd =   0.12272E+02
c	else
c	   write(6,*) 'cannot determine nucleus - crap!'
c	   stop
        endif

!       temp values for check.

!	ag =   0.75287E+01
!	bg =   0.18856E+03
!	cg =   0.22657E+00
!	dg =   0.14649E+02
!	eic incident energy in gev
!	ep final energy
!	thec scattering angle in radians
!	aa atomic number
!	zz number of protons
!	ann number of neutrons
!	esep spearation energy
!	pmax max p in integration of p > kf
!	innt integration steps over theta
!	integration steps over p
!       
!       
!	write(6,*)'E,EP,THETA, A, ZZ,ANN,ESP,PMAX,INNT,INNp,ag,bg,bigB, f0, alpha'
!	write(6,*)eic,ep,theta,aa,zz,ann,esep,pmax,innt,innp,ag,bg,bigB,
!	1 f0, alpha

	pi = 3.14159265
	thec=theta
	theta=theta*pi/180
	s2 = sin(theta/2.)**2
	anu = eic - ep
	q2 = 4.*eic*ep*s2
	
!	write (6,*) 'in bdis, got q2 of ', q2
	rq2=q2		! to pass to w1w2
!	write(6,*)eic,ep,thec,esep,q2, rq2,

!
	amp = 0.938273
	x = q2/2./amp/(eic-ep)
	amp2 = amp*amp
	amn = 0.939565 
	if(aa.eq.2) then
	  ama_1=0.93827
	elseif(aa .gt. 2 .and. aa .lt. 4)then !he3
	  ama_1 = 1.87609
	elseif (aa.eq.4) then
	  ama_1 = 2.8094
	elseif (aa.eq.9) then
	  ama_1 = 7.5402
	elseif (aa.eq.12) then
	  ama_1 = 10.25553
	elseif (aa.eq.27) then
	  ama_1 = 24.205
	elseif (aa.eq.56) then
	  ama_1 = 51.1743
c	elseif (aa.eq.63) then
	elseif (aa.eq.64) then
	  ama_1 = 57.6886
	elseif (aa.eq.197) then
	  ama_1 = 182.5394
	endif
	ama = ama_1 + amp - esep
!	write (6,*) 'ama_1 is ', ama_1
!
	q3 = sqrt(q2 +anu*anu)
	w1ap = 0.0
	w2ap = 0.0
	w1an = 0.0
	w2an = 0.0
	du = 2./innt
	akf = (1. - exp(-aa/8.))*.22 + 0.04
!
	dp = pmax/innp
!
!	calculate proton and neutron free structure functions as a check if
!	desired
!	inp = 1
!	w0 = sqrt(amp2 + 2.*amp*amu - q2)
!	call w1w2(q2,w0,w1,w2,inp)
!	w2p = w2
!	inp = 2
!	call w1w2(q2,w0,w1,w2,inp)
!	w2n = w2
	wp_max = 0.0
!	do smearing correction
	do ip = 1,innp
	  p = (ip - .5)*dp
	  w1ap1 = 0.0
	  w2ap1 = 0.0
	  w1an1 = 0.0
	  w2an1 = 0.0
	  
	  do 10 iu = 1,innt

	    u = (iu - .5)*du
	    u = u - 1.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	i removed this dependence on kf 12/4/87 just to see the
!	effect
!       if(p.le.akf)then
!       
	    ei = ama - sqrt(p*p+ama_1*ama_1)
!       
!       else
!       ei = amd - sqrt(p*p + amp**2)
!       endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    anup = (ei*anu - p*u*q3)/amp
!	note that w1 and w2 are zero if wp < M + mpi
!	
	    wpmin2 = (amp + .135)**2
cdg	    if (iu.eq.1) then
!       write (6,*)'This is ', this
!	      write (6,*) 'wpmin2', wpmin2, ' thingie ', (ei*ei - p*p +
!	1	2.*ei*anu - 2.*p*u*q3 - q2), ei, p, anu, u, q3, q2
!       write (6,*) 'Ei is,',ei
cdg	    endif
	    if((ei*ei - p*p + 2.*ei*anu - 2.*p*u*q3 - q2).le
     >       .wpmin2) goto 10
	    radical =(ei*ei - p*p + 2.*ei*anu - 2.*p*u*q3 - q2)
	    wp = sqrt(radical)
				!write (6,*) 'wp (nu) is ',wp
				!write (6,*) 'real nu is ', anu
	    inp = 1
				!write(6,*)'In if stmt, about to go to w1w2'
	    w11=0.
	    w22=0.
	    nu = (wp**2 - amp**2 + q2) / 2. / amp
!	    write(*,*) 'my nu is ', nu, wp**2, q2, rq2
cdg	    call w1w2(rq2,wp,w1,w2,inp)
c	    call F1F2IN06(1.,1.,rq2,wp**2,w11,w22)
	    call F1F2IN07(1.,1.,rq2,wp**2,w11,w22)
	    if((w11.eq.0).or.(w22.eq.0)) goto 10
	    w11=w11/amp
	    w22=w22/nu

	    w1=w11
	    w2=w22
!	    if ((ip.eq.1).and.(iu.eq.1)) then
!	      write(12,*) 'D P 1', w1, w2, w11, w22, wp**2, nu
!	    endif 
!	    if ((ip.eq.1).and.(iu.eq.1)) then
!	      write(12,*) 'D P 1', sngl(w1), sngl(w2), sngl(w11), sngl(w22), sngl(wp**2), sngl(nu), sngl(p), sngl(ei)
!	    endif 
	    if (wp .gt. wp_max) wp_max = wp

	    arg1 = w1 + (1 -u*u)*p*p*w2/2./amp2
	    arg21 = (1. + p*u*q2/amp/anup/q3)**2*anup*anup/anu/anu
!       errorinBR	arg21 = (1. - p*u*q2/amp/anup/q3)*anup*anup/anu/anu
	    arg22 = (1.- u*u)*p*p*q2/2./amp2/q3/q3
	    arg2 = w2*(arg21 + arg22)
	    w1ap1 = w1ap1 + arg1*du
	    w2ap1 = w2ap1 + arg2*du
	    inp = 2
cdg	    call w1w2(rq2,wp,w1,w2,inp)	    
	    w11=0.
	    w22=0.
c	    call F1F2IN06(0.,1.,rq2,wp**2,w11,w22)
	    call F1F2IN07(0.,1.,rq2,wp**2,w11,w22)
	    if((w11.eq.0).or.(w22.eq.0)) goto 10
	    w11=w11/amp
	    w22=w22/nu

	    w1=w11
	    w2=w22

!	    if ((ip.eq.1).and.(iu.eq.1)) then
!	      write(12,*) 'D P 0', sngl(w1), sngl(w2), sngl(w11), sngl(w22), sngl(wp**2), sngl(nu), sngl(p), sngl(ei)
!	    endif 
	    arg1 = w1 + (1 -u*u)*p*p*w2/2./amp2
!       errorinBR	arg21 = (1. - p*u*q2/amp/anup/q3)*anup*anup/anu/anu
	    arg21 = (1. + p*u*q2/amp/anup/q3)**2*anup*anup/anu/anu
	    arg22 = (1.- u*u)*p*p*q2/2./amp2/q3/q3
	    arg2 = w2*(arg21 + arg22)
	    w1an1 = w1an1 + arg1*du
	    w2an1 = w2an1 + arg2*du
!       
 10	  continue

!	  weight1 = (ag*bg/pi*exp(-bg*p**2) + cg*dg/pi*exp(-dg*p**2) )*p*p*dp
!	  p=p*1000.
!	  dp=dp*1000.

	  weight1 = ((exp(-(ag*p)**2)*(f0-bigB)*alpha**2)*(((ag**2)
     >     *(alpha**2+p**2)+1)/(alpha**2+p**2)**2))
c	  write(*,*) 'here ', weight1,ag,f0,bigB,alpha,p

	  
	  if (aa.eq.2) then
	    if(p.lt.0) then
	      weight1=weight1-(bg*bigB*exp(-bg*abs(p))/(2*p))
	    else
	      weight1=weight1+(bg*bigB*exp(-bg*abs(p))/(2*p))
	    endif
	  else
	      weight1=weight1+(bg**2*bigB*exp(-(bg*p)**2))
	  endif

	  
!	  p=p/1000
	  weight1=weight1*p*p*dp/pi!*2.48176
c	  write(*,*) 'here 2', weight1
	  
!	   deut factors  
!	  agd =   0.75287E+01
!	  bgd =   0.18856E+03
!	  cgd =   0.22657E+00
!	  dgd =   0.14649E+02

!carbon factors	  
!	  agd =   0.28757E+01
!	  bgd =   0.41922E+02
!	  cgd =   0.33801E+00
!	  dgd =   0.13824E+02
	  weight2 = (agd*bgd/pi*exp(-bgd*p**2) + cgd*dgd/pi*exp(-dgd*p
     >     **2) )*p*p*dp
!	  p=p/1000
!	  dp=dp/1000

c	  write(12,*) 'weight ', p, weight1, weight2
!	  weight1=weight2
!	  if (weight1.lt.0) write(*,*) 'uhoh', p, weight1, ag, bg,
!	1   bigB, g0, alpha
	  
	  w1ap = w1ap + w1ap1*weight1
	  w2ap = w2ap + w2ap1*weight1
	  w1an = w1an + w1an1*weight1
	  w2an = w2an + w2an1*weight1
	enddo
	
	fackf = 2.*pi
!	
	w1a = fackf*(zz*w1ap + ann*w1an)
	w2a = fackf*(zz*w2ap + ann*w2an)
!       
!	write(6,*)  fackf,zz,ann,w1ap,w1an,w2ap,w2an
	sigm = cos(theta/2.)/(2.*137.*eic*s2)
!        sigm=0.389379*sigm**2  !(GeV**2 mb )
	sigm = (0.1973*sigm)**2
	tt2 = tan(theta/2.)**2
	sigdeep = 1.e7*sigm*(w2a + 2.*w1a*tt2)
!	write(6,*)'End of Bdis, DIS, MOTT',sigdeep, sigm, tt2, w1a, w2a
!	write(50,*)x,q2,wp_max
	return
	end
	
