	subroutine sigmodel_calc(e1,e2,th,z,a,m_tgt,sig_dis,sig_qe,sig)
C       +______________________________________________________________________________
c	
c       Subroutine to compute model unradiated cross section from a combination
c       of Y-scaling and deep inelastic models. The DEEPSIG model can be smeared
c       with a gaussian in missing mass squared of width RESOL$ (one sigma).
c       If SMEAR$ is .TRUE., such smearing will occur. If SMEAR$ is .FALSE.,
c       smearing is suppresed. It should be noted that the smearing involves
c       computing a convolution integral, which significantly slows down
c       this subroutine!
c       
c       ARGUMENTS:
c       
c       E1:		-	Incident energy in GeV.
c       E2:		- Scattered energy in GeV.
c       TH:		- Scattering angle in Degrees.
c       A:		- 'A' of nucleus.
c       Z:		- Number of protons in nucleus.
c       M_TGT:	- Mass of target nucleus in GeV/c^2.
c       M_REC:	- Mass of recoiling nucleon in GeV/c^2.
c       E_SEP:	- Separation energy for target nucleus in GeV/c^2.
c       SIG  :	- Calculated cross section in nb/(MeV-ster).
C       ______________________________________________________________________________

        implicit none
	include 'math_physics.inc'
	include 'include.inc'

C       Declare arguments.

	real*8		e1,e2,th,a,z,m_tgt,aux(7)
	real*4          e1pass,e2pass,thpass,mpass
	real*8          sig,factpass,sig_dis_pass,sig_qe_pass
	integer         zpass,apass
	logical		modelhack

C       Declare locals.

	real*8	 	sig_qe,sig_dis,y,normfac,fact
        real*8		thr,cs,sn,tn,elastic_peak
	real*8          Q2,nu,WSQ, x, pmax,nn
        real*8          ld2_a,ld2_z,ld2_nn,ld2_aux(7),ld2_sig_dis
        real*8          emc_corr,ld2_inta
        real*8          x1,x2,sig_dis_emc,sig_donal,sig_dis_donal
	real*8          F1,F2,W1,W2,Rc,sigmott
        real*8          frac,corfac
	real*8          dum0,dum1,dum2,dum3,dum4,dum5,m_atom
	real*8          ff1,deltae,e1cc,e2cc
	integer         inta, intz
	integer         xflag !flag for which xsec to calculate 1=both 2=QE only 3=DIS only
	logical         first


	save

        real*8 emc_func_xem
	external emc_func_xem

	real*8 emc_func_slac
	external emc_func_slac

	data first/.true./

	xflag=1


C TEMP TEST TEMP TEST !!!!!!!
C Here - I'll try to make Coulomb corrections part of the Born Model
C First - redefine energies

c	deltae = 0.0102
	deltae=0.0
	
	ff1 = (e1+deltae)/e1

	e1cc = e1 + deltae
	e2cc = e2 + deltae


c	write(6,*) 'in sigmodel_calc',e1,e2,th,z,a,m_tgt,xflag
	if(first) then
	   first=.false.
C  Here, all we're concerned about are the aux parameters
	   call target_info(0.0,dum0,dum1,dum2,z,a,dum3,dum4,aux,dum5)


	   do nn=1,197
	      fermip(nn)=.2
	      EPSN(nn)=0.01
	   enddo

	   fermip(197)=.264
	   fermip(64)=.260
	   fermip(27)=.240
	   fermip(12)=.221
	   fermip(4)=.160
	   fermip(3)=.160D0
	   fermip(2)=.160D0

	   epsn(197)=.006
	   epsn(64)=0.01
	   epsn(27)=.0083
	   epsn(12)=.016
	   epsn(4)=.02
	   epsn(3)=.0055
	   epsn(2)=.0022




	   aa=1.
	   bb=0.
	   cc=0.
	   dd=0.
	   ee=0.
	   ff=0.
	   if (a.eq.1) then
	      m_atom=1.00794
	   elseif(a.eq.2) then
	      m_atom=2.0141
c	      aa=0.093235
c	      bb=3.24433
c	      cc=-0.598669
c	      dd=-0.83881
	      aa = 1.72816025139459
	      bb = 2.53114253906281
	      cc = -2.72505067059703
	      dd = -1.58637090989274
	      ee = -16.3972900673533
	   elseif(a.eq.3) then
	      m_atom=3.0160
c	      aa = 1.36909582093995
c	      bb = -0.397730607400245
c	      cc = -0.00610340423044023
c	      dd = 0.07861291737992
	      bb = 0.8
	      cc = 0.06864880408328
	      dd = -0.320972192919132
	      ee = 0
	      aa = 0.552199789237622
	   elseif(a.eq.4) then
	      m_atom=4.002602
c	      aa = 1.08000840429859
c	      bb = -0.0982300330454544
c	      cc = 0.000586613632274299
c	      dd = 0.0316207297635641
	      bb = 0.466102123780537
	      cc = 0.0156369553828335
	      dd = -0.122243059123825
	      aa = 0.682462282515971
	   elseif(a.eq.9) then
	      m_atom=9.012182
c	      aa = 1.025
c	      bb = -0.0287152887511119
c	      cc = 0.00445485764339834
c	      dd = 0.0		!032096131872732
	      bb = 0.463011918692135
	      cc = 0.0125252624899601
	      dd = -0.101843839213646
	      aa = 0.674455752091906
	   else if (a.eq.12) then
	      m_atom=12.0110
c	      aa = 1.27376778137511
c	      bb = -0.331919290386685
c	      cc = -0.000807971358556216
c	      dd = 0.0698843279899779
	      bb = 0.222436834975864
	      cc = 0.00769270345172033
	      dd = -0.060282702596254
	      aa = 0.840262866196151
	   else if (a.eq.27) then
	      m_atom=26.98
	   else if (a.eq.56) then
	      m_atom=55.8470
	   elseif(a.eq.64) then
	      m_atom=63.546
c	      aa= 0.668888        
c	      bb = 0.336285        
c	      cc = 0.00192735     
c	      dd= -0.00481499   

	      bb = 0.041323394008416
	      cc = 0.00447016532137148
	      dd = -0.0303635977582275
	      aa = 1.00675406673173
	   else if (a.eq.197) then
	      m_atom=196.9665
c	      aa = 0.735285        
c	      bb = 0.309688       
c	      cc  = -0.000453334   
c	      dd  = 0.00294095      
	      bb = 0.0667337559531751
	      cc = 0.00448892579200859
	      dd = -0.0334460588480325
	      aa = 0.981274819686673
	   else
	      write(6,*) 'CANT FIGURE OUT M_ATOM.' 
              write(6,*) 'GET RID OF THIS CRAP!!'
	   endif
	   
	endif

	sig =0.0
C       If at or beyond elastic peak, return ZERO!

	thr = th*d_r
	cs = cos(thr/2.)
	sn = sin(thr/2.)
	tn = tan(thr/2.)
	elastic_peak = e1cc/(1.+2.*e1cc*sn**2/m_tgt)
c       dg	if (e2.ge.elastic_peak) then
c       dg	   sig = 0.D00
c       dg	   return
c       dg	endif
c	write(6,*) 'got to 1'

	Q2 = 4.*e1cc*e2cc*sn**2
	nu=e1cc-e2cc
	WSQ = -Q2 + m_p**2 + 2.0*m_p*nu 
        x = Q2/2/m_p/nu


c	innt = 30.0
c	innp = 30.0

	nn=a-z
	inta=a
	intz=z

	pmax=1.0

	
	


	sig_qe=0.0
	fact = 0.0
	if((xflag.eq.1).or.(xflag.eq.2)) then
	   if(a.gt.1.5) then
	      call f_to_sig(e1cc,e2cc,th,a,z,m_tgt,aux,sig_qe,y,fact)
	   else
	      sig_qe=0.0
	   endif
	   sig_qe = sig_qe*ff1**2
	endif
	factpass = fact
c	write(6,*) 'cheesy poofs', sig_qe,fact

	sig_dis=0.0
	if((xflag.eq.1).or.(xflag.eq.3)) then
c----------------------------------------------------------------
c       
c       do inelastic stuff only for ld2, for nuc targets sig_dis=ld2*emc
c       for now set all the values  to ld2 
c       aux is the ld2 qe param  from target_info.f file


cccccccccccccccc
ccc     dont forget to change this  stuff when ld2 qe param changes
cccccccccccccccc
	
	   ld2_a=2
	   ld2_z=1
	   ld2_nn=1
	   ld2_aux(3)= 0.0091109
	   ld2_aux(4)= 0.00151273
	   ld2_aux(5)= 0.00937085
	   ld2_aux(6)= 0.0108193
	   ld2_aux(7)= 46.4
	   ld2_inta=ld2_a


c--------------------------------------------------------------------
	   x1=0.8		!x<x1 --> use emc corrected ld2
	   x2=0.9		!x1<=x<x2 --> smooth transition from emc corrected ld2 to donals smearing
				!x>=x2  --> donal's smearing

	   if(a.ge.2) then

	      if(WSQ.lt.2.25) then
		 innt=30
		 innp=30
	      else
		 innt=30
		 innp=10
	      endif
	      
	      call bdisnew4he3(e1cc,e2cc,th,a, z,nn, epsn(inta), 
     >  	   pmax, innp, innt,  aux(3),aux(4), aux(5), 
     >             aux(6),aux(7),sig_donal)

c	      write(6,*) 'back from bdisnew',sig_donal,e1cc,e2cc,th

	      if (x.lt.x1) then
		 emc_corr = emc_func_xem(x,a)
	      elseif ((x.ge.x1).and.(x.lt.x2)) then 
		 frac = (x-x1)/(x2-x1)
		 emc_corr = 1.0*frac + emc_func_xem(x,a)*(1.-frac) 
	      elseif(x.ge.x2) then
		 emc_corr = 1.0
	      endif

	      sig_dis_donal=sig_donal*emc_corr
	      sig_dis = sig_dis_donal

	   else  ! hydrogen
c	      call F1F2IN06(Z, A, Q2, WSQ, F1, F2)
	      call F1F2IN07(Z, A, Q2, WSQ, F1, F2)
C       Convert F1,F2 to W1,W2
	      W1 = F1/m_p
	      W2 = F2/nu


C       Mott cross section
	      sigmott=(19732.0/(2.0*137.0388*e1cc*sn**2))**2*cs**2/1.d6
	      sig_dis = 1d3*sigmott*(W2+2.0*W1*tn**2)
           endif

CDJGC Simple model
CDJG	   if(wsq.gt.1.151915) then
CDJG	      call F1F2IN06(1, 1, Q2, WSQ, F1, F2)
CDJG	   else
CDJG	      F1=0.0
CDJG	      F2=0.0
CDJG	   endif
CDJGC       Convert F1,F2 to W1,W2
CDJG	   W1 = F1/m_p
CDJG	   W2 = F2/nu
CDJG	   if(wsq.gt.1.151915) then
CDJG	      call F1F2IN06(0, 1, Q2, WSQ, F1, F2)
CDJG	   else
CDJG	      F1=0.0
CDJG	      F2=0.0
CDJG	   endif
CDJG	   W1 = W1 + F1/m_p
CDJG	   W2 = W2 + F2/nu
CDJG
CDJG	   if(x.lt.1.0) then
CDJG	      emc_corr= emc_func_xem(x,a)
CDJG	   else
CDJG	      emc_corr = emc_func_xem(1.0,a)
CDJG	   endif
CDJG
CDJG
CDJG	   W1=W1*emc_corr*a/2.0
CDJG	   W2=W2*emc_corr*a/2.0
CDJG
CDJGC       Mott cross section
CDJG	   sigmott=(19732.0/(2.0*137.0388*e1*sn**2))**2*cs**2/1.d6
CDJG	   sig_dis = 1.d3*sigmott*(W2+2.0*W1*tn**2)
CDJGC END simple model

c    do a high x tweak for the inelastic part of nuc targets
	   if ((x.gt.0.9).and.(a.gt.1.5)) then	   
cdgc	      call  highx_cor(a,x,corfac)
	      call  dis_highx_cor(a,x,corfac)
	      sig_dis = sig_dis*corfac	   
	   endif


	   sig_dis=sig_dis/1.0d3
	   if (sig_dis.lt.0.) sig_dis=0.0
	   
	   sig_dis = sig_dis*ff1**2

	endif ! test on xsec type
	sig = sig_qe + sig_dis

	sig = sig*1000.0 !to be consistent with Peter's units

c	sig=sig_qe

	sig_qe_pass = sig_qe
	sig_dis_pass = sig_dis

c	write(6,*) 'leaving sigmodel_calc',sig_qe,sig_dis,sig

	return
	end





c-------------------------------------------------------------------------------------------
	real*8 function emc_func_slac(x,A)
	real*8 x,A,atemp
	real*8 alpha,C

	atemp = A
				!	if(A.eq.4.or.A.eq.3) then  ! emc effect is more like C for these 2...
				!	   atemp = 12
				!	endif

	alpha = -0.070+2.189*x - 24.667*x**2 + 145.291*x**3
     >       -497.237*x**4 + 1013.129*x**5 - 1208.393*x**6
     >       +775.767*x**7 - 205.872*x**8

	C = exp( 0.017 + 0.018*log(x) + 0.005*log(x)**2)
	
	emc_func_slac = C*atemp**alpha
	return 
	end


c aji note:
c polynomial fit made to inelastic emc ratios from xem data
c at low x (x<0.3) the fit is constrained with world data (to get some sort of shadowing behaviour).

C DJG NOTE 7/11/2007
c Note that although this is called "emc_func", it is really a fit to the ratio of the emc effect to 
c a pure smearing claculation. So - this function is applied to the smeared n+p cross section
c to reproduce the correct emc effect. If you plot these functions, it will not look like
c the emc effect, so don't panic.
  
	real*8 function emc_func_xem(x,A) ! now compute the emc effect from our own fits.
	implicit none
        real*8 x,a,xtmp
	real*8 emc
	
	
c	if (x.le.1.0) then
	if(x.le.0.9) then
	   xtmp = x
	else
	   xtmp = 0.9
	endif

	emc =1.0
CDeuterium******************************************************************************
	if(A.eq.2) then
C it 2
c	   if(xtmp.lt.0.2) then
c	      emc=1.06
c	   else
c	      emc = 0.79515 +1.9777*xtmp - 3.9724*xtmp**2 -0.66967*xtmp**3
c	1	   +8.3082*xtmp**4 - 5.5263*xtmp**5
c	   endif
c	   emc = emc*0.96689
c it 3
c	   emc = 0.70443 +2.3742*xtmp - 4.6566*xtmp**2 -0.78540*xtmp**3
c     >   	+9.3838*xtmp**4 - 6.1256*xtmp**5
c it4
	   emc = 0.83818 +1.4201*xtmp - 2.8189*xtmp**2 -0.67289*xtmp**3
     >      +6.0543*xtmp**4 - 3.8575*xtmp**5
CHe3***********************************************************************************
	else if(A.eq.3) then
C it 2
c	      emc = 1.0118 +1.1029*xtmp -2.5081*xtmp**2 - 0.22033*xtmp**3
c	1	   + 4.8120*xtmp**4 - 3.2865*xtmp**5
C it 3
c	   emc = 0.92170 +1.7544*xtmp -3.7324*xtmp**2 - 0.24293*xtmp**3
c     >  	+ 6.7613*xtmp**4 - 4.6089*xtmp**5
c it4
	   emc = 0.95649 +1.2695*xtmp -2.3735*xtmp**2 - 0.62320*xtmp**3
     >        + 4.4999*xtmp**4 - 2.7871*xtmp**5
CHe4***********************************************************************************	      
	else if(A.eq.4) then
C it2
c            emc = 0.84622 + 2.2462*xtmp - 4.7909*xtmp**2
c	1	   + 0.065713*xtmp**3 + 7.6154*xtmp**4 - 5.2029*xtmp**5
C it3
c	   emc = 0.70050 + 3.1241*xtmp - 6.1738*xtmp**2
c     >  	- 0.049988*xtmp**3 + 9.3053*xtmp**4 - 6.1348*xtmp**5
C it4
	   emc = 0.85453 + 2.0067*xtmp - 4.0285*xtmp**2
     >       + 0.20116*xtmp**3 + 5.3353*xtmp**4 - 3.5246*xtmp**5
C Be**********************************************************************************
	else if(A.eq.9) then
C it 2
c	      emc = 0.80887 + 3.9354*xtmp - 8.6056*xtmp**2 -0.16342*xtmp**3
c	1	   + 14.074*xtmp**4 -9.3065*xtmp**5
C it 3
c	   emc = 0.46324 + 6.1220*xtmp - 12.184*xtmp**2 -1.0956*xtmp**3
c     >  	+ 20.316*xtmp**4 -12.899*xtmp**5
C it 4
	emc = 0.67858 + 4.5917*xtmp - 9.0413*xtmp**2 -1.4294*xtmp**3
     > 	   + 15.497*xtmp**4 -9.4333*xtmp**5
C Carbon**********************************************************************************
	else if(A.eq.12) then
C it 2
c         emc = 0.8279 + 3.5070*xtmp -7.5807*xtmp**2 
c	1	   -0.60935*xtmp**3 +13.081*xtmp**4 -8.5083*xtmp**5
C it 3
c	   emc = 0.63653 + 4.6458*xtmp -9.2994*xtmp**2 
c     >  	-1.2226*xtmp**3 +16.157*xtmp**4 -10.236*xtmp**5
C it 4
	   emc = 0.72341 + 3.9547*xtmp -7.7589*xtmp**2 
     >      -1.5734*xtmp**3 +14.090*xtmp**4 -8.6972*xtmp**5
C Al**********************************************************************************
	else if(a.eq.27) then
	   emc = 0.98645 + 3.0385*xtmp - 22.072*xtmp**2 + 74.981*xtmp**3
     >  	- 132.97*xtmp**4 + 113.06*xtmp**5 -35.612*xtmp**6
C Copper**********************************************************************************
	else if(A.eq.64) then 
C it 2
c	      emc = 1.1075 + 2.7709*xtmp - 6.5395*xtmp**2 -0.46848 *xtmp**3
c	1	   +10.534*xtmp**4 - 6.6257*xtmp**5
c it 3
c	   emc = 0.58372 + 6.0358*xtmp - 11.988*xtmp**2 -1.0211*xtmp**3
c     >  	+18.567*xtmp**4 - 11.482*xtmp**5
c it 4
	   emc = 0.70732 + 5.0635*xtmp - 9.5274*xtmp**2 -1.7399*xtmp**3
     >     +14.835*xtmp**4 - 8.4185*xtmp**5
C Gold**********************************************************************************	      
	else if(A.eq.197) then
C it 2
c	      emc = 1.1404 + 4.0660*xtmp -10.318*xtmp**2 -1.9036*xtmp**3
c	1	   + 21.969*xtmp**4 - 14.461*xtmp**5	   
C it 3
c	   emc = 0.44132 + 8.1232*xtmp -16.141*xtmp**2 -5.6562*xtmp**3
c     >  	+ 35.606*xtmp**4 - 22.008*xtmp**5
C it4
	      emc = 0.79377 + 5.8630*xtmp -11.149*xtmp**2 -5.8010*xtmp**3
     >	   + 26.139*xtmp**4 - 14.992*xtmp**5	
	      
	else  
	   write(*,*) '** in emc_func_xem, unknown target'
	   stop		
	endif
	
	emc_func_xem= emc
	return
	end

c-----------------------------------------------------------------------------------------------
c 

	subroutine highx_cor(anuc,x,cor)

	real*8 x,cor,anuc

         if(anuc.eq.3) then
	   if ((x.gt.0.9).and.(x.lt.1.4)) then
	    cor= -0.908273 + (4.13702*x) -(2.11462*x**2)
	   elseif (x.ge.1.4) then
	      cor=0.74
	   endif
        
        elseif(anuc.eq.4) then
	   if ((x.gt.0.9).and.(x.lt.1.17)) then
	      cor= 3.24553 - (3.47244*x) +  (1.11309*x**2)
	   elseif (x.ge.1.17) then
	      cor=0.7
	   endif
	elseif(anuc.eq.9) then

	   if ((x.gt.0.9).and.(x.lt.1.26)) then
	      cor= 0.779378 + (1.84808*x) - (1.7588*x**2)
	   elseif (x.ge.1.26) then
	      cor=0.3
	   endif

	elseif(anuc.eq.12) then
	   if ((x.gt.0.9).and.(x.lt.1.26)) then
	      cor=  1.09301 + (0.798708*x) - (0.939027*x**2)
	   elseif (x.ge.1.26) then
	      cor=0.55
	   endif

        elseif(anuc.eq.64) then
	   if ((x.gt.0.9).and.(x.lt.1.3)) then
	      cor=  3.3275 - (3.94771*x) + (1.496*x**2)
	   elseif (x.ge.1.3) then
	      cor=0.68
	   endif

          elseif(anuc.eq.197) then
	   if ((x.gt.0.9).and.(x.lt.1.3)) then
           cor= -0.389135+ (3.08065*x)- (1.66297*x**2)
	   elseif (x.ge.1.3) then
	      cor=0.8
	   endif

	else
	   cor=1.
	endif

	return
	end
c-----------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------

	subroutine dis_highx_cor_old(anuc,x,cor)
        implicit none
	real*8 x,cor,anuc, frac,xlow1,xhigh1,xlow2,xhigh2

	xlow1=0.9
	xhigh1=0.95

	xlow2=1.3
	xhigh2=1.4
	
	frac=1.
	if(anuc.eq.3) then
	   cor=-2.12112*x+3.03449
	   if((x.ge.xlow1).and.(x.le.xhigh1)) then
	     frac = (x-xlow1)/(xhigh1-xlow1)
	   endif
	   cor=frac*cor+1.-frac
	   
        
	 elseif(anuc.eq.4) then
	   cor=-1.76466*x+2.68897
	   if((x.ge.xlow1).and.(x.le.xhigh1)) then
	     frac = (x-xlow1)/(xhigh1-xlow1)
	   endif
	   cor=frac*cor+1.-frac
	   
	 elseif(anuc.eq.9) then
	   cor=-1.8383*x+2.77253
	   if((x.ge.xlow1).and.(x.le.xhigh1)) then
	     frac = (x-xlow1)/(xhigh1-xlow1)
	   endif
	   cor=frac*cor+1.-frac
	   
	 elseif(anuc.eq.12) then
	   cor=-1.32193*x+2.28754
	   if((x.ge.xlow1).and.(x.le.xhigh1)) then
	     frac = (x-xlow1)/(xhigh1-xlow1)
	   endif
!	   if((x.ge.xlow2).and.(x.le.xhigh2)) then
!	     frac = (x-xlow2)/(xhigh2-xlow2)
!	     frac=1.-frac
!	   endif
!	   if(x.gt.xhigh2) frac=0.
	   cor=frac*cor+1.-frac
!	   write(21,*) 'cor is ', cor, x
	 elseif(anuc.eq.64) then
!	   cor=-2.21331*x+3.02106
	   cor=-1.46912*x+2.31581
	   if((x.ge.xlow1).and.(x.le.xhigh1)) then
	     frac = (x-xlow1)/(xhigh1-xlow1)
	   endif
	   cor=frac*cor+1.-frac
	 elseif(anuc.eq.197) then
 	   cor= -1.72192*x+2.65671
	   if((x.ge.xlow1).and.(x.le.xhigh1)) then
	     frac = (x-xlow1)/(xhigh1-xlow1)
	   endif
	   cor=frac*cor+1.-frac
	   
	 else
	   cor=1.
	 endif
	 
	 if(cor.lt.0.4) cor=0.4

	return
	end
c--------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------

	subroutine dis_highx_cor(anuc,x,cor)
        implicit none
	real*8 x,cor,anuc, frac,xlow1,xhigh1,xlow2,xhigh2

	xlow1=0.9
	xhigh1=0.95

	xlow2=1.9
	xhigh2=2.0
	
	frac=1.
	cor=1.


	if(anuc.eq.2) then
	   cor=1.
	   cor=-3.30482*x+ 4.10442
	   if((x.ge.xlow1).and.(x.le.xhigh1)) then
	      frac = (x-xlow1)/(xhigh1-xlow1)
	   endif
	   cor=frac*cor+1.-frac
	elseif(anuc.eq.3) then
	   xlow1=1.
	   xhigh1=1.15
	   if (x.gt.xlow1) then
	      if (x.lt.1.15) then
		 cor=-4.80303*x+ 5.74758
	      else
		 cor=0.5
	      endif
	      if((x.ge.xlow1).and.(x.le.xhigh1)) then
		 frac = (x-xlow1)/(xhigh1-xlow1)
	      endif
	   endif
	   cor=frac*cor+1.-frac
	 elseif(anuc.eq.4) then
           cor=-1.78944*x+ 2.63272
	   if((x.ge.xlow1).and.(x.le.xhigh1)) then
	     frac = (x-xlow1)/(xhigh1-xlow1)
	   endif
	   cor=frac*cor+1.-frac
	 elseif(anuc.eq.9) then
	   cor=-1.7549631060248*x+ 2.6646629067298
	   if((x.ge.xlow1).and.(x.le.xhigh1)) then
	     frac = (x-xlow1)/(xhigh1-xlow1)
	   endif
	   cor=frac*cor+1.-frac
	   
	 elseif(anuc.eq.12) then
	   cor=-1.29213*x+ 2.2087
	   if((x.ge.xlow1).and.(x.le.xhigh1)) then
	     frac = (x-xlow1)/(xhigh1-xlow1)
	   endif
!	   if((x.ge.xlow2).and.(x.le.xhigh2)) then
!	     frac = (x-xlow2)/(xhigh2-xlow2)
!	     frac=1.-frac
!	   endif
!	   if(x.gt.xhigh2) frac=0.
	   cor=frac*cor+1.-frac
!	   write(21,*) 'cor is ', cor, x
	 elseif(anuc.eq.64) then
	   cor=-1.65829142599487*x+ 2.48174872208596
	   if((x.ge.xlow1).and.(x.le.xhigh1)) then
	     frac = (x-xlow1)/(xhigh1-xlow1)
	   endif
	   cor=frac*cor+1.-frac
!	   cor=1.0
	 elseif(anuc.eq.197) then
	   cor=-1.42430013496752*x+ 2.25789690593227

	   if((x.ge.xlow1).and.(x.le.xhigh1)) then
	     frac = (x-xlow1)/(xhigh1-xlow1)
	   endif
	   cor=frac*cor+1.-frac
	   
	 else
	   cor=1.
	 endif
	 
	 if(cor.lt.0.4) cor=0.4

	return
	end
c--------------------------------------------------------------------------------
