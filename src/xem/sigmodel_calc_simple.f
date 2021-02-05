	subroutine sigmodel_calc(e1,e2,th,z,a,mpass,sig_dis_pass,sig_qe_pass,sig)
C       +______________________________________________________________________________
c	
C       Calculate cross section using Peter's F1F209.f routine
c       
c       ARGUMENTS:
c       
c       E1:		-	Incident energy in GeV.
c       E2:		- Scattered energy in GeV.
c       TH:		- Scattering angle in Degrees.
c       A:		- 'A' of nucleus.
c       Z:		- Number of protons in nucleus.
c       M_TGT:	- Mass of target nucleus in GeV/c2.
c       M_REC:	- Mass of recoiling nucleon in GeV/c2.
c       E_SEP:	- Separation energy for target nucleus in GeV/c2.
c       SIG  :	- Calculated cross section in nb/(MeV-ster).
C       ______________________________________________________________________________

        implicit none
	include 'math_physics.inc'
	include 'include.inc'

C       Declare arguments.

	real*8		e1,e2,th,a,z,m_tgt
	real*8          e1pass,e2pass,thpass,mpass
	real*8          sig,factpass,sig_dis_pass,sig_qe_pass
	integer         zpass,apass
	logical		modelhack

C       Declare locals.

	real*8	 	sig_qe,sig_dis,y,normfac,fact
        real*8		thr,cs,sn,tn,elastic_peak
	real*8          Q2,nu,WSQ, x
	real*8          F1,F2,W1,W2,sigmott,r
	real*8          W1p,W1n,W1D,W2p,W2n,W2D
	real*8          inelastic_it
	real*8          f1mec,f2mec,W1mec,W2mec
	integer         xflag !flag for which xsec to calculate 1=both 2=QE only 3=DIS only
	logical         first


	save

        real*8 emc_func_xem
	external emc_func_xem

	real*8 emc_func_slac
	external emc_func_slac

	data first/.true./

	xflag=1
	m_tgt=mpass


	sig =0.0
	sig_qe=0.0
	sig_dis=0.0

C       If at or beyond elastic peak, return ZERO!

	thr = th*d_r
	cs = cos(thr/2.)
	sn = sin(thr/2.)
	tn = tan(thr/2.)
	elastic_peak = e1/(1.+2.*e1*sn**2/m_tgt)
      	if (e2.ge.elastic_peak) then
       	   sig = 0.0
       	   return
       	endif
c	write(6,*) 'got to 1'

	Q2 = 4.*e1*e2*sn**2
	nu=e1-e2
	WSQ = -Q2 + m_p**2 + 2.0*m_p*nu 
        x = Q2/2/m_p/nu



	F1=0
	F2=0
	r=0
	if((xflag.eq.1).or.(xflag.eq.3)) then
c----------------------------------------------------------------
c       
c       do inelastic stuff
c	   call F1F2IN09(Z, A, Q2, WSQ, F1, F2, r)
C Use old Bodek fit + SLAC EMC fit for now, b/c F1F2IN09 doesn't like large Q2,W2
	   if(wsq.gt.1.1664) then
	      call ineft(Q2,sqrt(wsq),W1p,W2p,1.0)
	      call ineft(Q2,sqrt(wsq),W1D,W2D,2.0)
	      W1n=2.0*W1D-W1p
	      W2n=2.0*W2D-W2p
C       Convert F1,F2 to W1,W2
c	   W1 = F1/m_p
c	   W2 = F2/nu
c	      call MEC2020(Z,A,wsq,Q2,f1mec)
c	      f2mec = 2.*x*f1mec/(1.+4.*x*x*m_p**2/q2)
c	      W1mec=f1mec/m_p
c	      W2mec=F2mec/nu

	      W1=Z*W1p+(A-Z)*W1n ! + W1mec
	      W2=Z*W2p+(A-Z)*W2n ! + W2mec
C       Mott cross section
	      sigmott=(19732.0/(2.0*137.0388*e1*sn**2))**2*cs**2/1.d6
	      if(A.gt.1.5 .and. A.lt.2.5) then
		 sig_dis = 1d3*sigmott*2.0*(W2D+2.0*W1D*tn**2) ! this is nb I think
	      else
		 sig_dis = 1d3*sigmott*(W2+2.0*W1*tn**2) ! this is nb I think
	      endif
	      if(A.gt.2.0) then
		 sig_dis = sig_dis*emc_func_slac(x, A)
	      endif
c	      write(6,*) 'cheesy poofs', e1,e2,th,q2,wsq,sig_dis
CDG apply "iteration" correction (used for XEM analysis)
CDG DO not use this for more "generic" stuff.
CDG	   sig_dis = sig_dis*inelastic_it(x,A)
	   else
	      W1=0.0
	      W2=0.0
	   endif
	endif


	if((xflag.eq.1).or.(xflag.eq.2)) then
	   call F1F2QE09(Z, A, Q2, WSQ, F1, F2)
C       Convert F1,F2 to W1,W2
	   W1 = F1/m_p
	   W2 = F2/nu
C       Mott cross section
	   sigmott=(19732.0/(2.0*137.0388*e1*sn**2))**2*cs**2/1.d6
	   sig_qe = 1d3*sigmott*(W2+2.0*W1*tn**2)
C Temp test - DJG May 23, 2013
c	   sig_qe=sig_qe/0.8
	endif


	sig = sig_qe + sig_dis !sig is already real*4
	

	sig_qe_pass = sig_qe ! pass back as real*4
	sig_dis_pass = sig_dis

	return
	end


	real*8 function inelastic_it(x,A)
C DJG: Correction to inelastic cross section to F1F209 from XEM
C DJG: 40 degree data. Just a simple one-pass iteration for use
C DJG: to check our model dependence.
	real*8 A,x
	real*8 p1,p2,p3,p4
	real*8 x1,x2,xit


	inelastic_it = 1.0 ! set to 1 by default

	if(A.gt.1.5.and.A.lt.2.5) then
	   x1=0.3172d0
	   x2=0.7275d0
	   if(x.lt.x1) then
	      xit=x1
	   elseif(x.gt.x2) then
	      xit=x2
	   elseif(x.ge.x1 .and. x.le.x2) then
	      xit=x
	   endif
	   inelastic_it=1.2394d0 - 2.1805d0*xit + 5.6853d0*xit**2
     >     -4.3908d0*xit**3
	endif

	if(A.gt.2.5 .and. A.lt.3.5) then
	   x1=0.3172
	   x2=1.055
	   if(x.lt.x1) then
	      xit=x1
	   elseif(x.gt.x2) then
	      xit=x2
	   elseif(x.ge.x1 .and. x.le.x2) then
	      xit=x
	   endif
	   inelastic_it=2.9235d0 - 16.075d0*xit + 46.426d0*xit**2
     >     -56.779d0*xit**3 + 25.007d0*xit**4
	endif

	if(A.gt.3.5 .and. A.lt.4.5) then
	   x1=0.3172
	   x2=1.0927
	   if(x.lt.x1) then
	      xit=x1
	   elseif(x.gt.x2) then
	      xit=x2
	   elseif(x.ge.x1 .and. x.le.x2) then
	      xit=x
	   endif
	   inelastic_it=1.505d0 - 4.8103d0*xit + 15.221d0*xit**2
     >     -19.713d0*xit**3 + 8.9693d0*xit**4
	endif

	if(A.gt.8.5 .and. A.lt.9.5) then
	   x1=0.6478
	   x2=1.0929
	   if(x.lt.x1) then
	      xit=x1
	   elseif(x.gt.x2) then
	      xit=x2
	   elseif(x.ge.x1 .and. x.le.x2) then
	      xit=x
	   endif
	   inelastic_it=0.89139d0 + 0.44009d0*xit -0.44163d0*xit**2
	endif

	if(A.gt.11.5 .and. A.lt.13.5) then
	   x1=0.3172
	   x2=1.005
	   if(x.lt.x1) then
	      xit=x1
	   elseif(x.gt.x2) then
	      xit=x2
	   elseif(x.ge.x1 .and. x.le.x2) then
	      xit=x
	   endif
	   inelastic_it=0.88964d0 + 0.39884d0*xit - 0.36051d0*xit**2
	endif

	if(A.gt.61.0 .and. A.lt.66.0) then
	   x1=0.3172
	   x2=1.005
	   if(x.lt.x1) then
	      xit=x1
	   elseif(x.gt.x2) then
	      xit=x2
	   elseif(x.ge.x1 .and. x.le.x2) then
	      xit=x
	   endif
	   inelastic_it=0.77646d0 + 0.90481d0*xit - 0.83815d0*xit**2
	endif


	if(A.gt.196.0 .and. A.lt.198.0) then
	   x1=0.3172
	   x2=1.005
	   if(x.lt.x1) then
	      xit=x1
	   elseif(x.gt.x2) then
	      xit=x2
	   elseif(x.ge.x1 .and. x.le.x2) then
	      xit=x
	   endif
	   inelastic_it=0.70439d0 + 1.0510d0*xit - 0.91679d0*xit**2
	endif

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
	1    -497.237*x**4 + 1013.129*x**5 - 1208.393*x**6
	2    +775.767*x**7 - 205.872*x**8

	C = exp( 0.017 + 0.018*log(x) + 0.005*log(x)**2)
	
	emc_func_slac = C*atemp**alpha
	return 
	end      
	   

	
      SUBROUTINE MEC2020(z,a,w2,q2,f1mec)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC   Subroutine for Transverse Enhancement new the QE and Delta due to meson   CCC
CCC   exchange currents and isobar excitations in the medium.  This is assumed  CCC
CCC   to be due to quasi-deuteron 2-body currents.  Shape is a distorted        CCC
CCC   Gaussian in W^2 with a cut-off near the 2-body threshold near x=2.        CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      
! fit to low q2 dip region: purefly empirical
! assume contribution is pure transverse
      implicit none
      real*8 z,a,q2,w2,mp/0.938272/,mp2,w,nu
      real*8 a1,b1,c1,t1,dw2
      real*8 x, f1mec
C use A=12 values
      real*8 xvalm(45) / 
     & 0.40388E+00,0.28053E+01,0.29253E+00,0.60248E+01,0.85436E+00,
     & 0.94254E-01,0.21100E+01,0.12762E+01,0.40249E+00,-.22073E+01,
     & 0.98110E+00,0.95991E+00,0.10305E+01,0.99548E+00,0.97845E+00,
     & 0.10000E+01,0.97772E+00,0.10062E+01,0.99096E+00,0.99684E+00,
     & 0.10051E+01,0.99674E+00,0.10032E+01,0.10070E+01,0.10054E+01,
     & 0.16075E+01,0.23972E+00,0.26369E+01,0.45830E+00,0.91716E+00,
     & 0.37683E-01,0.00000E+00,0.13844E+00,0.12306E+00,0.15551E+01,
     & 0.99849E+02,0.44784E-02,0.10040E-04,0.18568E-01,0.33686E-01,
     & 0.19603E+00,0.24499E+00,0.32459E+00,0.00000E+00,0.10000E+01 /

      mp2 = mp*mp
      f1mec = 0.0
      if(w2.le.0.0) return
      w  = sqrt(w2)
      nu = (w2 - mp2 + q2)/2./mp
      x  = q2/(2.0*mp*nu )

      if(A.lt.2.5) return

      a1 = A*q2**2*xvalm(1)*exp(-1.0*q2/xvalm(2))/
     &                                (xvalm(3)+q2)**xvalm(4)

      b1 = xvalm(5)+xvalm(6)*q2

      c1 = xvalm(33)+xvalm(34)*q2

      t1 = (w2-b1)**2/c1**2/2.

      dw2 = w2+q2*(1.-1./2.2)-1.0*mp2

      if(dw2.LT.0.0) dw2 = 0.0
      f1mec = a1*(exp(-1.*t1)*sqrt(dw2))

      
       if(f1mec.LE.1.0E-9 ) f1mec=0.0


      return
      end

!------------------------------------------------------------------------

       SUBROUTINE INEFT(QQ,W,W1,W2,amuM)                                      
                                                                        
C Modified 6feb87 by lww to accept target information passed through    
C common block /targt/.                                                 
                                                                        
C This program takes the old slac structure function model (Atwood,     
C Bodek, et.al.) and outputs values for W1 and W2 at given kinematics. 
! As of 11/3/95 this version is per NEUCLEON   ! Steve Rock

! amuM is atomic number, ie. 1. 2.xxx etc.

      Implicit None 
!      COMMON       /TARGT/ iZ, iA, avgN, avgA, avgM, amuM               
      REAL*8 QQ,W,W1,W2,amuM,WW,V,VV,OMEGAP,SP,UNIV,BRES,SLACF2,B
      REAL*8 VW2,X,EMCFAC
      REAL*8    C(24),CF(11),CD(24),CFD(11)                          
      REAL*8    EF(7) 
      REAL FITEMC
                                                  
      REAL*8         PM / .938256/,PMPM/.880324/,TPM/1.876512/            
      REAL*8         R /  .18/,ALPHAX/137.0388/,THCONST/0.0174533/
      LOGICAL GOODFIT
      common/testing/prttst
      logical prttst
      DATA         EF / -0.00136693,-.00510425,-.0375986,-.0946004,     
     +                  -.122435,-.0112751,0.406435/                    
                                                                        
C FINAL HYDROGEN COEFFS. FROM FITTING WITH NOTUSE=TRUE,ISLIDE=FALSE     
C CHISQ=3995,NO. OF POINTS=2533,FREE PARAMS=28,CHISQ/D.F.=1.59.         
C THIS DATA PROVIDED BY ARIE BODEK FOR E139 (10/3/83)                   
                                                                        
C C(24)=HYDROGEN COEFFICIENTS FOR B (BACKGROUND AND RESONANCE TERMS)    
                                                                        
      DATA   C(1) / 0.10741163E 01/,  C(2) / 0.75531124E 00/,           
     *       C(3) / 0.33506491E 01/,  C(4) / 0.17447015E 01/,           
     *       C(5) / 0.35102405E 01/,  C(6) / 0.10400040E 01/,           
     *       C(7) / 0.12299128E 01/,  C(8) / 0.10625394E 00/,           
     *       C(9) / 0.48132786E 00/,  C(10)/ 0.15101467E 01/,           
     *       C(11)/ 0.81661975E-01/,  C(12)/ 0.65587179E 00/,           
     *       C(13)/ 0.17176216E 01/,  C(14)/ 0.12551987E 00/,           
     *       C(15)/ 0.74733793E 00/,  C(16)/ 0.19538129E 01/,           
     *       C(17)/ 0.19891522E 00/,  C(18)/-0.17498537E 00/,           
     *       C(19)/ 0.96701919E-02/,  C(20)/-0.35256748E-01/,           
     *       C(21)/ 0.35185207E 01/,  C(22)/-0.59993696E 00/,           
     *       C(23)/ 0.47615828E 01/,  C(24)/ 0.41167589E 00/            
                                                                        
C CF(11)=HYDROGEN COEFFICIENTS FOR F2 (UNIVERSAL FUNCTION) OMEGAW FIT   
                                                                        
      DATA CF(1) / 0.25615498E 00/,  CF(2) / 0.21784826E 01/,           
     *     CF(3) / 0.89783738E 00/,  CF(4) /-0.67162450E 01/,           
     *     CF(5) / 0.37557472E 01/,  CF(6) / 0.16421119E 01/,           
     *     CF(7) / 0.37635747E 00/,  CF(8) / 0.93825625E 00/,           
     *     CF(9) / 0.10000000E 01/,  CF(10)/ 0.0           /,           
     *     CF(11)/ 0.50000000E 00/                                      
                                                                        
C FINAL DEUTERIUM COEFFS FROM FITTING WITH NOTUSE=TRUE,ISLIDE=FALSE     
C CHISQ=4456,NO. OF POINTS 2303,FREE PERAMS=26,CHISQ/D.F.=1.96.         
C THIS DATA PROVIDED BY ARIE BODEK FOR E139 (10/3/83)                   
                                                                        
C CD(24)=DEUTERIUM COEFFICIENTS FOR B (BACKGROUND AND RESONANT TERMS)   
                                                                        
      DATA  CD(1) / 0.10521935E 01/, CD(2) / 0.76111537E 00/,           
     *      CD(3) / 0.41469897E 01/, CD(4) / 0.14218146E 01/,           
     *      CD(5) / 0.37119053E 01/, CD(6) / 0.74847487E 00/,           
     *      CD(7) / 0.12399742E 01/, CD(8) / 0.12114898E 00/,           
     *      CD(9) / 0.11497852E-01/, CD(10)/ 0.14772317E 01/,           
     *      CD(11)/ 0.69579815E-02/, CD(12)/ 0.12662466E 00/,           
     *      CD(13)/ 0.15233427E 01/, CD(14)/ 0.84094736E-01/,           
     *      CD(15)/ 0.74733793E 00/, CD(16)/ 0.19538129E 01/,           
     *      CD(17)/ 0.19891522E 00/, CD(18)/-0.24480414E 00/,           
     *      CD(19)/ 0.14502846E-01/, CD(20)/-0.35256748E-01/,           
     *      CD(21)/ 0.35185207E 01/, CD(22)/-0.21261862E 00/,           
     *      CD(23)/ 0.69690531E 01/, CD(24)/ 0.40314293E 00/            
                                                                        
C CFD(11) ARE DEUTERIUM COEFFICIENTS FOR F2 (UNIVERSAL FUNCTION)        
C OMEGAW FIT                                                            
                                                                        
      DATA CFD(1) / 0.47708776E 00/, CFD(2) / 0.21601918E 01/,          
     *     CFD(3) / 0.36273894E 01/, CFD(4) /-0.10470367E 02/,          
     *     CFD(5) / 0.49271691E 01/, CFD(6) / 0.15120763E 01/,          
     *     CFD(7) / 0.35114723E 00/, CFD(8) / 0.93825625E 00/,          
     *     CFD(9) / 0.10000000E 01/, CFD(10)/ 0.0           /,          
     *     CFD(11)/ 0.50000000E 00/                                     
                                                                        
C COMPUTE SOME KINEMATIC QUANTITIES                                     
                                                                        
      WW     = W**2                                                     
      V      = (WW+QQ-PMPM)/2.D0/PM                                     
      VV     = V*V                                                      
      OMEGAP = TPM*V/QQ+PMPM/QQ                                         
                                                                        
C OVERCOME RISK OF UNDERFLOW IN THE EXPONENTIATION                      
      OMEGAP = DMIN1(20.0D0,OMEGAP)                                       
                                                                        
      SP = 1.0-EXP(-7.7*(OMEGAP-1.0))                                   
      IF (amuM.LE.1.5) THEN !hydrogen
C          UNIVERSAL AND RESONANCE FIT FOR HYDROGEN                     
           UNIV = SLACF2(W,QQ,CF)                                       
           BRES = B(W,QQ,C)                                             
      ELSE                                                              
C          UNIVERSAL AND RESONANCE FIT FOR DEUTERIUM                    
           UNIV = SLACF2(W,QQ,CFD)/SP
           BRES = B(W,QQ,CD)          
      ENDIF                                                             
                                                                        
C COMPUTE VW2,W2,W1                                                     
                                                                        
      VW2    = UNIV*BRES 
      IF (amuM.GE.1.5) VW2=VW2/2.  !*****  per nucleon 11/3/95   ***********
      W2     = VW2/V                                                    
      W1     = (1.0D0+VV/QQ)/(V*(1.0D0+R))*VW2                          
!      if(prttst) write(*,'(1x,''univ...='',6f10.4)') sp,univ,bres,
!     >  vw2,w2,w1
      IF (amuM.LE.2.5) RETURN                                               
      X      = QQ/2./PM/V
      EMCFAC= FITEMC((X),(amuM),GOODFIT)
cdg      EMCFAC= FITEMC(REAL(X),REAL(amuM),GOODFIT)
C$$      SUMEF  = EF(1)                                                    
C$$      DO 11 J=2,7                                                       
C$$      ZZ     = J-1.                                                     
C$$11    SUMEF  = SUMEF+EF(J)*X**ZZ                                        
C$$      EMCFAC = 1+SUMEF*DLOG(amuM)                                       
                                                                        
      W2     = W2*EMCFAC                                                
      W1     = W1*EMCFAC                                                
                                                                        
      RETURN                                                            
      END                                                               
   
                                                           
C-----------------------------------------------------------------------
                                                                        
      REAL*8 FUNCTION SLACF2(WM,QSQ,CF)                                        
                                                                        
C UNIVERSAL FUNCTION FOR ATWOOD'S FIT                                   

      Implicit none                                      
      REAL*8    WM,QSQ,CF(11)                                               
      REAL*8    PM2/1.876512/, PMSQ/.880324/, PHTR/.61993/
      REAL*8    V,OMEGA,XX,XPX,OMEGAW,ARG

                                                                        
C OMEGAW FIT...NO PHOTO-PRODUCTION COUPLING                             
                                                                        
      V      = (WM**2+QSQ-PMSQ)/PM2                                     
      OMEGA  = 2.*CF(8)*V/QSQ                                           
      XX     = 1./OMEGA                                                 
      XPX    = CF(9)+CF(10)*(XX-CF(11))**2                              
      OMEGAW = (2.D0*CF(8)*V+CF(6))/(QSQ+CF(7))                         
      ARG    = 1.-1./OMEGAW                                             
                                                                        
      SLACF2 = OMEGAW/OMEGA*ARG**3*(CF(1)+CF(2)*ARG+                    
     >         CF(3)*ARG**2+CF(4)*ARG**3+CF(5)*ARG**4)                  
      SLACF2 = SLACF2*XPX                                               
                                                                        
      RETURN                                                            
      END                           
!---------------------------------------------------------------------
      REAL FUNCTION FITEMC_N(X,A,Z,GOODFIT)                                         
!---------------------------------------------------------------------  
! Modified FITEMC.F with Neutron excess correction and proton=1 added 8/19/98
! Fit to EMC effect.  Steve Rock 8/3/94                                 
! Funciton returns value of sigma(A)/sigma(d) for isoscaler nucleus     
!  with A/2 protons=neutrons.                                           
! A= atomic number                                                      
! Z= number of protons
! x = Bjorken x.                                                        
!                                                                       
! Fit of sigma(A)/sigma(d) to form C*A**alpha where A is atomic number  
! First data at each x was fit to form C*A**alpha.  The point A=2(d)    
!  was included with a value of 1 and an error of about 2%.             
! For x>=.125 Javier Gomez fit of 7/93 to E139 data was used.           
! For .09 >=x>=.0085 NMC data from Amaudruz et al Z. Phys C. 51,387(91) 
!  Steve did the fit for alpha and C to the He/d. C/d and Ca/d NMC data.
! Alpha(x) was fit to a 9 term polynomial a0 +a1*x +a2*x**2 + a3*x**3 ..
! C(x) was fit to a 3 term polynomial in natural logs as                
!  Ln(C) = c0 + c1*Ln(x) + c2*[Ln(x)]**2.                               

! 6/2/98 *****  Bug (which set x= .00885 if x was out of range) fixed
!                    also gave value at x=.0085 if x>.88
! 8/19/98 **  If proton, return 1.
! 8/19/98 **  Add neutron excess correction.
!-----------------------------------------------------------------------
                                                                        
                                                                        
      IMPLICIT NONE                                                     
      INTEGER I                                                         
      REAL ALPHA, C,LN_C,X,A,Z ,X_U,SIG_N_P,F_IS
      REAL ATMP
      LOGICAL GOODFIT                                  
                                                                        
!Chisq=         19.   for 30 points
!Term    Coeficient     Error
      REAL*8 ALPHA_COEF(2,0:8)    /                                     
     > -6.98871401D-02,    6.965E-03,                                   
     >  2.18888887D+00,    3.792E-01,                                   
     > -2.46673765D+01,    6.302E+00,                                   
     >  1.45290967D+02,    4.763E+01,                                 
     > -4.97236711D+02,    1.920E+02,                                   
     >  1.01312929D+03,    4.401E+02,                                   
     > -1.20839250D+03,    5.753E+02,                                   
     >  7.75766802D+02,    3.991E+02,                                   
     > -2.05872410D+02,    1.140E+02 /                                  

             !
!Chisq=         22.    for 30 points
!Term    Coeficient     Error 
      REAL*8 C_COEF(2,0:2) /        ! Value and error for 6 term fit to 
     >  1.69029097D-02,    4.137D-03,                                   
     >  1.80889367D-02,    5.808D-03,                                   
     >  5.04268396D-03,    1.406D-03   /                                
                                                                        

CDG ATMP is a fake variable so I can iterate the EMC effect
      ATMP=A
      if(A.eq.4.0) then
         ATMP=5.0
      endif
      IF(A.LT.1.5) THEN    ! Added 8/19/98
       FITEMC_N=1.
       GOODFIT=.TRUE.
       RETURN
      ENDIF                                                                                
      IF( (X.GT.0.88).OR.(X.LT. 0.0085) ) THEN   !Out of range of fit   
       IF(X.LT. 0.0085) X_U =.0085
       IF(X.GT. 0.88) X_U =.88
       GOODFIT=.FALSE.
      ELSE
       X_U=X
       GOODFIT=.TRUE.
      ENDIF                                                            
                                                                        
      LN_C = C_COEF(1,0)                                                
      DO I =1,2                                                         
       LN_C = LN_C + C_COEF(1,I) * (ALOG(X_U))**I                         
      ENDDO                                                             
                                                                        
      C = EXP(LN_C)                                                     
                                                                        
      ALPHA = ALPHA_COEF(1,0)                                           
      DO I=1,8                                                          
       ALPHA = ALPHA + ALPHA_COEF(1,I) * X_U**I                           
      ENDDO                                                             
                                                                        
      FITEMC_N  =  C *ATMP**ALPHA    !isoscaler
      SIG_N_P = 1.-0.8*X_U
      F_IS = .5*(1.+ SIG_N_P)/(Z/A +(1.-Z/A)*SIG_N_P)
      FITEMC_N = FITEMC_N/F_IS
      RETURN                                                            
      END                                          
!---------------------------------------------------------------------

      REAL FUNCTION FITEMC(X,A,GOODFIT)                                         
!---------------------------------------------------------------------  
! Fit to EMC effect.  Steve Rock 8/3/94                                 
! Funciton returns value of sigma(A)/sigma(d) for isoscaler nucleus     
!  with A/2 protons=neutrons.                                           
! A= atomic number                                                      
! x = Bjorken x.                                                        
!                                                                       
! Fit of sigma(A)/sigma(d) to form C*A**alpha where A is atomic number  
! First data at each x was fit to form C*A**alpha.  The point A=2(d)    
!  was included with a value of 1 and an error of about 2%.             
! For x>=.125 Javier Gomez fit of 7/93 to E139 data was used.           
! For .09 >=x>=.0085 NMC data from Amaudruz et al Z. Phys C. 51,387(91) 
!  Steve did the fit for alpha and C to the He/d. C/d and Ca/d NMC data.
! Alpha(x) was fit to a 9 term polynomial a0 +a1*x +a2*x**2 + a3*x**3 ..
! C(x) was fit to a 3 term polynomial in natural logs as                
!  Ln(C) = c0 + c1*Ln(x) + c2*[Ln(x)]**2.                               

! 6/2/98 *****  Bug (which set x= .00885 if x was out of range) fixed
!                    also gave value at x=.0085 if x>.88
!-----------------------------------------------------------------------
                                                                        
                                                                        
      IMPLICIT NONE                                                     
      INTEGER I                                                         
      REAL*8 ALPHA, C,LN_C,X,A ,X_U
      LOGICAL GOODFIT                                  
                                                                        
!Chisq=         19.   for 30 points                                     
!Term    Coeficient     Error                                           

      REAL*8  ALPHA_COEF(2,0:8)    /                                     
     > -6.98871401D-02,    6.965D-03,                                   
     >  2.18888887D+00,    3.792D-01,                                   
     > -2.46673765D+01,    6.302D+00,                                   
     >  1.45290967D+02,    4.763D+01,                                   
     > -4.97236711D+02,    1.920D+02,                                   
     >  1.01312929D+03,    4.401D+02,                                   
     > -1.20839250D+03,    5.753D+02,                                   
     >  7.75766802D+02,    3.991D+02,                                   
     > -2.05872410D+02,    1.140D+02 /                                  
                                                     
                              
!Chisq=         22.    for 30 points                                   
!Term    Coeficient     Error                                          
      REAL*8 C_COEF(2,0:2) /        ! Value and error for 6 term fit to 
     >  1.69029097D-02,    4.137D-03,                                   
     >  1.80889367D-02,    5.808D-03,                                   
     >  5.04268396D-03,    1.406D-03   /                                
                                                                        
                                                                        
      IF( (X.GT.0.88).OR.(X.LT. 0.0085) ) THEN   !Out of range of fit   
       IF(X.LT. 0.0085) X_U =.0085
       IF(X.GT. 0.88) X_U =.88
       GOODFIT=.FALSE.
      ELSE
       X_U=X
       GOODFIT=.TRUE.
      ENDIF                                                            
                                                                        
      LN_C = C_COEF(1,0)                                                
      DO I =1,2                                                         
       LN_C = LN_C + C_COEF(1,I) * (ALOG(X_U))**I                         
      ENDDO                                                             
                                                                        
      C = EXP(LN_C)                                                     
                                                                        
      ALPHA = ALPHA_COEF(1,0)                                           
      DO I=1,8                                                          
       ALPHA = ALPHA + ALPHA_COEF(1,I) * X_U**I                           
      ENDDO                                                             
                                                                        
      FITEMC  =  C *A**ALPHA                                            
      RETURN                                                            
      END                                                               
C-----------------------------------------------------------------------
                                                                        
      REAL*8 FUNCTION B(WM,QSQ,C)                                              
                                                                        
C BACKGROUND AND RESONANCE CONTRIBUTION FOR ATWOOD'S FIT                

      Implicit none
      REAL*8  WM,QSQ,C(24),WSQ,OMEGA,X,XPX,PIEMSQ,B1,EB1,B2,BBKG,BRES
      REAL*8  RAM,RMA,RWD,QSTARN,QSTARO,TERM,TERMO,GAMRES,BRWIG,RES
      REAL*8  RESSUM,EB2,ressv(4)
      INTEGER   LSPIN(4),INDEX,J,K                                             
      REAL*8    PMSQ/.880324/, PM2/1.876512/, PM/.938256/            
      INTEGER   NRES/4/, NBKG/5/,I                                     
      common/testing/prttst
      logical prttst
      DATA      LSPIN/1,2,3,2/                                       
                                                                        
C KINEMATICS                                                            
                                                                        
      WSQ    = WM**2                                                    
      OMEGA  = 1.+(WSQ-PMSQ)/QSQ                                        
      X      = 1./OMEGA                                                 
      XPX    = C(22)+C(23)*(X-C(24))**2                                 
      PIEMSQ = (C(1)-PM)**2                                             
                                                                        
C COLLECT BACKGROUND TERMS AND CHECK FOR EXPONENTIAL UNDERFLOWS BEFORE  
C THEY HAPPEN                                                           
                                                                        
      B1 = 0.                                                           
      IF (WM.GT.C(1)) B1 = C(2)                                         
      EB1 = C(3)*(WM-C(1))                                              
      IF (EB1.LE.25.) B1 = B1*(1.-EXP(-EB1))                            
      B2 = 0.                                                           
      IF (WM.GT.C(4)) B2 = (1.-C(2))                                    
      EB2 = C(5)*(WSQ-C(4)**2)                                          
      IF (EB2.LE.25.0) B2 = B2*(1.-EXP(-EB2))                           
      BBKG = B1+B2                                                      
      BRES = C(2)+B2                                                    
                                                                        
C COLLECT RES. CONTRIBUTION                                             
                                                                        
      RESSUM = 0.                                                       
      DO 30 I=1,NRES                                                    
           INDEX  = (I-1)*3+1+NBKG                                      
           RAM    = C(INDEX)                                            
           IF (I.EQ.1) RAM=C(INDEX)+C(18)*QSQ+C(19)*QSQ**2              
           RMA    = C(INDEX+1)                                          
           IF (I.EQ.3) RMA=RMA*(1.D0+C(20)/(1.D0+C(21)*QSQ))            
           RWD    = C(INDEX+2)                                          
           QSTARN =SQRT(DMAX1(0.D0,((WSQ+PMSQ-PIEMSQ)/(2.*WM))**2-PMSQ)) 
           QSTARO = SQRT(DMAX1(0.D0,((RMA**2-PMSQ+PIEMSQ)/                
     >              (2.*RMA))**2-PIEMSQ))                               
                                                                        
           RES = 0.                                                     
           IF (QSTARO.NE.0.) THEN                                       
                TERM   = 6.08974*QSTARN                                 
                TERMO  = 6.08974*QSTARO                                 
                J      = 2*LSPIN(I)                                     
                K      = J+1                                            
                GAMRES = RWD*(TERM/TERMO)**K*(1.+TERMO**J)/(1.+TERM**J) 
                GAMRES = GAMRES/2.                                      
                BRWIG  = GAMRES/((WM-RMA)**2+GAMRES**2)/3.1415926       
                RES    = RAM*BRWIG/PM2                                  
           ENDIF                                                        
           ressv(i)=res
           RESSUM = RESSUM+RES                                          
30    CONTINUE                                                          
      if(prttst) write(*,'(1x,''w,q2,res='',6f7.3)') wm,qsq,
     >  ressv
                                                                   
C FORM VW2/F2                                                           
                                                                        
      B = BBKG*(1.+(1.-BBKG)*XPX)+RESSUM*(1.-BRES)                      
!      if(prttst) write(*,'(1x,''b...'',6f10.5)') b,bbkg,xpx,ressum                                                                  
      RETURN                                                            
      END                                                               
                  
