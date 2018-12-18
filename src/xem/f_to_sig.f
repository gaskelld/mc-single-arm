	subroutine f_to_sig(e1,e2,theta,a,z,m1,aux,sig,y,fact)
C+______________________________________________________________________________
!
! Subroutine to compute sigma (nb/sr/MeV) from Y-scaling model. This subroutine
!   is derived from the subroutine FTOSIG by R. McKeown. It uses the same
!   basic method, but calculates the scaling variable Y from a different
!   formula, in which the separation energy is used to calculate the mass
!   of the recoiling (A-1) system. F_TO_SIG also has an argument list
!   different from FTOSIG, and hence is NOT interchangeable with FTOSIG.
!
! INPUT ARGUMENTS:
!
!   E1:		(R*8) - Incident energy in GeV
!   E2:		(R*8) - Outgoing energy in GeV
!   THETA:	(R*8) - Scattering angle in degrees
!   A:		(R*8) - 'A' of target nucleus
!   Z:		(R*8) - 'Z' of target nucleus
!   M1:		(R*8) - Mass of target nucleus in GeV/c2.
!   MR:		(R*8) - Mass of recoiling nucleon in GeV/c2.
!   ES:		(R*8) - Separation energy in GeV.
!
! OUTPUT ARGUMENTS:
!
!   SIG:	(R*8) - Calculated cross section
!   Y:		(R*8) - Scaling variable Y in 1/GeV.
!
!   J. Arrington, 3/5/98
!
!   AUX(1) = RESOL parameter for DEEPSIG smearing (not used).
!   AUX(2) = ES (Separation energy)
!   AUX(3) = F(0)  parameter for F(y) function
!   AUX(4) = BigB  parameter for F(y)
!   AUX(5) = a     parameter for F(y)
!   AUX(6) = b     parameter for F(y)
!   AUX(7) = alpha parameter for F(y)
!
!   D. Potterveld, 9/16/86 - AUX DEFINITIONS FROM PREVIOUS F(y) MODEL>
!
!  ! AUX(1) = MR (Recoil mass)
!  ! AUX(2) = ES (Separation energy)
!  ! AUX(3) = Y0 parameter for F(y) function
!  ! AUX(4) = B  parameter for F(y)
!  ! AUX(5) = C  parameter for F(y)
!  ! AUX(6) = D  parameter for F(y)
!  ! AUX(7) = RESOL parameter for DEEPSIG smearing.
!
! D. Potterveld, 11-24-86
C-______________________________________________________________________________

	implicit none

C Math/physics constants.

	include 'include.inc'
        real*8 pi,d_r,mp,mp_sq
	parameter (pi = 3.141592654)
	parameter (d_r = pi/180.)
	parameter (mp = .93827231)		!Proton mass in GeV/c2.
	parameter (mp_sq = mp*mp)

C arguments.

	real*8 e1,e2,theta,a,z,m1,mr,es,sig,y
	real*8 aux(7)

C Local variables.

	real*8 th,nu,q3,m2,tau,dwdy,sig_mott
	real*8 fact
	real*8 f1p,f2p,f1n,f2n,WMA,kcm,kmin,Einit,Efinal
	real*8 prefac,qbar2,taubar,sig_p,sig_n
	real*8 ge2,gm2,w20,w10,pf,pp2,g2,f2,w2,w1

C Declare arguments for call to Hoehler(nform).

	real*8 q4_sq,gep,gmp,gen,gmn,temp,corfact

C Declare function data types.

	logical y_calc_ok
	real*8  fy,x_local, my_frac,x_high,x_low,johns_fac
	save

C ============================ Executable Code =================================

!!	mr = aux(1)	  		!Copy (R*8) args to int. storage
	mr = mp
	es = aux(2)

C Compute q2 and other items.
                            
	th = theta * d_r			!Angle in radians
  	q4_sq = 4.*e1*e2*sin(th/2.)**2 		!4-mom. transf. squared.
	nu = e1-e2				!Energy loss
	q3 = sqrt(q4_sq + nu**2)		!3-momentum transfer
	m2 = m1 - mr + es			!Mass of A-1 system.
      	tau = q4_sq/(4.*mp_sq)

	x_local=q4_sq/(2*mp*nu)



C Mott cross section.

	sig_mott = cos(th/2.)/(274.072*e1*sin(th/2.)**2.)
	sig_mott = (.019733*sig_mott)**2

!- rewritten starting here--------------------------------------

        call y_calc(e1,e2,theta,m1,mr,es,y,y_calc_ok)
	if (.not.y_calc_ok) then
	   sig = 0.
	   y = 1.0E20
	   goto 900
	endif

	if(a.gt.1.5) then
	   call sig_bar_df_dble(e1, e2, theta, y, 0.0, sig_p, sig_n)
	   sig_p=sig_p/1000.
	   sig_n=sig_n/1000.
	endif
C Calculate scaling variable.

	dwdy = q3/sqrt(mp**2+q3**2+y**2+2.0*q3*y)

c	dwdy = (q3+y)/sqrt(mp**2+q3**2+y**2+2.0*q3*y) + y/sqrt(m2**2+y**2)


	if(a.gt.1.5) then
	   fact = (Z*sig_p+(A-Z)*sig_n)/dwdy
	else
	   fact = 0.0
	endif
	sig = fact*fy(y,aux(3),aux(4),aux(5),aux(6),aux(7), A)*1.D6 !new f(y) fit

         johns_fac=1.0
        if (a.eq.3) then
          johns_fac=max(1.,1.+y*1.4*2.5)
        elseif (a.eq.4) then
          johns_fac=max(1.,1.+y*1.4*1.75)
        elseif ((a.eq.9).or.(a.eq.12)) then  
          johns_fac=max(1.,1.+y*1.4*2.5)
        elseif ((a.eq.56).or.(a.eq.64)) then
          johns_fac=max(1.,1.+y*1.4*3)
        elseif (a.eq.197) then
          johns_fac=max(1.,1.+y*1.4*4)
        endif
        sig=johns_fac*sig

c	if(x_local.gt.3.) then	! this is to avoid the negative cor factor; aji
c	   x_local=3.
c	endif
c
c         if(a.eq.2) then
c	  x_low=1.1
c	  x_high=1.2
c	endif
c
c	if(a.gt.2) then
c	  if((a.eq.64).or.(a.eq.4).or.(a.eq.197)) then
c	     x_low=1.2
c	     x_high=1.4
c	  else
c	     x_low=1.4
c	     x_high=1.6
c	  endif
c	endif
c
c	if(x_local.gt.2.0) then
c	  x_local=2.0
c	endif
c
c	if(a.eq.2) then
c	  x_low=1.
c	  x_high=1.05
c	endif
c
c	if(a.gt.2) then
c	  if((a.eq.64).or.(a.eq.4).or.(a.eq.197)) then
c	    x_low=1.2
c	    x_high=1.4
c	  else
c	    x_low=1.4
c	    x_high=1.6
c	  endif
c	endif
	
	if(x_local.gt.2.0) then
	  x_local=2.0
	endif

	if(a.eq.2) then
	  x_low=1.4
	  x_high=1.45
	endif

	if(a.gt.2) then
	  if((a.eq.64).or.(a.eq.4).or.(a.eq.197)) then
	    x_low=1.2
	    x_high=1.4
	  else
	    x_low=1.4
	    x_high=1.6
	  endif
	endif



	if((x_local.ge.(x_low))) then !.and.(A.eq.2)) then
	   corfact=(aa*exp(bb*x_local)+cc*x_local**6+dd*x_local**4+ee*x_local**2+ff)
	   if(x_local.lt.(x_high)) then
	      my_frac=(x_local-x_low)/(x_high-x_low)	      
	      corfact=my_frac*corfact+(1.-my_frac)
	   endif
	   sig=sig*corfact	   
	endif



900	continue 

	return
	end

