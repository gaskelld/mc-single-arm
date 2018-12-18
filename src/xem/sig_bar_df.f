	subroutine sig_bar_df(e1,e2,theta,pl,pt,sig_p,sig_n)
C+______________________________________________________________________________
!
! DESCRIPTION:
!
!   Calculate the DeForest based off-shell cross sections averaged over the
!   angle PHI, which is the angle between the electron scattering plane and
!   the plane containing the initial and final momentum of the struck nucleon.
!
!   We use CIT form factors, based on the four-momentum transfer.
!
!   All energy/momentum variables are in GeV; cross sections are in mB/ster.
!
! ARGUMENTS: (All R*4 data type)
!
!   E1:    Energy of incident electron.
!   E2:    Energy of scattered electron.
!   THETA: Electron scattering angle in deg.
!   PL:    Projection of init. state nucleon's momentum along 3-mom. transfer.
!   PT:    Component of init. state nucleons momentum perp. to 3-mom. trnsfr.
!   SIG_P: Cross section for proton.  (OUTPUT)
!   SIG_N: Cross section for neutron. (OUTPUT)
C-______________________________________________________________________________

	implicit real*4(a-z)

	include		'math_physics.inc'

C Argument data types.

	real*4		e1,e2,theta,pl,pt,sig_p,sig_n

C ============================ Executable Code =================================

C Compute some kinematics.

1	th = theta*d_r
	tan_2 = tan(th/2)**2
	nu = e1 - e2
	q4_2 = 4*e1*e2*sin(th/2)**2
	tau = q4_2/(4*m_p**2)
	qv_2 = q4_2 + nu**2
	qv = sqrt(qv_2)
	pt_2 = pt*pt
	p_2 = pl**2 + pt_2

! Mott cross section in mB/sR.

      	sig_mott = hc_2*(alpha*cos(th/2))**2/(2*e1*sin(th/2)**2)**2

! Final energy of struck nucleon, assuming it is on shell.

	en_f = sqrt(m_p**2+qv_2+p_2+2*qv*pl)

C 'BAR'-ed quantities of DeForest.

	e_bar = sqrt(p_2 + m_p**2)
	q4_bar_2 = qv_2 - (en_f - e_bar)**2
	tau_bar = q4_bar_2/(4*m_p**2)

C Get form factors. Convert to F1, F2.

!	call nuc_form(q4_2,gep,gmp,gen,gmn)

!	write(6,*) 'q4_2,gep,gen,gmp,gmn=',q4_2,gep,gen,gmp,gmn

	call nform_xem (12.0,q4_2,gep,gen,gmp,gmn)
c	write(6,*) 'q4_2,gep,gen,gmp,gmn=',q4_2,gep,gen,gmp,gmn
	tmp = gmp

	call nform_xem (1.0,q4_2,gep,gen,gmp,gmn)
	gmp = tmp
c	write(6,*) 'q4_2,gep,gen,gmp,gmn=',q4_2,gep,gen,gmp,gmn

	f1p = (gep + tau*gmp)/(1 + tau)
	f1n = (gen + tau*gmn)/(1 + tau)
	f2p = (gmp - gep)/(1 + tau)
	f2n = (gmn - gen)/(1 + tau)

C DeForest cross sections.

	qrat = q4_2/qv_2
	t1 = q4_bar_2*tan_2/2 + .25*qrat*(q4_bar_2 - q4_2)
	t2 = .25*(qrat**2)*(e_bar + en_f)**2 + (qrat/2 + tan_2)*pt_2
	sig_p = (t1*(f1p + f2p)**2 + t2*(f1p**2 + tau_bar*f2p**2))
	sig_n = (t1*(f1n + f2n)**2 + t2*(f1n**2 + tau_bar*f2n**2))
	sig_p = sig_mott*sig_p/(en_f*e_bar)
	sig_n = sig_mott*sig_n/(en_f*e_bar)

	return
	end

C ############################ SIG_BAR_DF_DBLE #################################

	subroutine sig_bar_df_dble(e1,e2,theta,pl,pt,sig_p,sig_n)
C+______________________________________________________________________________
!
! Double precision version of above.
C-______________________________________________________________________________

	implicit	real*8(a-z)
	include		'math_physics.inc'

C Argument data types.

	real*8		e1,e2,theta,pl,pt,sig_p,sig_n
	real*8		gep,gmp,gen,gmn

C ============================ Executable Code =================================

C Compute some kinematics.
!	write (*,*) 'stuff passed into john', e1, e2, theta, pl, pt

1	th = theta*d_r
 	tan_2 = tan(th/2)**2
	nu = e1 - e2
	q4_2 = 4*e1*e2*sin(th/2)**2
	tau = q4_2/(4*m_p**2)
	qv_2 = q4_2 + nu**2
	qv = sqrt(qv_2)
	pt_2 = pt*pt
	p_2 = pl**2 + pt_2

! Mott cross section in mB/sR.

      	sig_mott = hc_2*(alpha*cos(th/2))**2/(2*e1*sin(th/2)**2)**2

!	write(*,*) 'sig mott is ', sig_mott

! Final energy of struck nucleon, assuming it is on shell.

	en_f = sqrt(m_p**2+qv_2+p_2+2*qv*pl)

C 'BAR'-ed quantities of DeForest.

	e_bar = sqrt(p_2 + m_p**2)
	q4_bar_2 = qv_2 - (en_f - e_bar)**2
	tau_bar = q4_bar_2/(4*m_p**2)

!	write(*,*) 'next set ', en_f, e_bar, q4_bar_2, tau_bar, q4_2
C Get form factors. Convert to F1, F2.

!	call nuc_form(real(q4_2),gep,gmp,gen,gmn)

!	gep=0.0
!	gen=0.0
!	gmp=0.0
!	gmn=0.0
	q4_sq=q4_2
        call nform_xem (12.0,q4_sq,gep,gen,gmp,gmn)
        tmp = gmp
c	write(6,*) 'q4_2,gep,gen,gmp,gmn=',q4_2,gep,gen,gmp,gmn
!	write(*,*) 'ff in john ', gep, gen, gmp, gmn
        call nform_xem (1.0,q4_sq,gep,gen,gmp,gmn)
        gmp = tmp

!	write(6,*) 'q4_2,gep,gen,gmp,gmn=',q4_2,gep,gen,gmp,gmn
!	write(*,*) 'ff in john 2 ', gep, gen, gmp, gmn
	f1p = (gep + tau*gmp)/(1 + tau)
	f1n = (gen + tau*gmn)/(1 + tau)
	f2p = (gmp - gep)/(1 + tau)
	f2n = (gmn - gen)/(1 + tau)
	
!	write(*,*) 'others are ', f1p, f1n, f2p, f2n

C DeForest cross sections.

	qrat = q4_2/qv_2
	t1 = q4_bar_2*tan_2/2 + .25*qrat*(q4_bar_2 - q4_2)
	t2 = .25*(qrat**2)*(e_bar + en_f)**2 + (qrat/2 + tan_2)*pt_2
	sig_p = (t1*(f1p + f2p)**2 + t2*(f1p**2 + tau_bar*f2p**2))
	sig_n = (t1*(f1n + f2n)**2 + t2*(f1n**2 + tau_bar*f2n**2))
	sig_p = sig_mott*sig_p/(en_f*e_bar)
	sig_n = sig_mott*sig_n/(en_f*e_bar)
	
!	write(*,*) 'about to return ', sig_p, sig_n

	return
	end
