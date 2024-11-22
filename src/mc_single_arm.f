	program mc_single_arm

C+______________________________________________________________________________
!
! Monte-Carlo of SHMS spectrometer using uniform illumination.
!   This version uses TRANSPORT right-handed coordinate system.
C-______________________________________________________________________________

	implicit none

	include 'hms/struct_hms.inc'
	include 'shms/struct_shms.inc'
	include 'spectrometers.inc'
	include 'constants.inc'
	include 'hbook.inc'

c Vector (real*4) for hut ntuples - needs to match dimension of variables
	real*8		shms_hut(23)
	real*8          shms_spec(59)

	real*8          hms_hut(23)
c
	real*8 xs_num,ys_num,xc_sieve,yc_sieve
	real*8 xsfr_num,ysfr_num,xc_frsieve,yc_frsieve
        logical use_front_sieve /.false./
        logical use_sieve /.false./           
c
        common /sieve_info/  xs_num,ys_num,xc_sieve,yc_sieve
     > ,xsfr_num,ysfr_num,xc_frsieve,yc_frsieve,use_sieve, use_front_sieve


C Local declarations.
	integer*4	i,
     >			chanin	/1/,
     >			chanout	/2/,
     >			n_trials,trial,
     >			tmp_int

	integer*4 Itrial                        ! TH - add this for gfortran: forces integer type cast
	logical*4	iss

	real*8 th_nsig_max                      ! TH - add this for gfortran
	parameter(th_nsig_max=3.0d0)            !max #/sigma for gaussian ran #s

C Event limits, topdrawer limits, physics quantities
	real*8 gen_lim(8)			!M.C. phase space limits.
	real*8 gen_lim_up(3)
	real*8 gen_lim_down(3)
        real*8 gen_mom                         !local variable for track momentum in elastic event

	real*8 cut_dpp,cut_dth,cut_dph,cut_z	!cuts on reconstructed quantities
	real*8 xoff,yoff,zoff                   !Beam offsets
        real*8 spec_xoff,spec_yoff,spec_zoff    !Spectrometer offsets
	real*8 spec_xpoff, spec_ypoff           !Spectrometer angle offsets
	real*8 th_ev,cos_ev,sin_ev		!cos and sin of event angle

	real*8 cos_ts,sin_ts			!cos and sin of spectrometer angle
	real*8 x,y,z,dpp,dxdz,dydz,t1,t2,t3,t4	!temporaries

	real*8 x_a,y_a,z_a,dydz_a,dif_a,dydz_aa,dif_aa   ! TH - for target aperture check
	real*8 musc_targ_len			!target length for multiple scattering
        real*8 foil_nm, foil_tk                 !multifoil target
	real*8 foil_zcent
        parameter (foil_tk=0.02)
	real*8 m2				!particle mass squared.
	real*8 rad_len_cm			!conversion r.l. to cm for target
	real*8 pathlen				!path length through spectrometer.
	logical*4 ok_spec			!indicates whether event makes it in MC
	integer*4 hit_calo                      !flag for hitting the calorimeter
	integer*4 armSTOP_successes,armSTOP_trials
        real *8 beam_energy, el_energy, theta_sc!elastic calibration
        real *8 tar_mass, tar_atom_num          !elastic calibration
C Initial and reconstructed track quantities.
	real*8 dpp_init,dth_init,dph_init,xtar_init,ytar_init,ztar_init
	real*8 dpp_recon,dth_recon,dph_recon,ztar_recon,ytar_recon
	real*8 x_fp,y_fp,dx_fp,dy_fp		!at focal plane
	real*8 fry,fr1,fr2
	real*8 p_spec,th_spec,phi_spec			!spectrometer setting
	real*8 resmult

	real*8 ptot
	real*8 px_gen,py_gen,pz_gen,px_vert,py_vert,pz_vert,zgen
	real*8 ux,uy,uz,ux0,uy0,uz0,cos_dtheta,y_event
	real*8 pe_recon,px_spec,py_spec,pz_spec,px_lab,py_lab,pz_lab
	real*8 hms_stop_id
	real*8 bproton,bair,bal,bkevlar,bmylar,bt_tot
	real*8 Lair,Lal,Lkevlar,Lmylar
	real*8 ranphot,deltaE

C Control flags (from input file)
	integer*4 ispec
	integer*4 p_flag			!particle identification
	logical*4 ms_flag
	logical*4 wcs_flag
	logical*4 store_all

c	common /hutflag/ cer_flag,vac_flag
C Hardwired control flags.
	logical*4 hut_ntuple	/.true./
        logical*4 spec_ntuple   /.false./
	logical*4 decay_flag	/.false./

	real*8	dpp_var(2),dth_var(2),dph_var(2),ztg_var(2)
	integer*8	stime,etime

	character*132	str_line
C Local  spectrometer varibales
	real*8 x_s,y_s,z_s
	real*8 dxdz_s,dydz_s,dpp_s

C Function definitions.

	integer*4	last_char
	logical*4	rd_int,rd_real
	real*8          grnd,gauss1
	INTEGER irnd
	REAL rnd(99)
        integer      itime,ij
        character	timestring*30

        character*80 rawname, filename, hbook_filename
	real*4  secnds,zero

	parameter(zero=0.0)

	integer ivar

	integer GetNumBranches
	external GetNumBranches
	integer GetNtEntries
	external GetNtEntries 

	save		!Remember it all!

C ================================ Executable Code =============================

C Initialize
C using SIMC unstructured version
C
C SHMS
	shmsSTOP_trials	= 0
        shmsSTOP_targ_hor	= 0
        shmsSTOP_targ_vert	= 0
        shmsSTOP_targ_oct	= 0
        shmsSTOP_FRONTSLIT_hor	= 0
        shmsSTOP_FRONTSLIT_vert	= 0
	shmsSTOP_HB_in	= 0
        shmsSTOP_HB_men = 0
	shmsSTOP_HB_mex	= 0
	shmsSTOP_HB_out	= 0	
	shmsSTOP_DOWNSLIT	= 0
	shmsSTOP_COLL_hor	= 0
	shmsSTOP_COLL_vert	= 0
	shmsSTOP_COLL_oct	= 0
	shmsSTOP_Q1_in	= 0
	shmsSTOP_Q1_men	= 0
	shmsSTOP_Q1_mid	= 0
	shmsSTOP_Q1_mex	= 0
	shmsSTOP_Q1_out	= 0
	shmsSTOP_Q2_in	= 0
	shmsSTOP_Q2_men	= 0
	shmsSTOP_Q2_mid	= 0
	shmsSTOP_Q2_mex	= 0
	shmsSTOP_Q2_out	= 0
	shmsSTOP_Q3_in	= 0
	shmsSTOP_Q3_men	= 0
	shmsSTOP_Q3_mid	= 0
	shmsSTOP_Q3_mex	= 0
	shmsSTOP_Q3_out	= 0
c	shmsSTOP_Q3_out1	= 0
c	shmsSTOP_Q3_out2	= 0
c	shmsSTOP_Q3_out3	= 0
c	shmsSTOP_Q3_out4	= 0
c	shmsSTOP_Q3_out5	= 0
c	shmsSTOP_Q3_out6	= 0
	shmsSTOP_D1_in	= 0
        shmsSTOP_D1_flr = 0
	shmsSTOP_D1_men = 0
	shmsSTOP_D1_mid1 = 0
	shmsSTOP_D1_mid2 = 0
	shmsSTOP_D1_mid3 = 0
	shmsSTOP_D1_mid4 = 0
	shmsSTOP_D1_mid5 = 0
	shmsSTOP_D1_mid6 = 0
	shmsSTOP_D1_mid7 = 0
	shmsSTOP_D1_mex = 0
	shmsSTOP_D1_out	= 0
	shmsSTOP_BP_in  = 0
	shmsSTOP_BP_out = 0
	shmsSTOP_hut	= 0
	shmsSTOP_dc1	= 0
	shmsSTOP_dc2	= 0
	shmsSTOP_s1	= 0
	shmsSTOP_s2	= 0
	shmsSTOP_s3	= 0
	shmsSTOP_cal	= 0
	shmsSTOP_successes	= 0
	shmsSTOP_id = 0
C HMS
	hSTOP_trials	= 0
	hSTOP_fAper_hor	= 0
	hSTOP_fAper_vert= 0
	hSTOP_fAper_oct	= 0
	hSTOP_bAper_hor	= 0
	hSTOP_bAper_vert= 0
	hSTOP_bAper_oct	= 0
	hSTOP_slit	= 0
	hSTOP_Q1_in	= 0
	hSTOP_Q1_mid	= 0
	hSTOP_Q1_out	= 0
	hSTOP_Q2_in	= 0
	hSTOP_Q2_mid	= 0
	hSTOP_Q2_out	= 0
	hSTOP_Q3_in	= 0
	hSTOP_Q3_mid	= 0
	hSTOP_Q3_out	= 0
	hSTOP_D1_in	= 0
	hSTOP_D1_out	= 0
	hSTOP_VP1	= 0
	hSTOP_VP2	= 0
	hSTOP_VP3	= 0
	hSTOP_VP4	= 0
	hSTOP_hut	= 0
	hSTOP_dc1	= 0
	hSTOP_dc2	= 0
	hSTOP_scin	= 0
	hSTOP_cal	= 0
	hSTOP_successes	= 0
	hSTOP_id        = 0

C Open setup file.

	write(*,*)'Enter input filename (assumed to be in infiles dir)'
	read(*,1968) rawname
 1968	format(a)
	filename = '../infiles/'//rawname(1:last_char(rawname))//'.inp'
	print *,filename,'opened'
	open(unit=chanin,status='old',file=filename)

C Define HBOOK/NTUPLE filename if used.
c	if (hut_ntuple) then
c	  hbook_filename = '../worksim/'//rawname(1:last_char(rawname))//'.bin'
c	endif
C Open Output file.
	filename = '../outfiles/'//rawname(1:last_char(rawname))//'.out'
	open (unit=chanout,status='unknown',file=filename)
	
C Read in real*8's from setup file

	str_line = '!'

C Strip off header

	do while (str_line(1:1).eq.'!')
	  read (chanin,1001) str_line
	enddo

! Read data lines.

c	write(*,*),str_line(1:last_char(str_line))
c	iss = rd_int(str_line,n_trials)
c	if (.not.iss) stop 'ERROR (ntrials) in setup!'

! Spectrometer flag:
c	read (chanin,1001) str_line
	write(*,*) str_line(1:last_char(str_line))
	iss = rd_int(str_line,ispec)
	if (.not.iss) stop 'ERROR (Spectrometer selection) in setup!'
! Open HBOOK/NTUPLE file here
c	if(hut_ntuple) then
c	   if(ispec.eq.2) then
c	      call shms_hbook_init(hbook_filename,spec_ntuple)
c	   elseif(ispec.eq.1) then
c	      call hms_hbook_init(hbook_filename,spec_ntuple)
c	   else
c	      write(6,*) 'Uknown spectrometer, stopping.'
c	      stop
c	   endif
c	endif

! Spectrometer momentum:
	read (chanin,1001) str_line
	write(*,*) str_line(1:last_char(str_line))
	iss = rd_real(str_line,p_spec)
	if (.not.iss) stop 'ERROR (Spec momentum) in setup!'

! Spectrometer angle:
	read (chanin,1001) str_line
	write(*,*) str_line(1:last_char(str_line))
	iss = rd_real(str_line,th_spec)
	if (.not.iss) stop 'ERROR (Spec theta) in setup!'
	th_spec = abs(th_spec) / degrad
	cos_ts = cos(th_spec)
	sin_ts = sin(th_spec)

! M.C. limits (half width's for dp,th,ph, full width's for x,y,z)
C DJG: should not be needed - events already generated
c	do i=1,3
c	  read (chanin,1001) str_line
c	  write(*,*),str_line(1:last_char(str_line))
c	  iss = rd_real(str_line,gen_lim(i))
c	  if (.not.iss) stop 'ERROR (M.C. limits) in setup!'
c	  gen_lim_down(i) = gen_lim(i)
c	  read (chanin,1001) str_line
c	  write(*,*),str_line(1:last_char(str_line))
c	  iss = rd_real(str_line,gen_lim(i))
c	  if (.not.iss) stop 'ERROR (M.C. limits) in setup!'
c	  gen_lim_up(i) = gen_lim(i)
c	enddo
c
	do i = 4,6
	  read (chanin,1001) str_line
	  write(*,*) str_line(1:last_char(str_line))
	  iss = rd_real(str_line,gen_lim(i))
	  if (.not.iss) stop 'ERROR (M.C. limits) in setup!'
	enddo
c
c! Raster size
	do i=7,8
	   read (chanin,1001) str_line
	   write(*,*) str_line(1:last_char(str_line))
	   iss = rd_real(str_line,gen_lim(i))
	   if (.not.iss) stop 'ERROR (Fast Raster) in setup'
	enddo
c
c! Cuts on reconstructed quantities
c	read (chanin,1001) str_line
c	write(*,*),str_line(1:last_char(str_line))
c	if (.not.rd_real(str_line,cut_dpp)) 
c     > stop 'ERROR (CUT_DPP) in setup!'
c
c	read (chanin,1001) str_line
c	write(*,*),str_line(1:last_char(str_line))
c	if (.not.rd_real(str_line,cut_dth)) 
c     > stop 'ERROR (CUT_DTH) in setup!'
c
c	read (chanin,1001) str_line
c	write(*,*),str_line(1:last_char(str_line))
c	if (.not.rd_real(str_line,cut_dph)) 
c     > stop 'ERROR (CUT_DPH) in setup!'
c
c	read (chanin,1001) str_line
c	write(*,*),str_line(1:last_char(str_line))
c	if (.not.rd_real(str_line,cut_z)) 
c     > stop 'ERROR (CUT_Z) in setup!'
c
! Read in radiation length of target material in cm
	read (chanin,1001) str_line
	write(*,*) str_line(1:last_char(str_line))
	if (.not.rd_real(str_line,rad_len_cm)) 
     > stop 'ERROR (RAD_LEN_CM) in setup!'

! Beam and target offsets
	read (chanin, 1001) str_line
	write(*,*) str_line(1:last_char(str_line))
	iss = rd_real(str_line,xoff)
	if(.not.iss) stop 'ERROR (xoff) in setup!'

	read (chanin, 1001) str_line
	write(*,*) str_line(1:last_char(str_line))
	iss = rd_real(str_line,yoff)
	if(.not.iss) stop 'ERROR (yoff) in setup!'

	read (chanin, 1001) str_line
	write(*,*) str_line(1:last_char(str_line))
	iss = rd_real(str_line,zoff)
	if(.not.iss) stop 'ERROR (zoff) in setup!'

! Spectrometer offsets
	read (chanin, 1001) str_line
	write(*,*) str_line(1:last_char(str_line))
	iss = rd_real(str_line,spec_xoff)
	if(.not.iss) stop 'ERROR (spect. xoff) in setup!'

	read (chanin, 1001) str_line
	write(*,*) str_line(1:last_char(str_line))
	iss = rd_real(str_line,spec_yoff)
	if(.not.iss) stop 'ERROR (spect. yoff) in setup!'

	read (chanin, 1001) str_line
	write(*,*) str_line(1:last_char(str_line))
	iss = rd_real(str_line,spec_zoff)
	if(.not.iss) stop 'ERROR (spect. zoff) in setup!'

	read (chanin, 1001) str_line
	write(*,*) str_line(1:last_char(str_line))
	iss = rd_real(str_line,spec_xpoff)
	if(.not.iss) stop 'ERROR (spect. xpoff) in setup!'

	read (chanin, 1001) str_line
	write(*,*) str_line(1:last_char(str_line))
	iss = rd_real(str_line,spec_ypoff)
	if(.not.iss) stop 'ERROR (spect. ypoff) in setup!'

! read in flag for particle type.
	read (chanin,1001) str_line
	write(*,*) str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,p_flag)) 
     > stop 'ERROR: p_flag in setup file!'


! Read in flag for multiple scattering.
	read (chanin,1001) str_line
	write(*,*) str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,tmp_int)) 
     > stop 'ERROR: ms_flag in setup file!'
	if (tmp_int.eq.1) ms_flag = .true.

! Read in flag for wire chamber smearing.
	read (chanin,1001) str_line
	write(*,*) str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,tmp_int)) 
     > stop 'ERROR: wcs_flag in setup file!'
	if (tmp_int.eq.1) wcs_flag = .true.

! Read in flag to keep all events - success or not
	read (chanin,1001) str_line
C	write(*,*),str_line(1:last_char(str_line))
	write(*,*) str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,tmp_int)) 
     > stop 'ERROR: store_all in setup file!'
	if (tmp_int.eq.1) store_all = .true.

!     Read in flag for 'beam energy(MeV)' to trigger on elastic event if present
      beam_energy=-0.1  !by default do not use elastic event generator
      tar_atom_num=12.  !by default it is carbon
c      read (chanin,1001,end=1000,err=1000) str_line
c      write(*,*),str_line(1:last_char(str_line))
c      iss = rd_real(str_line,beam_energy)
      
! Read in flag to use sieve
c	read (chanin,1001,end=1000,err=1000) str_line
c	write(*,*) str_line(1:last_char(str_line))
c	if (.not.rd_int(str_line,tmp_int)) 
c     > stop 'ERROR: use_sieve in setup file!'
c	if (tmp_int.eq.1) then
c	  if (ispec.eq.1) use_sieve=.true.
c	  if (ispec.eq.2) use_sieve=.true.
c	  if (ispec.eq.2) use_front_sieve=.false.
c	endif
	use_sieve=.false.
	
!     Read in flag for 'target atomic number (Z+N)' for elastic event if present
c      read (chanin,1001,end=1000,err=1000) str_line
c      write(*,*) str_line(1:last_char(str_line))
c      iss = rd_real(str_line,tar_atom_num)


 1000	continue

	write(*,*)'Enter root tree name with generated events (assumed to be in generated_events dir)'
	read(*,1969) rawname
 1969	format(a)
	filename = '../generated_events/'//rawname(1:last_char(rawname))//'.root'
	call InitRootNT(filename,'UPDATE')
	write(6,*) 'opening root tree with generated events'
	write(6,*) GetNumBranches(),' branches'

	
C Set particle masses.
	m2 = me2			!default to electron
	if(p_flag.eq.0) then
	  m2 = me2
	else if(p_flag.eq.1) then
	  m2 = mp2
	else if(p_flag.eq.2) then
	  m2 = md2
	else if(p_flag.eq.3) then
	  m2 = mpi2
	else if(p_flag.eq.4) then
	  m2 = mk2
	endif

C------------------------------------------------------------------------------C
C                           Top of Monte-Carlo loop                            C
C------------------------------------------------------------------------------C
        call TagNtBranch(px_gen,'GSE_px ')
	call TagNtBranch(py_gen,'GSE_py ')
	call TagNtBranch(pz_gen,'GSE_pz ')
	call TagNtBranch(zgen,'GV_z ')
	call AddNtBranch(dpp_recon,'hsdelta')
	call AddNtBranch(dph_recon,'hsxptar')
	call AddNtBranch(dth_recon,'hsyptar')
	call AddNtBranch(ytar_recon,'hsytar')
	call AddNtBranch(x_fp,'hsxfp')
	call AddNtBranch(dx_fp,'hsxpfp')
	call AddNtBranch(y_fp,'hsyfp')
	call AddNtBranch(dy_fp,'hsypfp')
	call AddNtBranch(px_lab,'MC_RSE_px')
	call AddNtBranch(py_lab,'MC_RSE_py')
	call AddNtBranch(pz_lab,'MC_RSE_pz')
	call AddNtBranch(hms_stop_id,'hms_stop_ID')
	
	Write(6,*) GetNumBranches(),' branches'
C   open root tree, and cycle over events
	n_trials = GetNtEntries()
	write(6,*) 'Number of generated events: ', n_trials
	do Itrial = 1,n_trials
	   if(ispec.eq.1) then
	      armSTOP_successes=hSTOP_successes
	   elseif(ispec.eq.2) then
	      armSTOP_successes=shmsSTOP_successes
	   endif
	  if(mod(Itrial,5000).eq.0) write(*,*)'event #: ',
     >Itrial,'       successes: ',armSTOP_successes

c get electron info from root tree	  
	  call ReadNtBranch(Itrial)
c       reset reconstructed components
	  px_lab=-1.0d10
	  py_lab=-1.0d10
	  pz_lab=-1.0d10
	  dpp_recon=-1.0d10
	  dph_recon=-1.0d10
	  dth_recon=-1.0d10
	  ytar_recon=-1.0d10

c convert to MeV	  
	  px_gen=px_gen*1000.0
	  py_gen=py_gen*1000.0
	  pz_gen=pz_gen*1000.0
	  z=zgen
	  
c generate x/y from instrimcic beam size	  
	  x = gauss1(th_nsig_max) * gen_lim(4) / 6.0 !beam width
	  y = gauss1(th_nsig_max) * gen_lim(5) / 6.0 !beam height

C Generate x/y from raster	 
C DJG Assume flat raster
	  fr1 = (grnd() - 0.5) * gen_lim(7)   !raster x
	  fr2 = (grnd() - 0.5) * gen_lim(8)   !raster y

	  fry = -fr2  !+y = up, but fry needs to be positive when pointing down

	  x = x + fr1
	  y = y + fr2
	  
	  x = x + xoff
	  y = y + yoff
	  z = z + zoff

c       convert to spectrometer coordinates
	  ptot = sqrt(px_gen**2+py_gen**2+pz_gen**2)
	  px_vert =-py_gen
	  py_vert = px_gen
	  pz_vert = pz_gen

	  ux=px_vert/ptot
	  uy=py_vert/ptot
	  uz=pz_vert/ptot

	  if(ispec.eq.1) then	!HMS
	     phi_spec=3.0*pi/2.0
	  elseif (ispec.eq.2) then
	     phi_spec=pi/2.0
	  endif
	  ux0=sin(th_spec)*cos(phi_spec)
	  uy0=sin(th_spec)*sin(phi_spec)
	  uz0=cos(th_spec)

	  cos_dtheta = ux*ux0 + uy*uy0 + uz*uz0
	  dxdz = ux/cos_dtheta
	  dydz = sqrt(1.0/cos_dtheta**2-1.0-dxdz**2)
	  y_event = uy/cos_dtheta
	  if (y_event .lt. uy0) dydz = -dydz
	  
	  dpp = (ptot/p_spec-1.0)*100.0

C Pick scattering angles and DPP from independent, uniform distributions.
C dxdz and dydz in HMS TRANSPORT coordinates.

c	  dpp  = grnd()*(gen_lim_up(1)-gen_lim_down(1))
c     &             + gen_lim_down(1)
c	  dydz = grnd()*(gen_lim_up(2)-gen_lim_down(2))
c     &          /1000.   + gen_lim_down(2)/1000.
c	  dxdz = grnd()*(gen_lim_up(3)-gen_lim_down(3))
c     &          /1000.   + gen_lim_down(3)/1000.

	  if(ispec.eq.2) then ! SHMS
C Transform from target to SHMS (TRANSPORT) coordinates.
C Version for a spectrometer on the left-hand side: (i.e. SHMS)
	     x_s    = -y
	     y_s    = x * cos_ts - z * sin_ts
	     z_s    = z * cos_ts + x * sin_ts
	  elseif(ispec.eq.1) then ! HMS
C Below assumes that HMS is on the right-hand side of the beam
C line (looking downstream).
	     x_s    = -y
	     y_s    = x * cos_ts + z * sin_ts
	     z_s    = z * cos_ts - x * sin_ts
	  else
	     write(6,*) 'unknown spectrometer: stopping'
	     stop
	  endif

C DJG Apply spectrometer offsets
C DJG If the spectrometer if too low (positive x offset) a particle
C DJG at "x=0" will appear in the spectrometer to be a little high
C DJG so we should subtract the offset

	  x_s = x_s - spec_xoff
	  y_s = y_s - spec_yoff
	  z_s = z_s - spec_zoff

	  dpp_s  = dpp
	  dxdz_s = dxdz
	  dydz_s = dydz

C DJG Apply spectrometer angle offsets
	  dxdz_s = dxdz_s - spec_xpoff/1000.0
	  dydz_s = dydz_s - spec_ypoff/1000.0

C Drift back to zs = 0, the plane through the target center
	  x_s = x_s - z_s * dxdz_s
	  y_s = y_s - z_s * dydz_s
	  z_s = 0.0

C Save init values for later.
	  xtar_init = x_s
	  ytar_init = y_s
	  ztar_init = z
	  dpp_init = dpp
	  dth_init = dydz_s*1000.		!mr
	  dph_init = dxdz_s*1000.		!mr

C Calculate multiple scattering length of target
	  if(ispec.eq.1) then ! spectrometer on right
	     cos_ev = (cos_ts+dydz_s*sin_ts)/sqrt(1+dydz_s**2+dxdz_s**2)
	  elseif(ispec.eq.2) then ! spectrometer on left
	     cos_ev = (cos_ts-dydz_s*sin_ts)/sqrt(1+dydz_s**2+dxdz_s**2)
	  endif
	  th_ev = acos(cos_ev)
	  sin_ev = sin(th_ev)

C Case 1 : extended cryo target:
C Choices: 
C 1. cryocylinder: Basic cylinder(2.65 inches diameter --> 3.37 cm radius) w/flat exit window (5 mil Al)
C 2. cryotarg2017: Cylinder (1.32 inches radisu)  with curved exit window (same radius) 5 mil sides/exit
C 3. Tuna can: shaped like a tuna can - 4 cm diameter (usually)  - 5 mil window. 
	  if (gen_lim(6).gt.3.) then ! anything longer than 3 cm assumed to be cryotarget
c	    call cryotuna(z,th_ev,rad_len_cm,gen_lim(6),musc_targ_len)
c	    call cryocylinder(z,th_ev,rad_len_cm,gen_lim(6),musc_targ_len)
	     call cryotarg2017(z,th_ev,rad_len_cm,gen_lim(6),musc_targ_len)
C Simple solid target
	 else
	    if (gen_lim(6).gt.0) then
	    musc_targ_len = abs(gen_lim(6)/2. - z)/rad_len_cm/cos_ev
	    else ! using multifoil target
	    musc_targ_len = abs(foil_tk/2. - (z-foil_zcent))/rad_len_cm/cos_ev	       
	    endif
	 endif

C Scattering before magnets:  Approximate all scattering as occuring AT TARGET.
C SHMS
C  20 mil Al scattering chamber window (X0=8.89cm)
C  57.27 cm air (X0=30420cm)
C spectrometer entrance window
C  10 mil Al s (X0=8.89cm)

C  HMS
C  20 mil Al scattering chamber window (X0=8.89cm)
C  24.61 cm air (X0=30420cm)
C spectrometer entrance window
C  15 mil Kevlar (X0=74.6 cm)
C   5 mil Mylar (X0=28.7 cm)

c       b's for RC
	 bproton = 1.35742
	 bair=1.36472
	 bal=1.36676
	 bkevlar=1.3642
	 bmylar=1.3632

c	 bt_tot= 0.0
	 bt_tot=bproton*musc_targ_len
	  if(ispec.eq.2) then
	     musc_targ_len = musc_targ_len + .020*2.54/8.89 +
     >          57.27/30420. +  .010*2.54/8.89
	  elseif(ispec.eq.1) then
	     Lal=.020*2.54/8.89 ! rad length of exit window
	     Lair=24.61/30420.0
	     Lkevlar=.015*2.54/74.6
	     Lmylar= .005*2.54/28.7
	     musc_targ_len = musc_targ_len + Lal + Lair+Lkevlar + Lmylar
	     bt_tot = bt_tot + bal*Lal + bair*Lair + bkevlar*Lkevlar + bmylar*Lmylar
c       generate radiated photon energy for scattered electron (external)
	     ranphot=grnd()
	     deltaE = ptot*ranphot**(1.0/bt_tot)
	     ptot = ptot-deltaE
	     dpp = (ptot/p_spec-1.0)*100.0
	     dpp_s  = dpp
	  endif

c
	  if (ms_flag ) call musc(m2,p_spec*(1.+dpp_s/100.),
     > musc_targ_len,dydz_s,dxdz_s)

! ----------------------------------------------------------------------------
	  if(ispec.eq.2) then

	     call mc_shms(p_spec, th_spec, dpp_s, x_s, y_s, z_s, 
     >          dxdz_s, dydz_s,
     >		x_fp, dx_fp, y_fp, dy_fp, m2, shms_spec,
     >		ms_flag, wcs_flag, decay_flag, resmult, xtar_init, ok_spec, 
     >          pathlen, 5)

c	     if (spec_ntuple) then
c		shms_spec(59) = shmsSTOP_id
c		shms_spec(58)= ytar_init
c		shms_spec(53)= dpp_init
c		shms_spec(54)= dph_init ! dx/dz (mr)
c		shms_spec(55)= dth_init ! dy/dz (mr)
cc            if (ok_spec) spec(58) =1.
c		do ivar=1,SpecNtupleSize
c		   write(SpecNtupleIO) shms_spec(ivar)
c		enddo
c	     endif
	  elseif(ispec.eq.1) then
	     call mc_hms(p_spec, th_spec, dpp_s, x_s, y_s, z_s, 
     >          dxdz_s, dydz_s,
     >          x_fp, dx_fp, y_fp, dy_fp, m2,
     >          ms_flag, wcs_flag, decay_flag, resmult, xtar_init, ok_spec, 
     >          pathlen)
	  else
	     write(6,*) 'Unknown spectrometer! Stopping..'
	     stop
	  endif

	  if (ok_spec) then !Success, increment arrays
	    dpp_recon = dpp_s
            dth_recon = dydz_s			!mr
	    dph_recon = dxdz_s
	    ztar_recon = + y_s / sin_ts 
            ytar_recon = y_s

C Compute sums for calculating reconstruction variances.
	    dpp_var(1) = dpp_var(1) + (dpp_recon - dpp_init)
	    dth_var(1) = dth_var(1) + (dth_recon*1000.0 - dth_init)
	    dph_var(1) = dph_var(1) + (dph_recon*1000.0 - dph_init)
	    ztg_var(1) = ztg_var(1) + (ztar_recon - ztar_init)

	    dpp_var(2) = dpp_var(2) + (dpp_recon - dpp_init)**2
	    dth_var(2) = dth_var(2) + (dth_recon*1000.0 - dth_init)**2
	    dph_var(2) = dph_var(2) + (dph_recon*1000.0 - dph_init)**2
	    ztg_var(2) = ztg_var(2) + (ztar_recon - ztar_init)**2

C       calculate the reconstructed quantities (this is the same algorithm as hcana)
	    pe_recon = p_spec*(1.0+dpp_recon/100.0)
	    pe_recon = pe_recon/1000.0 ! convert to GeV
	    pz_spec = pe_recon/sqrt(1.0+dydz_s**2+dxdz_s**2)
	    px_spec = pz_spec*dxdz_s
	    py_spec = pz_spec*dydz_s

	    px_lab = pz_spec*(dydz_s*cos_ts - sin_ts)
	    py_lab = -px_spec
	    pz_lab = pz_spec*(dydz_s*sin_ts + cos_ts)

	 endif			!Incremented the arrays

	 hms_stop_id = hSTOP_ID
	 call FillNTBranch('hsdelta')
	 call FillNTBranch('hsxptar')
	 call FillNTBranch('hsyptar')
	 call FillNTBranch('hsytar')
	 call FillNTBranch('hsxfp')
	 call FillNTBranch('hsxpfp')
	 call FillNTBranch('hsyfp')
	 call FillNTBranch('hsypfp')
	 call FillNTBranch('MC_RSE_px')
	 call FillNTBranch('MC_RSE_py')
	 call FillNTBranch('MC_RSE_pz')
	 call FillNTBranch('hms_stop_ID')



C Output NTUPLE entry.
C This is ugly, but want the option to have different outputs
C for spectrometer ntuples
c	 if(ispec.eq.2) then
c	    if (store_all.OR.(hut_ntuple.AND.ok_spec)) then
c	       shms_hut(1) = x_fp
c	       shms_hut(2) = y_fp
c	       shms_hut(3) = dx_fp
c	       shms_hut(4) = dy_fp
c	       shms_hut(5) = ztar_init
c	       shms_hut(6) = ytar_init
c	       shms_hut(7) = dpp_init
c	       shms_hut(8) = dth_init/1000.
c	       shms_hut(9) = dph_init/1000.
c	       shms_hut(10) = ztar_recon
c	       shms_hut(11) = ytar_recon
c	       shms_hut(12)= dpp_recon
c	       shms_hut(13)= dth_recon
c	       shms_hut(14)= dph_recon
c	       shms_hut(15)= xtar_init
c	       shms_hut(16)= fry
c	       if (use_front_sieve) then
c		  shms_hut(17)= xsfr_num
c		  shms_hut(18)= ysfr_num
c		  shms_hut(19)= xc_frsieve
c		  shms_hut(20)= yc_frsieve
c	       endif
c	       if (use_sieve) then
c		  shms_hut(17)= xs_num
c		  shms_hut(18)= ys_num
c		  shms_hut(19)= xc_sieve
c		  shms_hut(20)= yc_sieve
c	       endif
c	       shms_hut(21)= shmsSTOP_id
c	       shms_hut(22)= x
c	       shms_hut(23)= y
c	       do ivar=1,NtupleSize
c		  write(NtupleIO) shms_hut(ivar)
c	       enddo
c	    endif
c	 endif
c
c	 if(ispec.eq.1) then
c	    if (store_all.OR.(hut_ntuple.AND.ok_spec)) then
c	       hms_hut(1) = x_fp
c	       hms_hut(2) = y_fp
c	       hms_hut(3) = dx_fp
c	       hms_hut(4) = dy_fp
c	       hms_hut(5) = xtar_init
c	       hms_hut(6) = ytar_init
c	       hms_hut(7) = dph_init/1000.
c	       hms_hut(8) = dth_init/1000.
c	       hms_hut(9) = ztar_init
c	       hms_hut(10)= dpp_init
c	       hms_hut(11)= ytar_recon
c	       hms_hut(12)= dph_recon
c	       hms_hut(13)= dth_recon
c	       hms_hut(14)= ztar_recon
c	       hms_hut(15)= dpp_recon
c	       hms_hut(16)= fry
c	       if (use_sieve) then
c		  hms_hut(17)= xs_num
c		  hms_hut(18)= ys_num
c		  hms_hut(19)= xc_sieve
c		  hms_hut(20)= yc_sieve
c	       endif
c              hms_hut(21)=hSTOP_id
c	       hms_hut(22)= x
c	       hms_hut(23)= y
c	       do ivar=1,NtupleSize
c		  write(NtupleIO) hms_hut(ivar)
c	       enddo
c	    endif
c	 endif


	 
C We are done with this event, whether GOOD or BAD.
C Loop for remainder of trials.

500	  continue

	  
	enddo			!End of M.C. loop
	call PrintNT()
	call RootNTOutp()      !close root file

C------------------------------------------------------------------------------C
C                           End of Monte-Carlo loop                            C
C------------------------------------------------------------------------------C

C Close NTUPLE file.
c
c	close(NtupleIO)
c	if (spec_ntuple) close(SPecNtupleIO)


	write (chanout,1002)
	write (chanout,1003) p_spec,th_spec*degrad
        write (chanout,1004) (gen_lim(i),i=1,6)

	write (chanout,1005) n_trials

	if(ispec.eq.1) then
	   armSTOP_successes=hSTOP_successes
	   armSTOP_trials=hSTOP_trials
	elseif(ispec.eq.2) then
	   armSTOP_successes=shmsSTOP_successes
	   armSTOP_trials=shmsSTOP_trials
	endif

C Indicate where particles are lost in spectrometer.
	if(ispec.eq.2) then
	   write (chanout,1015)
     >        shmsSTOP_targ_hor,shmsSTOP_targ_vert,shmsSTOP_targ_oct,
     >        shmsSTOP_FRONTSLIT_hor,shmsSTOP_FRONTSLIT_vert,
     >        shmsSTOP_HB_in,shmsSTOP_HB_men,shmsSTOP_HB_mex,
     >        shmsSTOP_HB_out,shmsSTOP_DOWNSLIT,
     >        shmsSTOP_COLL_hor,shmsSTOP_COLL_vert,shmsSTOP_COLL_oct,
     >        shmsSTOP_Q1_in,shmsSTOP_Q1_men,
     >        shmsSTOP_Q1_mid,shmsSTOP_Q1_mex,shmsSTOP_Q1_out,
     >        shmsSTOP_Q2_in,shmsSTOP_Q2_men,shmsSTOP_Q2_mid,
     >        shmsSTOP_Q2_mex,shmsSTOP_Q2_out,
     >        shmsSTOP_Q3_in,shmsSTOP_Q3_men,shmsSTOP_Q3_mid,
     >        shmsSTOP_Q3_mex,shmsSTOP_Q3_out,
     >        shmsSTOP_D1_in,shmsSTOP_D1_flr,shmsSTOP_D1_men,
     >        shmsSTOP_D1_mid1,shmsSTOP_D1_mid2,shmsSTOP_D1_mid3,
     >        shmsSTOP_D1_mid4,shmsSTOP_D1_mid5,shmsSTOP_D1_mid6,
     >        shmsSTOP_D1_mid7,shmsSTOP_D1_mex,shmsSTOP_D1_out,
     >        shmsSTOP_BP_in, shmsSTOP_BP_out

	   write (chanout,1006)
     >	   shmsSTOP_trials,shmsSTOP_hut,shmsSTOP_dc1,shmsSTOP_dc2,
     >     shmsSTOP_s1,shmsSTOP_s2,shmsSTOP_s3,shmsSTOP_cal,
     >     shmsSTOP_successes,shmsSTOP_successes

	elseif(ispec.eq.1) then
	   write (chanout,1016)
     >        hSTOP_fAper_hor,hSTOP_fAper_vert,hSTOP_fAper_oct,
     >        hSTOP_bAper_hor,hSTOP_bAper_vert,hSTOP_bAper_oct,
     >        hSTOP_slit,
     >        hSTOP_Q1_in,hSTOP_Q1_mid,hSTOP_Q1_out,
     >        hSTOP_Q2_in,hSTOP_Q2_mid,hSTOP_Q2_out,
     >        hSTOP_Q3_in,hSTOP_Q3_mid,hSTOP_Q3_out,
     >        hSTOP_D1_in,hSTOP_D1_out,
     >        hSTOP_VP1,hSTOP_VP2,hSTOP_VP3,hSTOP_VP4

	   write (chanout,1007)
     >	   hSTOP_trials,hSTOP_hut,hSTOP_dc1,hSTOP_dc2,hSTOP_scin,hSTOP_cal,
     >     hSTOP_successes,hSTOP_successes
	endif


C Compute reconstruction resolutions.

	if (armSTOP_successes.eq.0) armSTOP_successes=1
	t1 = sqrt(max(0.,dpp_var(2)/armSTOP_successes 
     > - (dpp_var(1)/armSTOP_successes)**2))
	t2 = sqrt(max(0.,dth_var(2)/armSTOP_successes 
     > - (dth_var(1)/armSTOP_successes)**2))
	t3 = sqrt(max(0.,dph_var(2)/armSTOP_successes 
     > - (dph_var(1)/armSTOP_successes)**2))
	t4 = sqrt(max(0.,ztg_var(2)/armSTOP_successes 
     > - (ztg_var(1)/armSTOP_successes)**2))

	write (chanout,1011) dpp_var(1)/armSTOP_successes,t1,
     > dth_var(1)/armSTOP_successes,
     >		t2,dph_var(1)/armSTOP_successes,t3,
     > ztg_var(1)/armSTOP_successes,t4

	write(6,*) armSTOP_trials,' Trials',armSTOP_successes
     > ,' Successes'
	write (6,1011) dpp_var(1)/armSTOP_successes,t1,
     > dth_var(1)/armSTOP_successes,
     >		t2,dph_var(1)/armSTOP_successes,t3,
     > ztg_var(1)/armSTOP_successes,t4

C ALL done!

	stop ' '

C =============================== Format Statements ============================

1001	format(a)
1002	format('!',/,'! Uniform illumination Monte-Carlo results')
1003	format('!',/'! Spectrometer setting:',/,'!',/,
     >g11.5,' =  P  spect (MeV)',/,
     >g11.5,' =  TH spect (deg)')

1004	format('!',/'! Monte-Carlo limits:',/,'!',/,
     >  g11.5,'= GEN_LIM(1) - DP/P   (half width,% )',/,
     >  g11.5,'= GEN_LIM(2) - Theta  (half width,mr)',/,
     >  g11.5,'= GEN_LIM(3) - Phi    (half width,mr)',/,
     >  g11.5,'= GEN_LIM(4) - HORIZ (full width of 3 sigma cutoff,cm)',/,
     >  g11.5,'= GEN_LIM(5) - VERT  (full width of 3 sigma cutoff,cm)',/,
     >  g11.5,'= GEN_LIM(6) - Z      (Full width,cm)')

!inp     >	,/,
!inp     >	g18.8,' =  Hor. 1/2 gap size (cm)',/,
!inp     >	g18.8,' =  Vert. 1/2 gap size (cm)')

1005	format('!',/,'! Summary:',/,'!',/,
!     >	i,' Monte-Carlo trials:')
     >  i11,' Monte-Carlo trials:')

1006	format(i11,' Initial Trials',/
     >  i11,' Trials made it to the hut',/
     >  i11,' Trial cut in dc1',/
     >  i11,' Trial cut in dc2',/
     >  i11,' Trial cut in s1',/
     >  i11,' Trial cut in s2',/
     >  i11,' Trial cut in s3',/
     >  i11,' Trial cut in cal',/
     >  i11,' Trials made it thru the detectors and were reconstructed',/
     >  i11,' Trials passed all cuts and were histogrammed.',/
     >  )

 1007	format(i11,' Initial Trials',/
     >  i11,' Trials made it to the hut',/
     >  i11,' Trial cut in dc1',/
     >  i11,' Trial cut in dc2',/
     >  i11,' Trial cut in scin',/
     >  i11,' Trial cut in cal',/
     >  i11,' Trials made it thru the detectors and were reconstructed',/
     >  i11,' Trials passed all cuts and were histogrammed.',/
     >  )


!1008	format(8i)
!1009	format(1x,i4,g,i)
!1010	format(a,i)
1011	format(
     >  'DPP ave error, resolution = ',2g18.8,' %',/,
     >  'DTH ave error, resolution = ',2g18.8,' mr',/,
     >  'DPH ave error, resolution = ',2g18.8,' mr',/,
     >  'ZTG ave error, resolution = ',2g18.8,' cm')

1012	format(1x,16i4)

1015	format(/,
     >     i11,' stopped in the TARG APERT HOR',/
     >     i11,' stopped in the TARG APERT VERT',/
     >     i11,' stopped in the TARG APERT OCTAGON',/
     >     i11,' stopped in the FIXED FRONT SLIT HOR (id=-1)',/
     >     i11,' stopped in the FIXED FRONT SLIT VERT (id=-1)',/
     >     i11,' stopped in HB ENTRANCE (id=1)',/
     >     i11,' stopped in HB MAG ENTRANCE (id=2)',/
     >     i11,' stopped in HB MAG EXIT (id=3)',/
     >     i11,' stopped in HB EXIT (id=4)',/
     >     i11,' stopped in the DOWN SIEVE SLIT (id=99)',/
     >     i11,' stopped in the COLLIMATOR HOR (id=5)',/
     >     i11,' stopped in the COLLIMATOR VERT (id=5)',/
     >     i11,' stopped in the COLLIMATOR OCTAGON (id=5)',/
     >     i11,' stopped in Q1 ENTRANCE (id=6)',/
     >     i11,' stopped in Q1 MAG ENTRANCE (id=7)',/
     >     i11,' stopped in Q1 MIDPLANE (id=8)',/
     >     i11,' stopped in Q1 MAG EXIT (id=9)',/
     >     i11,' stopped in Q1 EXIT (id=10)',/
     >     i11,' stopped in Q2 ENTRANCE (id=11)',/
     >     i11,' stopped in Q2 MAG ENTRANCE (id=12)',/
     >     i11,' stopped in Q2 MIDPLANE (id=13)',/
     >     i11,' stopped in Q2 MAG EXIT (id=14)',/
     >     i11,' stopped in Q2 EXIT (id=15)',/
     >     i11,' stopped in Q3 ENTRANCE (id=16)',/
     >     i11,' stopped in Q3 MAG ENTRANCE (id=17)',/
     >     i11,' stopped in Q3 MIDPLANE (id=18)',/
     >     i11,' stopped in Q3 MAG EXIT (id=19)',/
     >     i11,' stopped in Q3 EXIT (id=20)',/
     >     i11,' stopped in D1 ENTRANCE (id=21)',/
     >     i11,' stopped in D1 FLARE (id=22)',/
     >     i11,' stopped in D1 MAG ENTRANCE (id=23)',/
     >     i11,' stopped in D1 MID-1 (id=24)',/
     >     i11,' stopped in D1 MID-2 (id=25)',/
     >     i11,' stopped in D1 MID-3 (id=26)',/
     >     i11,' stopped in D1 MID-4 (id=27)',/
     >     i11,' stopped in D1 MID-5 (id=28)',/
     >     i11,' stopped in D1 MID-6 (id=29)',/
     >     i11,' stopped in D1 MID-7 (id=30)',/
     >     i11,' stopped in D1 MAG EXIT (id=31)',/
     >     i11,' stopped in D1 EXIT (id=32)',/
     >     i11,' stopped in BP ENTRANCE',/
     >     i11,' stopped in BP EXIT',/
     >     )

 1016 format(/,
     >     i11,' stopped in the Front-end Aperture HOR (id=5)',/
     >     i11,' stopped in the Front-end Aperture VERT (id=5)',/
     >     i11,' stopped in the Front-end Aperture OCTAGON (id=5)',/
     >     i11,' stopped in the Back-end Aperture HOR (id=6)',/
     >     i11,' stopped in the Back-end Aperture VERT (id=6)',/
     >     i11,' stopped in the Back-end Aperture OCTAGON (id=6)',/
     >     i11,' stopped in Sieve Slit (id=99)',/
     >     i11,' stopped in Q1 MAG ENTRANCE (id=7)',/
     >     i11,' stopped in Q1 MIDPLANE (id=8)',/
     >     i11,' stopped in Q1 MAG EXIT (id=9)',/
     >     i11,' stopped in Q2 MAG ENTRANCE (id=12)',/
     >     i11,' stopped in Q2 MIDPLANE (id=13)',/
     >     i11,' stopped in Q2 MAG EXIT (id=14)',/
     >     i11,' stopped in Q3 MAG ENTRANCE (id=17)',/
     >     i11,' stopped in Q3 MIDPLANE (id=18)',/
     >     i11,' stopped in Q3 MAG EXIT (id=19)',/
     >     i11,' stopped in D1 ENTRANCE (id=23)',/
     >     i11,' stopped in D1 EXIT (id=31)',/
     >     i11,' stopped in Vacuum Pipe Plane-1 (id=32)',/
     >     i11,' stopped in Vacuum Pipe Plane-2 (id=33)',/
     >     i11,' stopped in Vacuum Pipe Plane-3 (id=34)',/
     >     i11,' stopped in Vacuum Pipe Plane-4 (id=35)',/
     >     )

1100	format('!',79('-'),/,'! ',a,/,'!')
1200	format(/,'! ',a,' Coefficients',/,/,
     >  (5(g18.8,','))
     >  )
1300	format(/,'! ',a,' Coefficient uncertainties',/,/,
     >  (5(g18.8,','))
     >  )

	end
