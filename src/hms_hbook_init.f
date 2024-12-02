      subroutine hms_hbook_init(filename,spec_ntuple)
C Initialize output names differently based on which spectrometer we
C are using
C
      implicit none
      include 'hbook.inc'
      character*80 filename
      logical spec_ntuple
      character*16	hut_nt_names(13)/
     >    'eventnumber','hsxfp', 'hsyfp', 'hsxpfp', 'hsypfp',
     >     'hsytar','hsxptar','hsyptar','hsdelta',
     >     'hms_stop_id', 'MC_RSE_px', 'MC_RSE_py', 'MC_RSE_pz'/

      integer*4 i

      NtupleIO=30
      NtupleSize=13

      open(NtupleIO,file=filename,form="unformatted",access="sequential")

      write(NtupleIO) NtupleSize
      do i=1,NtupleSize
         write(NtupleIO) hut_nt_names(i)
      enddo

      return
      end
      
