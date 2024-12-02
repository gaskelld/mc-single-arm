      subroutine hms_hbook_init(filename,spec_ntuple)
C Initialize output names differently based on which spectrometer we
C are using
C
      implicit none
      include 'hbook.inc'
      character*80 filename
      logical spec_ntuple
      character*16	hut_nt_names(10)/
     >    'eventnumber','hsxfp', 'hsyfp', 'hsxpfp', 'hsypfp',
     >     'hsytar','hsxptar','hsyptar','hsdelta',
     >     'hms_stop_id'/

      integer*4 i

      NtupleIO=30
      NtupleSize=10

      open(NtupleIO,file=filename,form="unformatted",access="sequential")

      write(NtupleIO) NtupleSize
      do i=1,NtupleSize
         write(NtupleIO) hut_nt_names(i)
      enddo

      return
      end
      
