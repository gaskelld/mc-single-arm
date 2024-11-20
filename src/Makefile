RM        = rm -f 
SHELL     = /bin/sh

H      = hms/
SH     = shms/
SHARE  = shared/
CH     = cern/

my_objs = mc_single_arm.o hms_hbook_init.o shms_hbook_init.o \
          target_cans.o \
          $(H)mc_hms.o $(H)mc_hms_hut.o $(H)mc_hms_recon.o \
          $(SH)mc_shms.o $(SH)mc_shms_hut.o $(SH)mc_shms_recon.o \
	  $(SHARE)project.o $(SHARE)transp.o $(SHARE)rotate_haxis.o \
	  $(SHARE)rotate_vaxis.o $(SHARE)locforunt.o $(SHARE)musc.o \
          $(SHARE)musc_ext.o $(SHARE)stringlib.o mt19937.o gauss1.o \
	  $(CH)lfit.o $(CH)ranlux.o $(CH)fint.o $(CH)kerset.o $(CH)abend.o \
          loren.o

MYOS := $(subst -,,$(shell uname))
FFLAGSA=-O -ffixed-line-length-132 -ff2c -fno-automatic -fdefault-real-8
FFLAGS= $(INCLUDES) $(FFLAGSA)
OTHERLIBS = -L/usr/lib
FC  := gfortran
F77 :=gfortran

ROOTSYS=/cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/root/6.30.04-gcc11.4.0
#ROOTSYS=/home/gaskelld/root/
CXX=g++
CXXFLAGS0 = -O3 -std=c++17
CXXFLAGS=$(CXXFLAGS0) $(DROOT)
ROOTDIR=$(ROOTSYS)
DMYROOT= -DMYROOT
ROOTLIBS     := $(shell $(ROOTDIR)/bin/root-config --prefix=$(ROOTDIR)  --libs)
ROOTINCLUDE  := -I $(shell $(ROOTDIR)/bin/root-config --prefix=$(ROOTDIR) --incdir)

all: mc_single_arm 

none: mc_single_arm

mc_single_arm:  $(my_objs) froot.co Makefile
	$(F77) -o $@ $(FFLAGS) $(my_objs) froot.co $(ROOTLIBS) $(OTHERLIBS) -L/usr/lib64 -lstdc++

#%.o: %.f
#	$(F77) $(FFLAGS) -c $<

%.co: %.c
	$(CXX) -c $(CXXFLAGS) $(DMYROOT) $(ROOTINCLUDE) -o $@ $<

clean:
	rm -f *.o $(H)*.o $(SH)*.o $(SHARE)*.o $(CH)*.o *.co mc_single_arm 

