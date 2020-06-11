# This is the main makefile.
# Use with GNU make.
# Relies on makedepf90 to do dependency lists:
# You can obtain this free program from:
# http://www.helsinki.fi/~eedelman/makedepf90.html
# New $ARCH entries may be required, depending on your compiler.
# Also requires fpx3 fortran preprocessor available at
# http://wwwuser.gwdg.de/~jbehren/fpx3.html

SHELL = /bin/sh
EXE = non_periodic_lagrange

LIBPATH = -L/glade/u/home/mostafa/softwares/fftw-2.1.5/lib
LIBS = $(LIBPATH) -lrfftw -lfftw -lm -limf

USE_MPI = yes
USE_LVLSET = no
FPP = ./fpx3
ifeq ($(USE_MPI), yes)
  FPP += -DMPI
   #LIBS += -lmpichf90 -lfmpich -lmpich
  # LIBS +=  -lmpi_f90 -lmpi_f77 -lmpi
endif
ifeq ($(USE_LVLSET), yes)
  FPP += -DLVLSET
endif

# Directory for the .o files
OPATH = obj
# Directory for the .mod files, if your compiler generates them
# May want to just make this 'obj' as well
MPATH = mod
FCOMP = ifort

ifeq ($(FCOMP),ifort)
  
  FPP += -DIFORT
  
  ifeq ($(USE_MPI), yes)
    FC = mpif90
  else
    FC = ifort
  endif

  #FFLAGS = -O0 -check bounds -g -debug all -traceback
  #FFLAGS = -fast
  #FFLAGS = -O3 -ipo
  FFLAGS = -O2
  #FFLAGS = -O3 -ip  -ftz
  #FFLAGS = -axSSE4.2 -xS -ftz -ip -ipo -O3 
  #FFLAGS += -warn all
  #FDEBUG = -g -debug all
  FPROF = -p
  #LDFLAGS = -threads -shared-intel
  MODDIR = -I$(MPATH) -module $(MPATH)
  FFLAGS += $(MODDIR)
endif

ifeq ($(FCOMP),gfortran)
  FPP += -DGFORTRAN
  ifeq ($(USE_MPI), yes)
    FC = mpif90
  else
    FC = gfortran
  endif
  FFLAGS = -O0
#  FFLAGS = -O2 -ffree-form -ffixed-line-length-none
#  FFLAGS += -Wall
  FDEBUG = -g
  FPROF = -p
#  LDFLAGS = -static 
  MODDIR = -I$(MPATH) -J$(MPATH)  
  FFLAGS += $(MODDIR)  
endif

SRCS =  \
        bottombc.f90 \
        convec.f90 \
        calc_tke.f90 \
        calc_spectra.f90 \
        ddx.f90 ddxy.f90 ddy.f90 ddz_uv.f90 ddz_w.f90 \
        dealias1.f90 dealias2.f90 debug_mod.f90 \
        divstress_uv.f90 divstress_w.f90 dns_stress.f90 \
        energy.f90 \
        fft.f90 filt_da.f90 forcing.f90 \
        ic.f90 ic_dns.f90 immersedbc.f90 initial.f90 \
        interpolag_Sdep.f90 interpolag_Ssim.f90 io.f90 \
        interpolag_scalar_Sdep.f90 \
        lagrange_Sdep.f90 lagrange_Ssim.f90 \
        lagrange_scalar_Sdep.f90 \
	lagrange_scalar_newSdep.f90 \
        main.f90 messages.f90 \
        new_scaledep_dynamic.f90 \
        padd.f90 param.f90 press_stag_array.f90 \
        ran3.f90 rmsdiv.f90 \
        scaledep_dynamic.f90 scalars_module.f90 scalars_module2.f90 \
        scalar_scaledep_dynamic.f90 \
        sgs_stag.f90 sgsmodule.f90 sim_param.f90 \
        spectra.f90 \
        std_dynamic.f90 \
        test_filtermodule.f90 tecplot.f90 topbc.f90 \
        tridag.f90 tridag_array.f90 types.f90 \
        unpadd.f90 \
        wallstress.f90 wallstress_dns.f90 \
        write_spectra.f90 \

LVLSET_SRCS = level_set_base.f90 level_set.f90 linear_simple.f90

ifeq ($(USE_MPI), yes)
  SRCS += tridag_array_pipelined.f90
endif
ifeq ($(USE_LVLSET), yes)
  SRCS += $(LVLSET_SRCS)
endif

COMPSTR = '$(FPP) $$< > t.$$<; $$(FC) -c -o $$@ $$(FFLAGS) t.$$<; rm -f t.$$<'

include .depend

.depend: $(SRCS)
	mkdir -p $(OPATH) $(MPATH);
	/glade/u/home/mostafa/softwares/makedepf90-2.7.0/makedepf90 -r $(COMPSTR) -d $(OPATH) -o $(EXE) $(SRCS) > .depend
	# /tigress-hsm/jingh/bin/makedepf90

debug:
	$(MAKE) $(EXE) "FFLAGS = $(FDEBUG) $(FFLAGS)"

prof:
	$(MAKE) $(EXE) "FFLAGS = $(FPROF) $(FFLAGS)"

# Other support programs are listed below this point
#interp: interp.f90
interp: interp_scalars.f90
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $<
#interp: interp_scalars.f90
#	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $<


# This doesn't remove .mod files--should be OK as long a dependency list 
# for the .o files is correct.
# FOBJ is defined in .depend
.PHONY : clean
clean :
	#echo \0>./total_time.dat
	rm -rf $(FOBJ) .depend* $(MPATH)/*.mod
	#rm $(EXE)
	#mkdir ./output/fields_3d
	#mkdir output/code_used
	#cp $(SRCS) param.f90 makefi:le ./output/code_used/
