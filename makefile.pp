#
#echo "Make sure MPICH/GM params set in both global.h and makefile!"
#echo "Make sure MPICH/GM params set in both global.h and makefile!"
USEMPI=0

USECCC=1
USEEGCS=1
USEGCC=1

ifeq ($(USEGCC),1)
COMP = gcc
endif
ifeq ($(USEEGCS),1)
COMP = egcs
endif
# ccc overrides above
ifeq ($(USECCC),1)
COMP = ccc
endif

ifeq ($(USEMPI),1)
MCC=mpicc
endif

#General env. vars
RM = rm -f
SUFF = c
#
# Define a cpp source file directory, a fortran file subdirectory,
# an object file subdirectory and an executable file subdirectory
#
SRCD = .
OBJD = .
OBJD2 = ./Obj
BIND = ./bin

CMD = $(BIND)/postproc

PREP = prep
FINISH = finish

#
# Define preprocessor and compile flags, and linked libraries
#
# linux on a pentium
ifeq ($(HOSTTYPE),i386)
# override since no ccc
LOCALLINUX=0
ifeq ($(HOSTNAME),metric.physics.uiuc.edu)
LOCALLINUX=1
endif
ifeq ($(HOSTNAME),kerr.physics.uiuc.edu)
LOCALLINUX=1
endif


# override for now
#LOCALLINUX=1


ifeq ($(LOCALLINUX),1)
ifeq ($(USEGCC),1)
COMP = gcc
endif
ifeq ($(USEEGCS),1)
COMP = egcs
endif
#CFLAGS = -Wall -mpentium -O3 -pipe  -malign-loops=2 -malign-jumps=2 -malign-functions=2 -DCPU=686 -DNEED_GETOPT -DLINUX -ffast-math
#CFLAGS = -Wall -mpentium -O3 -pipe  -malign-loops=2 -malign-jumps=2 -malign-functions=2 -DCPU=686 -DNEED_GETOPT -DLINUX -ffast-math -pg
#CFLAGS = -Wall -mcpu=pentiumpro -march=pentiumpro -O4 -pipe  -malign-loops=2 -malign-jumps=2 -malign-functions=2 -malign-double -mstack-align-double -ffast-math -pg
#CFLAGS = -Wall -mcpu=pentiumpro -march=pentiumpro -O4 -malign-loops=2 -malign-jumps=2 -malign-functions=2 -malign-double -mstack-align-double -ffast-math -finline-functions -pg
#CFLAGS = -Wall -mcpu=pentiumpro -march=pentiumpro -O4 -malign-loops=2 -malign-jumps=2 -malign-functions=2 -malign-double -ffast-math -finline-functions -pg
CFLAGS = -Wall -mcpu=pentiumpro -march=pentiumpro -O3 -malign-loops=2 -malign-jumps=2 -malign-functions=2 -ffast-math -finline-functions -pg -g -a
#CFLAGS = -Wall -mcpu=pentiumpro -march=pentiumpro -O4 -malign-loops=2 -malign-jumps=2 -malign-functions=2 -ffast-math -finline-functions -g
#-pg
#-pg -g  source lines
#-pg -g -a   line count
# gprof -l <file> > out.txt
# gprof -A -I<sourcedir>
# gprof -l -A -x s

#below is typical flags for double precision...can take -pg off for no profile
#add -mstack-align-double if using pgcc
#CFLAGS = -Wall -O0 -g
#  -fomit-frame-pointer



#CFLAGS = -Wall -O0
#CFLAGS = -O6 -g
#CFLAGS = -O0 -pg -g
LDFLAGS = -lm

#CC = cc
#AR	=	ar r
#RANLIB	=	ranlib
endif

ifeq ($(LOCALLINUX),0)
# testing linux cluster
# after compile, to submit do:
# submitjob 8 2 hh:mm:ss myrinet myjob myjob.out "/u/ncsa/user/mycommand arg1 arg2" 
COMP = pgCC
CFLAGS = -fast -Minline=levels:10 --no_exceptions
LDFLAGS = -lm
# for mpi:
ifeq ($(USEMPI),1)
CFLAGS =  -I/usr/local/mpich/include -fast -Minline=levels:10 --no_exceptions
LDFLAGS = -L/usr/local/mpich/lib -lmpich -lvmi -ldl -lpthread  -lm
endif

endif

CC=$(COMP)


endif # endif i386



# linux on a pentium686(pgcc)
ifeq ($(HOSTTYPE),i686)
# override since no ccc
ifeq ($(USEGCC),1)
COMP = gcc
endif
ifeq ($(USEEGCS),1)
COMP = egcs
endif

CFLAGS = -Wall -mcpu=pentiumpro -march=pentiumpro -O4 -malign-loops=2 -malign-jumps=2 -malign-functions=2 -malign-double -ffast-math -finline-functions -pg

#below is typical flags for double precision...can take -pg off for no profile
#add -mstack-align-double if using pgcc
#CFLAGS = -Wall -O0 -g
#  -fomit-frame-pointer

LDFLAGS = -lm

#CC = cc
#AR	=	ar r
#RANLIB	=	ranlib

CC=$(COMP)

endif # endif i686





# Define preprocessor and compile flags, and linked libraries

ifeq ($(HOSTTYPE),alpha)
LDFLAGS = -lm

ifeq ($(USECCC),1)
CDEBUG = -g3 # -g3 for higher opts than -O0
#CDEBUG = -g
# do profile
#CDEBUG = -pg
# production level
CFLAGS3 = -Wall -O4 -fast -msg_disable badsubscript -msg_disable subscrbounds -msg_disable unreachcode -msg_disable noparmlist -msg_disable subscrbounds2
# super annoying develop level
#CFLAGS3 = -Wall -O2 -fast
#CFLAGS3 = -fast -arch ev67
# debug level
#CFLAGS3 = -Wall -O0 -msg_disable badsubscript -msg_disable subscrbounds -msg_disable unreachcode -msg_disable noparmlist -msg_disable subscrbounds2
#CFLAGS3 = -Wall -O0
#CFLAGS3 = -Wall -O2

ifeq ($(USEMPI),1)

CC=$(MCC) -cc=$(COMP) $(CFLAGS3) $(CFLAGS2) $(CDEBUG)

ifeq ($(HOSTNAME),photon.physics.uiuc.edu)
CFLAGS2 = -arch ev6
endif

ifeq ($(HOSTNAME),rainman.physics.uiuc.edu)
CFLAGS2 = -arch ev6
endif


ifeq ($(HOSTNAME),wiseguy.physics.uiuc.edu)
CFLAGS2 = -arch ev56
endif

ifeq ($(HOSTNAME),alphadog.physics.uiuc.edu)
CFLAGS2 = -arch ev56
endif

endif # endif usempich==1

ifeq ($(USEMPI),0)

ifeq ($(HOSTNAME),photon.physics.uiuc.edu)
CFLAGS2 = -arch ev67
endif

ifeq ($(HOSTNAME),rainman.physics.uiuc.edu)
CFLAGS2 = -arch ev6
endif

ifeq ($(HOSTNAME),wiseguy.physics.uiuc.edu)
CFLAGS2 = -arch ev56
endif

ifeq ($(HOSTNAME),alphadog.physics.uiuc.edu)
CFLAGS2 = -arch ev56
endif

CC=$(COMP)  $(CFLAGS3) $(CFLAGS2) $(CDEBUG)

endif # endif usempich==0

endif # endif useccc==1

ifeq ($(USECCC),0)
#CDEBUG = -g
CDEBUG = -pg -g
#-pg
#-pg -g  source lines
#-pg -g -a   line count (doesn't work on photon)
# gprof -l <file> > out.txt
# gprof -A -I<sourcedir>
# gprof -l -A -x s
#CDEBUG =
CFLAGS3 = -O4 -Wall
#CFLAGS3 = -O0 -Wall


ifeq ($(USEMPI),1)

CC=$(MCC) -cc=$(COMP) $(CFLAGS3) $(CFLAGS2) $(CDEBUG)

endif

ifeq ($(USEMPI),0)

CC=$(COMP) $(CFLAGS3) $(CFLAGS2) $(CDEBUG)

endif


endif # endif use usempich==1



endif # endif alpha


#
# Define source files
#
SRCS = \
$(SRCD)/analsol.$(SUFF) \
$(SRCD)/boundrect1.$(SUFF) \
$(SRCD)/boundgen1.$(SUFF) \
$(SRCD)/boundgensimple1.$(SUFF) \
$(SRCD)/boundmpi.$(SUFF) \
$(SRCD)/diag.$(SUFF) \
$(SRCD)/diagflux.$(SUFF) \
$(SRCD)/gpar.$(SUFF) \
$(SRCD)/init.$(SUFF) \
$(SRCD)/postproc.$(SUFF) \
$(SRCD)/numerics.$(SUFF) \
$(SRCD)/ranc.$(SUFF) \
$(SRCD)/stepgen.$(SUFF) \
$(SRCD)/step2d.$(SUFF) \
$(SRCD)/stepmag.$(SUFF) \
$(SRCD)/sweep.$(SUFF) \
$(SRCD)/dqcalc.$(SUFF) \
$(SRCD)/timestep.$(SUFF) \
$(SRCD)/utilfun.$(SUFF) \
$(SRCD)/steptvdlf.$(SUFF)

#
# Define object files
#                                         
OBJS = \
$(OBJD)/analsol.o \
$(OBJD)/boundrect1.o \
$(OBJD)/boundgen1.o \
$(OBJD)/boundgensimple1.o \
$(OBJD)/boundmpi.o \
$(OBJD)/diag.o \
$(OBJD)/diagflux.o \
$(OBJD)/gpar.o \
$(OBJD)/init.o \
$(OBJD)/postproc.o \
$(OBJD)/numerics.o \
$(OBJD)/ranc.o \
$(OBJD)/stepgen.o \
$(OBJD)/step2d.o \
$(OBJD)/stepmag.o \
$(OBJD)/sweep.o \
$(OBJD)/dqcalc.o \
$(OBJD)/timestep.o \
$(OBJD)/utilfun.o \
$(OBJD)/steptvdlf.o





#
all: 	$(PREP) $(OBJD) $(BIND) $(CMD) $(FINISH)
#

$(PREP):
	( sh ./makedecs.h.sh )
#	( sh ./makenopp.sh ) # just do manually

$(OBJD):
	mkdir $(OBJD)

$(BIND):
	mkdir $(BIND)

$(CMD): $(OBJS) $(SRCS) makefile.pp
	$(CC) $(CFLAGS) -o $(CMD) $(OBJS) $(LDFLAGS)

$(FINISH):
	#( touch twodcode.tgz )
	#( rm twodcode.tgz )
	#tar cvzf $(OBJD)/twodcode.tgz *
#	mv -f ./*.o $(OBJD2)

# dependencies
$(OBJS)	: $(SRCD)/defs.h $(SRCD)/global.h $(SRCD)/global2dsup.h $(SRCD)/global3dsup.h $(SRCD)/makefile.pp

$(OBJD)/postproc.o : $(SRCD)/postproc.h
$(OBJD)/timestep.o : $(SRCD)/timestep.h
$(OBJD)/boundrect1.o : $(SRCD)/boundrect.h $(SRCD)/boundpos.h $(SRCD)/bound.h
$(OBJD)/boundgen1.o : $(SRCD)/boundgen.h $(SRCD)/bound.h
$(OBJD)/boundgensimple1.o : $(SRCD)/boundgen.h $(SRCD)/bound.h
$(OBJD)/boundmpi.o : $(SRCD)/bound.h
$(OBJD)/diag.o : $(SRCD)/imageloophead.h $(SRCD)/imageloopinside.h $(SRCD)/diag.h
$(OBJD)/stepgen.o : $(SRCD)/step.h
$(OBJD)/step2d.o : $(SRCD)/step.h
$(OBJD)/stepmag.o : $(SRCD)/step.h


cleandat:
	$(RM) $(BIND)/*.dat $(BIND)/*.dat.ras
 
clean:
	$(RM) *.o $(CMD)

