FFLAGS= -fpp -mkl -openmp -O3 -mtune=native -march=native -vec-report
CCFLAGS = -fopenmp  -mcmodel=medium -Wall -pedantic -Ofast -funroll-loops -finline-functions -ftree-vectorize  

LDFLAGS = -lstdc++ -lgomp
LDFLAGS_SEQ = -lstdc++
LINKS=/Soft/intelf/composer_xe_2013.3.163/compiler/lib/intel64/libiomp5.a $(LDFLAGS)
LINKS_SEQ= $(LDFLAGS)

FC=/Soft/intelf/composer_xe_2013.3.163/bin/intel64/ifort
CXX=/Soft/gcc/6.2.0/bin/g++

KIND_OBJS=  
F90_LEVELS_NEIG_OBJS= quicksort.o random.o glb.o acorni.o beyond.o chknam.o cova_1D.o cova3.o exhaustive_srch.o gauinv.o getindx.o grid_info.o aniso_d.o graph_subs.o kdtree3.o ksol.o locate.o ordrel.o picksupr.o powint.o setrot.o setsupr.o sortem.o sqdist.o sisim_lva_levels_neig.o Scaling-ISOMAP_landmark.o 
FOR_OBJS= 
LAPACK_OBJS=LAPACK/eig.o LAPACK/solve.o
OBJS=   eig.o solve.o
CPP_OBJS= boostdijkstra.o boostdijkstra_capi.o Boostdijkstra.o
ALL_LEVELS_NEIG_OBJS= $(CPP_OBJS) $(FOR_OBJS) $(F90_LEVELS_NEIG_OBJS)

OUT_LEVELS_NEIG=SISIM_LVA_levels_neig

default: $(KIND_OBJS) $(LAPACK_OBJS) $(ALL_LEVELS_NEIG_OBJS)
	$(FC) $(FFLAGS) $(LINKS) $(KIND_OBJS) $(OBJS) $(ALL_LEVELS_NEIG_OBJS) -o $(OUT_LEVELS_NEIG)
	cp $(OUT_LEVELS_NEIG) ../bin/

clean:
	rm *.mod *.o $(OUT_LEVELS_NEIG) ../bin/$(OUT_LEVELS_NEIG)

.SUFFIXES: .o .F90 .f90 .for .f .fpp .cpp

.F90.o : ; $(FC) -c $(FFLAGS) $*.F90
.f90.o : ; $(FC) -c $(FFLAGS) $*.f90
.for.o : ; $(FC) -c $(FFLAGS) $*.for
.f.o : ; $(FC) -c $(FFLAGS) $*.f
.fpp.o : ; $(FC) -c $(FFLAGS) $*.fpp
.cpp.o : ; ${CXX} -c ${CCFLAGS} $*.cpp
