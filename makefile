 ALL:

ifdef TEC360HOME
CFLAGS		 = -I${TEC360HOME}/include/ -DTECIO=1
LIBS		 = ${TEC360HOME}/lib/tecio64.a -lstdc++
else
CFLAGS		 = -Wall -g 
LIBS             = 
endif
LDFLAGS		 = 
FFLAGS		 =
CPPFLAGS 	 =  
FPPFLAGS         =
LOCDIR		 = 
MANSEC           = SNES
LIBFLAG          =

SOURCEC = bcs.c bmv.c compgeom.c ibm.c ibm_io.c init.c\
          main.c metrics.c poisson.c rhs.c rheology.c\
	  variables.c fsi.c implicitsolver.c\
	  fsi_move.c solvers.c copepod.c fish.c cstart.c spline.c\
          les.c k-omega.c wallfunction.c rhs2.c platlet.c interpolation.c  poisson_hypre.c

OBJSC =  bcs.o bmv.o compgeom.o ibm.o ibm_io.o init.o\
         main.o metrics.o poisson.o rhs.o rheology.o\
         variables.o fsi.o implicitsolver.o\
         fsi_move.o solvers.o copepod.o fish.o cstart.o spline.o\
         les.o k-omega.o wallfunction.o rhs2.o  platlet.o interpolation.o poisson_hypre.o 

LIBBASE = libpetscmat

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

testt: ${OBJSC}
	-$(CLINKER) -o testt ${OBJSC} ${PETSC_LIB}

	rm main.o

data: variables.o  compgeom.o data_ibm.o ibm_io.o fsi.o fsi_move.o fish.o  data.o
	-${CLINKER} -o data variables.o compgeom.o data_ibm.o ibm_io.o fsi.o fsi_move.o fish.o data.o ${PETSC_LIB} ${PETSC_SNES_LIB} ${PETSC_TS_LIB} ${LIBS}


interpolation: interpolation.o
	  -${CLINKER} ${CFLAGS} -o interpolation interpolation.o ${PETSC_LIB}

swarm_interp: swarm_interp.o variables.o walkingsearch.o ParticleSwarm.o logging.o grid.o
	  -${CLINKER} ${CFLAGS} -o swarm_interp swarm_interp.o variables.o walkingsearch.o ParticleSwarm.o logging.o grid.o ${PETSC_LIB}


itfcsearch: itfcsearch.o variables.o compgeom.o
	-${CLINKER} -o itfcsearch itfcsearch.o  variables.o compgeom.o ${PETSC_SNES_LIB} ${PETSC_TS_LIB}

data_vtk: data_surface.o 
	-${CLINKER} -o data_vtk data_surface.o  ${PETSC_SNES_LIB} ${PETSC_TS_LIB}

datalis: data_file2lis.o
	-${CLINKER} -o datalis data_file2lis.o ${PETSC_SNES_LIB} ${PETSC_TS_LIB} ${LIBS}

datafile: data_list2file.o
	-${CLINKER} -o datafile data_list2file.o ${PETSC_SNES_LIB} ${PETSC_TS_LIB} ${LIBS}
cleanobj:
	rm -f *.o

# Swarm Test Executable
swarm_test: swarm_test.o variables.o walkingsearch.o
	-${CLINKER} -o swarm_test swarm_test.o variables.o walkingsearch.o ${PETSC_LIB} ${LIBS}

# Compilation rule for swarm_test.c
swarm_test.o: swarm_test.c
	${CLINKER} ${CFLAGS} -c swarm_test.c -o swarm_test.o ${PETSC_INCLUDE}

include ${PETSC_DIR}/lib/petsc/conf/test
