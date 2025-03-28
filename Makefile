# Makefile for your_project
# This Makefile builds the project executables, organizes the build process,
# and includes custom targets for cleaning and generating TAGS for code navigation.

# -----------------------------------------------------
# Directory Variables
# -----------------------------------------------------
SRCDIR    = src
INCDIR    = include
OBJDIR    = obj
BINDIR    = bin
SCRIPTDIR = scripts

# -----------------------------------------------------
# Compiler and Flags
# -----------------------------------------------------
CC        = mpicc
CFLAGS    = -Wall -g -I$(INCDIR)
# Update LDFLAGS with correct rpath and rpath-link flags based on module paths
LDFLAGS   = -Wl,-rpath,/sw/eb/sw/Hypre/2.28.0-foss-2022b \
            -Wl,-rpath,/sw/eb/sw/SuperLU_DIST/8.1.2-foss-2022b \
            -Wl,-rpath-link,/sw/eb/sw/Hypre/2.28.0-foss-2022b \
            -Wl,-rpath-link,/sw/eb/sw/SuperLU_DIST/8.1.2-foss-2022b
LIBS      =
CLINKER   = $(CC)
PETSC_LIB =

ifdef TEC360HOME
CFLAGS   += -I${TEC360HOME}/include/ -DTECIO=1
LIBS     += ${TEC360HOME}/lib/tecio64.a -lstdc++
endif

# -----------------------------------------------------
# PETSc Integration
# -----------------------------------------------------
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

# -----------------------------------------------------
# Source and Object Files
# -----------------------------------------------------
SOURCEC   = $(wildcard $(SRCDIR)/*.c)
OBJSC     = $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(SOURCEC))

# -----------------------------------------------------
# Executable Names
# -----------------------------------------------------
TESTT_EXE         = $(BINDIR)/testt
DATA_EXE          = $(BINDIR)/data
INTTEST_EXE       = $(BINDIR)/inttest
ITFCSEARCH_EXE    = $(BINDIR)/itfcsearch
DATA_VTK_EXE      = $(BINDIR)/data_vtk
DATALIS_EXE       = $(BINDIR)/datalis
DATAFILE_EXE      = $(BINDIR)/datafile
SWARM_TEST_EXE    = $(BINDIR)/swarm_test

# New Executable
POSTPROCESS_EXE   = $(BINDIR)/postprocess

# -----------------------------------------------------
# Phony Targets
# -----------------------------------------------------
.PHONY: all cleanobj clean_all tags \
        testt data inttest \
        itfcsearch data_vtk datalis datafile swarm_test \
        postprocess

# -----------------------------------------------------
# Default Target
# -----------------------------------------------------
all: testt data inttest \
     itfcsearch data_vtk datalis datafile swarm_test \
     postprocess

# -----------------------------------------------------
# Create Necessary Directories
# -----------------------------------------------------
dirs: | $(OBJDIR) $(BINDIR)

$(OBJDIR):
	@mkdir -p $(OBJDIR)

$(BINDIR):
	@mkdir -p $(BINDIR)

# -----------------------------------------------------
# Compilation Rule
# -----------------------------------------------------
$(OBJDIR)/%.o: $(SRCDIR)/%.c | $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@

# -----------------------------------------------------
# 1. testt executable
# -----------------------------------------------------
testt: dirs $(TESTT_EXE)

$(TESTT_EXE): $(OBJSC)
	$(CLINKER) $(CFLAGS) -o $@ $(OBJSC) $(PETSC_LIB) $(LDFLAGS) $(LIBS)

# -----------------------------------------------------
# 2. data executable
# -----------------------------------------------------
data: dirs $(DATA_EXE)

$(DATA_EXE): $(OBJDIR)/variables.o $(OBJDIR)/compgeom.o $(OBJDIR)/data_ibm.o \
             $(OBJDIR)/ibm_io.o $(OBJDIR)/fsi.o $(OBJDIR)/fsi_move.o \
             $(OBJDIR)/fish.o $(OBJDIR)/data.o
	$(CLINKER) $(CFLAGS) -o $@ $^ $(PETSC_LIB) $(PETSC_SNES_LIB) $(PETSC_TS_LIB) $(LDFLAGS) $(LIBS)

# -----------------------------------------------------
# 3. inttest (swarm_interp) executable
# -----------------------------------------------------
inttest: dirs $(INTTEST_EXE)

$(INTTEST_EXE): $(OBJDIR)/inttest.o $(OBJDIR)/interpolation.o $(OBJDIR)/walkingsearch.o \
                $(OBJDIR)/ParticleSwarm.o $(OBJDIR)/logging.o $(OBJDIR)/setup.o \
                $(OBJDIR)/AnalyticalSolution.o $(OBJDIR)/grid.o $(OBJDIR)/io.o $(OBJDIR)/ParticleMotion.o
	$(CLINKER) $(CFLAGS) -o $@ $^ $(PETSC_LIB) $(LDFLAGS) $(LIBS)

# -----------------------------------------------------
# 4. itfcsearch
# -----------------------------------------------------
itfcsearch: dirs $(ITFCSEARCH_EXE)

$(ITFCSEARCH_EXE): $(OBJDIR)/itfcsearch.o \
                   $(OBJDIR)/walkingsearch.o $(OBJDIR)/logging.o \
                   $(OBJDIR)/grid.o $(OBJDIR)/io.o
	$(CLINKER) $(CFLAGS) -o $@ $^ $(PETSC_LIB) $(LDFLAGS) $(LIBS)

# -----------------------------------------------------
# 5. data_vtk
# -----------------------------------------------------
data_vtk: dirs $(DATA_VTK_EXE)

$(DATA_VTK_EXE): $(OBJDIR)/data_vtk.o $(OBJDIR)/logging.o \
                 $(OBJDIR)/io.o $(OBJDIR)/grid.o
	$(CLINKER) $(CFLAGS) -o $@ $^ $(PETSC_LIB) $(LDFLAGS) $(LIBS)

# -----------------------------------------------------
# 6. datalis
# -----------------------------------------------------
datalis: dirs $(DATALIS_EXE)

$(DATALIS_EXE): $(OBJDIR)/datalis.o $(OBJDIR)/logging.o $(OBJDIR)/io.o
	$(CLINKER) $(CFLAGS) -o $@ $^ $(PETSC_LIB) $(LDFLAGS) $(LIBS)

# -----------------------------------------------------
# 7. datafile
# -----------------------------------------------------
datafile: dirs $(DATAFILE_EXE)

$(DATAFILE_EXE): $(OBJDIR)/datafile.o $(OBJDIR)/logging.o $(OBJDIR)/io.o
	$(CLINKER) $(CFLAGS) -o $@ $^ $(PETSC_LIB) $(LDFLAGS) $(LIBS)

# -----------------------------------------------------
# 8. swarm_test
# -----------------------------------------------------
swarm_test: dirs $(SWARM_TEST_EXE)

$(SWARM_TEST_EXE): $(OBJDIR)/swarm_test.o $(OBJDIR)/logging.o \
                   $(OBJDIR)/io.o $(OBJDIR)/grid.o
	$(CLINKER) $(CFLAGS) -o $@ $^ $(PETSC_LIB) $(LDFLAGS) $(LIBS)

# -----------------------------------------------------
# NEW: 9. postprocess executable
# -----------------------------------------------------
postprocess: dirs $(POSTPROCESS_EXE)

$(POSTPROCESS_EXE): $(OBJDIR)/postprocess.o $(OBJDIR)/interpolation.o \
			$(OBJDIR)/walkingsearch.o $(OBJDIR)/grid.o $(OBJDIR)/ParticleSwarm.o \
			$(OBJDIR)/logging.o $(OBJDIR)/io.o $(OBJDIR)/setup.o 
	$(CLINKER) $(CFLAGS) -o $@ $^ $(PETSC_LIB) $(LDFLAGS) $(LIBS)

# -----------------------------------------------------
# Custom Targets
# -----------------------------------------------------
.PHONY: tags
tags:
	find $(SRCDIR) $(INCDIR) $(SCRIPTDIR) -type f \( -name "*.c" -o -name "*.h" -o -name "*.py" \) -print | etags -f TAGS -

.PHONY: clean_tags
clean_tags:
	rm -f TAGS

cleanobj:
	rm -f $(OBJDIR)/*.o

clean_all: cleanobj
	rm -f $(BINDIR)/* TAGS

# Include PETSc test rules
include ${PETSC_DIR}/lib/petsc/conf/test
