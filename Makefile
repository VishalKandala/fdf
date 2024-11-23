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
LDFLAGS   =
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

# -----------------------------------------------------
# Phony Targets
# -----------------------------------------------------
.PHONY: all cleanobj clean_all tags \
        testt data inttest \
        itfcsearch data_vtk datalis datafile swarm_test

# -----------------------------------------------------
# Default Target
# -----------------------------------------------------
all: testt data inttest \
     itfcsearch data_vtk datalis datafile swarm_test

# -----------------------------------------------------
# Create Necessary Directories
# -----------------------------------------------------
# Ensure directories exist using order-only prerequisites
dirs: | $(OBJDIR) $(BINDIR)

$(OBJDIR):
	@mkdir -p $(OBJDIR)

$(BINDIR):
	@mkdir -p $(BINDIR)

# -----------------------------------------------------
# Compilation Rule
# -----------------------------------------------------
# Compile .c files in src/ to .o files in obj/
# $< refers to the source file
# $@ refers to the target object file
$(OBJDIR)/%.o: $(SRCDIR)/%.c | $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@

# -----------------------------------------------------
# Executable Targets
# -----------------------------------------------------
# 1. testt executable
testt: dirs $(TESTT_EXE)

$(TESTT_EXE): $(OBJSC)
	$(CLINKER) $(CFLAGS) -o $@ $(OBJSC) $(PETSC_LIB)

# 2. data executable
data: dirs $(DATA_EXE)

$(DATA_EXE): $(OBJDIR)/variables.o $(OBJDIR)/compgeom.o $(OBJDIR)/data_ibm.o \
             $(OBJDIR)/ibm_io.o $(OBJDIR)/fsi.o $(OBJDIR)/fsi_move.o \
             $(OBJDIR)/fish.o $(OBJDIR)/data.o
	$(CLINKER) -o $@ $^ $(PETSC_LIB) $(PETSC_SNES_LIB) $(PETSC_TS_LIB) $(LIBS)

# 4. swarm_interp executable
inttest: dirs $(INTTEST_EXE)

$(INTTEST_EXE): $(OBJDIR)/inttest.o $(OBJDIR)/interpolation.o $(OBJDIR)/walkingsearch.o \
                     $(OBJDIR)/ParticleSwarm.o $(OBJDIR)/logging.o \
                     $(OBJDIR)/grid.o  $(OBJDIR)/io.o 
	$(CLINKER) $(CFLAGS) -o $@ $^ $(PETSC_LIB)

# ... (Other executables defined similarly)

# -----------------------------------------------------
# Custom Targets
# -----------------------------------------------------
# Target for generating TAGS
.PHONY: tags
tags:
	find $(SRCDIR) $(INCLUDEDIR) $(SCRIPTSDIR) -type f \( -name "*.c" -o -name "*.h" -o -name "*.py" \) -print | etags -f TAGS -

# Clean TAGS file
.PHONY: clean_tags
clean_tags:
	rm -f TAGS

cleanobj:
	rm -f $(OBJDIR)/*.o

clean_all: cleanobj
	rm -f $(BINDIR)/* TAGS

# Include PETSc test rules
include ${PETSC_DIR}/lib/petsc/conf/test
