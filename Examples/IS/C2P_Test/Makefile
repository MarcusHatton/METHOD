# Make file for the main function. Builds all modules
# and links for main.cc

# Compiler
CC = mpic++
##CC = g++

USE_MPI = 0
USE_OMP = 0
USE_HDF = 1

PROJECT_TYPE = CPU

# Module directory
MODULE_DIR = ./../../../Project/${PROJECT_TYPE}/Src

# Rootfinder directory
RTFIND_DIR = ./../../../Project/${PROJECT_TYPE}/CminpackLibrary

# Include directory
INC_DIR = ./../../../Project/${PROJECT_TYPE}/Include

# Cminpack include directory
RTFIND_INC_DIR = ./../../../Project/${PROJECT_TYPE}/CminpackLibrary/Include

# Cminpack source directory
RTFIND_SRC_DIR = ./../../../Project/${PROJECT_TYPE}//CminpackLibrary/Src

ifeq ($(USE_OMP), 1)
    OMP_FLAGS = -fopenmp
endif

# C++ compiler flags
CXXFLAGS = -Wall -std=c++11 -O0 $(OMP_FLAGS) -Wno-unknown-pragmas -g

# Sources
SRCS = simulation.cc \
			 IS.cc \
			 initFunc.cc \
			 simData.cc \
			 boundaryConds.cc

SERIAL_SRCS = serialSaveData.cc \
			  serialEnv.cc

ifeq ($(USE_HDF), 1)
    SERIAL_SRCS += serialSaveDataHDF5.cc
endif

PARALLEL_SRCS = parallelSaveData.cc \
				parallelBoundaryConds.cc \
				parallelEnv.cc

ifeq ($(USE_HDF), 1)
    PARALLEL_SRCS += parallelSaveDataHDF5.cc
endif

ifeq ($(USE_MPI), 1)
    SRCS += ${PARALLEL_SRCS}
else
    SRCS += ${SERIAL_SRCS}
endif

ifeq ($(USE_HDF), 1)
    export HDF5_CXX := $(CC)
    export HDF5_CLINKER := $(CC)

    ifeq ($(USE_OMP), 1)
		CC = h5pcc.openmpi
    else
    	CC = h5pcc
    endif
endif

# Headers
HDRS = ${SRCS:.cc=.h} cudaErrorCheck.h

# Objects
OBJS = ${SRCS:.cc=.o}

# Rootfinder objects
RTFIND_OBJS = $(RTFIND_SRC_DIR)/dogleg.o \
							$(RTFIND_SRC_DIR)/dpmpar.o \
							$(RTFIND_SRC_DIR)/enorm.o \
							$(RTFIND_SRC_DIR)/fdjac1.o \
							$(RTFIND_SRC_DIR)/hybrd.o \
							$(RTFIND_SRC_DIR)/hybrd1.o \
							$(RTFIND_SRC_DIR)/qform.o \
							$(RTFIND_SRC_DIR)/qrfac.o \
							$(RTFIND_SRC_DIR)/r1mpyq.o \
							$(RTFIND_SRC_DIR)/r1updt.o

# Programmes
PROGS = ${SRCS:.cc=}

# Main programme executable
EXEC = main

# Rootfinder executables
RTFIND = buildRootfinder


# Build everything
build : $(RTFIND) $(EXEC)

# Make and run all tests on simulation scripts and run main
all : $(RTFIND) $(EXEC)
	cd $(TEST_DIR) && $(MAKE) test
	@$(MAKE) run

# Build and run main (does not run any tests if called explicitly)
run : $(RTFIND) $(EXEC)
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	@echo ''
	@echo ''
	@echo '######################################################'
	@echo '#              Executing main programme              #'
	@echo '######################################################'
	@./$(EXEC)

clean :
	rm	-f	$(EXEC)	*.o	*.gch


#################
# Build objects #
#################

simData.o : $(MODULE_DIR)/simData.cc $(INC_DIR)/simData.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

initFunc.o : $(MODULE_DIR)/initFunc.cc $(INC_DIR)/initFunc.h	$(INC_DIR)/simData.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

IS.o : $(MODULE_DIR)/IS.cc $(INC_DIR)/IS.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR) -I$(RTFIND_INC_DIR)

simulation.o : $(MODULE_DIR)/simulation.cc	$(INC_DIR)/simulation.h $(INC_DIR)/model.h	$(INC_DIR)/simData.h $(INC_DIR)/saveData.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

main.o : ./main.cc $(INC_DIR)/simulation.h $(INC_DIR)/model.h $(INC_DIR)/simData.h $(INC_DIR)/initFunc.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR) -I$(RTFIND_INC_DIR) $(CXXFLAGS)

boundaryConds.o : $(MODULE_DIR)/boundaryConds.cc $(INC_DIR)/boundaryConds.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

serialEnv.o : $(MODULE_DIR)/serialEnv.cc $(INC_DIR)/serialEnv.h
	$(CC)   $<  -c  $(CXXFLAGS) -I$(INC_DIR)

# Executable
main : main.o $(OBJS) $(RTFIND_OBJS)
	$(CC)	$^ -o $@ $(CXXFLAGS) $(CXXFLAGS)

buildRootfinder:
	@cd $(RTFIND_DIR) && $(MAKE)	objects
