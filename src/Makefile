DIR_INC_FASTP := ./inc
DIR_SRC_FASTP := src_fastp
DIR_OBJ_FASTP := obj_fastp
DIR_INC_GENCORE := ./inc
DIR_SRC_GENCORE := src_gencore
DIR_OBJ_GENCORE := obj_gencore

RHTSLIB_LIBS=$(shell echo 'cat(file.path(dirname(system.file("usrlib",package = "Rhtslib", mustWork = TRUE)),"include",""))'|\
    "$(R_HOME)/bin/R" --vanilla --slave)
RHTSLIB_CPPFLAGS=$(shell echo 'cat(file.path(system.file("usrlib",package = "Rhtslib", mustWork = TRUE),"libhts.a"))'|\
    "$(R_HOME)/bin/R" --vanilla --slave)


PKG_LIBS=$(RHTSLIB_LIBS)
PKG_CPPFLAGS=$(RHTSLIB_CPPFLAGS)

PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin
INCLUDE_DIRS ?= 
LIBRARY_DIRS ?= 

SRC_FASTP := $(wildcard ${DIR_SRC_FASTP}/*.cpp)
OBJ_FASTP := $(patsubst %.cpp,${DIR_OBJ_FASTP}/%.o,$(notdir ${SRC_FASTP}))
SRC_GENCORE = $(wildcard ${DIR_SRC_GENCORE}/*.cpp)  
OBJ_GENCORE = $(patsubst %.cpp,${DIR_OBJ_GENCORE}/%.o,$(notdir ${SRC_GENCORE})) 

TARGET_FASTP := fastp
TARGET_GENCORE = gencore

BIN_TARGET_FASTP := ${TARGET_FASTP}
FINAL_TARGET_FASTP := ../inst/${TARGET_FASTP}
BIN_TARGET_GENCORE = ${TARGET_GENCORE}
FINAL_TARGET_GENCORE := ../inst/${TARGET_GENCORE}


CXX = g++
CXXFLAGS := -std=c++11 -g -O3 -I${DIR_INC_FASTP} $(foreach includedir,$(INCLUDE_DIRS),-I$(includedir)) ${CXXFLAGS}
LIBS := -lz -lpthread
LD_FLAGS := $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir)) $(LIBS) $(LD_FLAGS)
CC = g++
CFLAGS = -std=c++11 -g -I${DIR_INC_GENCORE} $(foreach includedir,$(INCLUDE_DIRS),-I$(includedir))

${FINAL_TARGET_GENCORE}: ${BIN_TARGET_GENCORE}
	mv ${BIN_TARGET_GENCORE} ${FINAL_TARGET_GENCORE}

${BIN_TARGET_GENCORE}:${OBJ_GENCORE} ${FINAL_TARGET_FASTP}
	$(CXX) $(OBJ_GENCORE) -L. $(RHTSLIB_CPPFLAGS) -lz -lm -lbz2 -llzma -lcurl -lz  -lpthread -o $@
    
${DIR_OBJ_GENCORE}/%.o:${DIR_SRC_GENCORE}/%.cpp make_obj_dir2 ${FINAL_TARGET_FASTP}
	$(CXX) $(CFLAGS) -I$(RHTSLIB_LIBS) -O3 -c $< -o $@
	
${FINAL_TARGET_FASTP}: ${BIN_TARGET_FASTP}
	mv ${BIN_TARGET_FASTP} ${FINAL_TARGET_FASTP}
	
${BIN_TARGET_FASTP}:${OBJ_FASTP}
	$(CXX) $(OBJ_FASTP) -o $@ $(LD_FLAGS)

${DIR_OBJ_FASTP}/%.o:${DIR_SRC_FASTP}/%.cpp make_obj_dir
	$(CXX) -c $< -o $@ $(CXXFLAGS)

make_obj_dir:
	@if test ! -d $(DIR_OBJ_FASTP) ; \
	then \
		mkdir $(DIR_OBJ_FASTP) ; \
	fi

make_obj_dir2:
	@if test ! -d $(DIR_OBJ_GENCORE) ; \
	then \
		mkdir $(DIR_OBJ_GENCORE) ; \
	fi
