DIR_INC := ./inc
DIR_SRC := ./
DIR_OBJ := ./obj

PREFIX ?= ./
BINDIR ?= $(PREFIX)/bin
INCLUDE_DIRS ?=
LIBRARY_DIRS ?=

SRC := $(wildcard ${DIR_SRC}/*.cpp)
OBJ := $(patsubst %.cpp,${DIR_OBJ}/%.o,$(notdir ${SRC}))

TARGET := ../inst/fastp

BIN_TARGET := ${TARGET}

CXX ?= g++
CXXFLAGS := -std=c++11 -g -O3 -I${DIR_INC} $(foreach includedir,$(INCLUDE_DIRS),-I$(includedir)) ${CXXFLAGS}
LIBS := -lz -lpthread
LD_FLAGS := $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir)) $(LIBS) $(LD_FLAGS)

FASTP_BIN_LIST = fastp

${BIN_TARGET}:${OBJ}
	$(CXX) $(OBJ) -o $@ $(LD_FLAGS)

${DIR_OBJ}/%.o:${DIR_SRC}/%.cpp make_obj_dir
	$(CXX) -c $< -o $@ $(CXXFLAGS)

.PHONY:clean
clean:
	rm obj/*.o
	rm $(TARGET)

make_obj_dir:
	@if test ! -d $(DIR_OBJ) ; \
	then \
		mkdir $(DIR_OBJ) ; \
	fi


