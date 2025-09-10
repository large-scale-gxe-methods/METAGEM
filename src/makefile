CXX=g++
CXXFLAGS = -std=c++11 -O3 -g -Wall

# Uncomment line below for LAPACK/BLAS only and comment other lines beginning with LIBS.
#LIBS= -llapack -lblas -ldl -lm 

LIBS= -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -ldl -lm 

CXX_SOURCES= matrix_utils.cpp read_parameters.cpp file.cpp metagem.cpp
CXX_OBJECTS= matrix_utils.o read_parameters.o file.o metagem.o
OBJS= $(CXX_OBJECTS)

#targets
all: METAGEM

METAGEM: $(CXX_OBJECTS) $(C_OBJECTS)
	$(CXX) $^ $(CXXFLAGS) -o $@ $(LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f $(OBJS)
