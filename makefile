# Compiler
CXX = mpic++
# Compiler flags
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra

# HDF5 include and library paths
HDF5_INC = /usr/local/hdf5/include
HDF5_LIB = /usr/local/hdf5/lib
HDF5_LIBS = -lhdf5 -lhdf5_cpp

# FFTW library (add if used)
FFTW_LIBS = -lfftw3 -lm

# Output binary name
TARGET = solver

# Source file
SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)

# Build target
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -I$(HDF5_INC) -L$(HDF5_LIB) $(OBJS) -o $(TARGET) $(HDF5_LIBS) $(FFTW_LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -I$(HDF5_INC) -c $< -o $@

# Clean target
clean:
	rm -f $(TARGET) $(OBJS)