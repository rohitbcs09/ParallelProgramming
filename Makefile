CXX=g++

CXXFLAGS=-std=c++11 -Iinclude

TARGET=src/rec_mat_mul.cpp

BINARY=out
out=src/rec_mat_mul.o
all: $(TARGET)

$(TARGET): $(TARGET).cpp
	$(CXX) $(CXXFLAGS) -o $(BINARY) $(TARGET)

clean: $(RM) $(BINARY)
