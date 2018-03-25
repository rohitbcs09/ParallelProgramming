CXX=g++

CXXFLAGS=-std=c++11 

TARGET=mat_mul_sharing.cpp

INCLUDE=mat_mul_sharing.h

BINARY=mat_mul_sharing

mat_mul_sharing:mat_mul_sharing.o

all: $(TARGET)

$(BINARY): $(BINARY).o
	@$(CXX) $(CXXFLAGS) -o $(BINARY) $(BINARY).o

$(BINARY).o: $(TARGET) $(INCLUDE)
	@$(CXX) $(CXXFLAGS) -c $(TARGET)

clean: 
	@rm *.o mat_mul_sharing
