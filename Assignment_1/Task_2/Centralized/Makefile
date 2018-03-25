CXX=g++

CXXFLAGS=-std=c++11 

TARGET=mat_mul_central.cpp

INCLUDE=mat_mul_central.h

BINARY=mat_mul_central

mat_mul_central:mat_mul_central.o

all: $(TARGET)

$(BINARY): $(BINARY).o
	@$(CXX) $(CXXFLAGS) -o $(BINARY) $(BINARY).o

$(BINARY).o: $(TARGET) $(INCLUDE)
	@$(CXX) $(CXXFLAGS) -c $(TARGET)

clean: 
	@rm *.o mat_mul_central
