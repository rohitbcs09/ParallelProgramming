CXX=g++

CXXFLAGS=-std=c++11 -Iinclude

TARGET=src/rec_mat_mul.cpp

INCLUDE=include/rec_mat_mul.h

BINARY=rec_mat_mul

rec_mat_mul:rec_mat_mul.o

all: $(TARGET)

$(BINARY): $(BINARY).o
	@$(CXX) $(CXXFLAGS) -o $(BINARY) $(BINARY).o

$(BINARY).o: $(TARGET) $(INCLUDE)
	@$(CXX) $(CXXFLAGS) -c $(TARGET)

clean: 
	@rm *.o rec_mat_mul
