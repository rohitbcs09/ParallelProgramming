CC=icpc
CFLAGS=-O0 -std=c++11

all: 
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o steal_L1 mat_mul_steal_L1.cpp ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o steal_L2 mat_mul_steal_L2.cpp ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o steal_L3 mat_mul_steal_L3.cpp ${TACC_PAPI_LIB}/libpapi.a
clean:
	rm -rf steal*
