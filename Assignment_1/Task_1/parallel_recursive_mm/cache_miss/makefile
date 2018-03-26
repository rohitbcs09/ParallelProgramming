CC=icpc
CFLAGS=-O0

all: 
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_par_rec_mm_L1 par_rec_mm_L1.cpp ${TACC_PAPI_LIB}/libpapi.a -std=c++11
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_par_rec_mm_L2 par_rec_mm_L2.cpp ${TACC_PAPI_LIB}/libpapi.a -std=c++11
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_par_rec_mm_L3 par_rec_mm_L3.cpp ${TACC_PAPI_LIB}/libpapi.a -std=c++11

clean:
	rm -rf a* *.o output*
