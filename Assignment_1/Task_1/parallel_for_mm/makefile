CC=icpc
CFLAGS=-O0

all: 
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_ikj_i iter_mm_ikj_i.cpp -std=c++11 ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_ikj_j iter_mm_ikj_j.cpp -std=c++11 ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_ikj_ij iter_mm_ikj_ij.cpp -std=c++11 ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_kij_i iter_mm_kij_i.cpp -std=c++11 ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_kij_j iter_mm_kij_j.cpp -std=c++11 ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_kij_ij iter_mm_kij_ij.cpp -std=c++11 ${TACC_PAPI_LIB}/libpapi.a
		
clean:
	rm -rf a_ikj*  a_kij* output* *.o
