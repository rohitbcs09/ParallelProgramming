CC=icpc
CFLAGS=-O0

all: 
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_ikj_i iter_mm_ikj_i.cpp -std=c++11 ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_ikj_j iter_mm_ikj_j.cpp -std=c++11 ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_ikj_k iter_mm_ikj_k.cpp -std=c++11 ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_ikj_ij iter_mm_ikj_ij.cpp -std=c++11 ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_ikj_jk iter_mm_ikj_jk.cpp -std=c++11 ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_ikj_ik iter_mm_ikj_ik.cpp -std=c++11 ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_ikj_ijk iter_mm_ikj_ijk.cpp -std=c++11 ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_kij_i iter_mm_kij_i.cpp -std=c++11 ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_kij_j iter_mm_kij_j.cpp -std=c++11 ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_kij_k iter_mm_kij_k.cpp -std=c++11 ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_kij_ij iter_mm_kij_ij.cpp -std=c++11 ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_kij_jk iter_mm_kij_jk.cpp -std=c++11 ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_kij_ik iter_mm_kij_ik.cpp -std=c++11 ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_kij_ijk iter_mm_kij_ijk.cpp -std=c++11 ${TACC_PAPI_LIB}/libpapi.a
		
clean:
	rm -rf a_ikj*  a_kij* output* *.o
