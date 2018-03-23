CC=icpc
CFLAGS=-O0

all: 
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_ijk iter_mm_ijk.cpp ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_ikj iter_mm_ikj.cpp ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_jik iter_mm_jik.cpp ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_jki iter_mm_jki.cpp ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_kij iter_mm_kij.cpp ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_kji iter_mm_kji.cpp ${TACC_PAPI_LIB}/libpapi.a
		
clean:
	rm -rf a_ijk a_ijk_output a_ikj a_ikj_output a_jik a_jik_output  a_jki a_jki_output a_kij a_kij_output a_kji a_kji_output *.o
