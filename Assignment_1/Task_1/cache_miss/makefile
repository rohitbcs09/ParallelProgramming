CC=icpc
CFLAGS=-O0

all: 
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_ijk_L1 iter_mm_ijk_L1.cpp ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_ijk_L2 iter_mm_ijk_L2.cpp ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_ijk_L3 iter_mm_ijk_L3.cpp ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_ikj_L1 iter_mm_ikj_L1.cpp ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_ikj_L2 iter_mm_ikj_L2.cpp ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_ikj_L3 iter_mm_ikj_L3.cpp ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_jik_L1 iter_mm_jik_L1.cpp ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_jik_L2 iter_mm_jik_L2.cpp ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_jik_L3 iter_mm_jik_L3.cpp ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_jki_L1 iter_mm_jki_L1.cpp ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_jki_L2 iter_mm_jki_L2.cpp ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_jki_L3 iter_mm_jki_L3.cpp ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_kij_L1 iter_mm_kij_L1.cpp ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_kij_L2 iter_mm_kij_L2.cpp ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_kij_L3 iter_mm_kij_L3.cpp ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_kji_L1 iter_mm_kji_L1.cpp ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_kji_L2 iter_mm_kji_L2.cpp ${TACC_PAPI_LIB}/libpapi.a
	$(CC) -I${TACC_PAPI_INC} $(CFLAGS) -o a_kji_L3 iter_mm_kji_L3.cpp ${TACC_PAPI_LIB}/libpapi.a

clean:
	rm -rf a* *.o output*
