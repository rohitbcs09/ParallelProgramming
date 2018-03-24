CC=icpc
CFLAGS=-O0

all: 
	$(CC) $(CFLAGS) -o a_ijk iter_mm_ijk.cpp
	$(CC) $(CFLAGS) -o a_ikj iter_mm_ikj.cpp 
	$(CC) $(CFLAGS) -o a_jik iter_mm_jik.cpp
	$(CC) $(CFLAGS) -o a_jki iter_mm_jki.cpp 
	$(CC) $(CFLAGS) -o a_kij iter_mm_kij.cpp
	$(CC) $(CFLAGS) -o a_kji iter_mm_kji.cpp 
		
clean:
	rm -rf a_ijk a_ijk_output a_ikj a_ikj_output a_jik a_jik_output  a_jki a_jki_output a_kij a_kij_output a_kji a_kji_output *.o
