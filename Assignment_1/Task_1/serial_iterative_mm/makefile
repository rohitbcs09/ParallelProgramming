CC=icpc
CFLAGS=-O0

all: 
	$(CC) $(CFLAGS) -o a_ijk iter_mm_ijk.cpp -std=c++11
	$(CC) $(CFLAGS) -o a_ikj iter_mm_ikj.cpp -std=c++11
	$(CC) $(CFLAGS) -o a_jik iter_mm_jik.cpp -std=c++11
	$(CC) $(CFLAGS) -o a_jki iter_mm_jki.cpp -std=c++11
	$(CC) $(CFLAGS) -o a_kij iter_mm_kij.cpp -std=c++11
	$(CC) $(CFLAGS) -o a_kji iter_mm_kji.cpp -std=c++11
		
clean:
	rm -rf a_* output *.o
