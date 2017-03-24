CC=g++
#CFLAGS=    -ggdb -fPIC -march=native -mtune=native -std=c++11 -O3 -fopenmp
CFLAGS= -fPIC -march=native -mtune=native -std=c++11 -O3 -fopenmp
DEBUGFLAGS=-ggdb -fPIC -march=native -mtune=native -std=c++11

all: main3
	./main3

main3: takagi_consize_refactor.cpp takagi.h vector3d.hpp main2.cpp
	$(CC) $(CFLAGS) main2.cpp takagi_consize_refactor.cpp -o main3
	#$(CC) $(DEBUGFLAGS) main2.cpp takagi_consize_refactor.cpp -o maind3


clean:
	rm takagi.o *.data main3 maind3
