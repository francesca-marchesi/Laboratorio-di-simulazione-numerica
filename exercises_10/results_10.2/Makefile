CC = gcc
MPI = mpic++
CFLAGS = -Wall -O3 --std=c++11

MPI.exe : MPI.o random.o
	mpic++ MPI.o random.o -o MPI.exe
MPI.o : MPI.cpp
	mpic++ -c MPI.cpp -o MPI.o $(CFLAGS)
random.o : random.cpp random.h
	gcc -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o *.exe
