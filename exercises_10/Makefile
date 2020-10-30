CC = g++
CFLAGS = -Wall -O3 --std=c++11

simulatedannealing.exe : simulatedannealing.o random.o
	$(CC) random.o simulatedannealing.o -o simulatedannealing.exe
simulatedannealing.o : simulatedannealing.cpp
	$(CC) -c simulatedannealing.cpp -o simulatedannealing.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o simulatedannealing.exe
