CC = g++
CFLAGS = -Wall -O3 --std=c++11

genetic.exe : genetic.o random.o
	$(CC) random.o genetic.o -o genetic.exe
genetic.o : genetic.cpp
	$(CC) -c genetic.cpp -o genetic.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o genetic.exe seed.out
