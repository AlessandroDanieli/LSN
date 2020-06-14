CC = g++#mpicxx
CFLAGS = -Wall -O3 --std=c++11
LIBS = main.o random.o simulannealing.o vector.o

main.exe: ${LIBS}
	$(CC) -o $@ $^ ${CFLAGS}
%.o: %.cpp %.h
	$(CC) -c $< ${CFLAGS}
	
clean:
	rm -rf *.out *.o
	rm -rf PathLengths/*.out
final:
	rm -rf *.o
all:
	@make -s
	@./main.exe
