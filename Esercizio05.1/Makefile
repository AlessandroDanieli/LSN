CC = g++
CFLAGS = -Wall -O3 --std=c++17
LIBS = main.o random.o stat.o vector.o orbital.o metropolisorbitals.o auxiliary.o

main.exe: ${LIBS}
	$(CC) -o $@ $^ ${CFLAGS}
%.o: %.cpp %.h
	$(CC) -c $< ${CFLAGS}
clean:
	rm -rf *.o main.exe *.out *.xyz
	rm -rf OutputData/*.out OutputData/*.xyz
all:
	@make -s
	@./main.exe
